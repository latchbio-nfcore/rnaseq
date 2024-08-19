import csv
import os
import shutil
import subprocess
import sys
import time
import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Annotated, Any, List, Optional

import requests
from flytekit.core.annotation import FlyteAnnotation
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.utils import urljoins

sys.stdout.reconfigure(line_buffering=True)


@dataclass
class SampleSheet:
    """
    Represents a sample in the RNA-seq analysis.

    Attributes:
        sample (str): The name or identifier of the sample.
        fastq_1 (LatchFile): The first FASTQ file for the sample.
        fastq_2 (Optional[LatchFile]): The second FASTQ file for paired-end data (optional).
        strandedness (str): The strandedness of the library preparation.
    """

    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    strandedness: Optional[str] = None


class Reference_Type(Enum):
    """
    Enumeration of supported reference genomes.

    Each enum value represents a different species and its corresponding reference genome.
    """

    homo_sapiens = "Homo sapiens (RefSeq GRCh38.p14)"
    mus_musculus = "Mus musculus (RefSeq GRCm39)"
    rattus_norvegicus = "Rattus norvegicus (RefSeq GRCr8)"
    # drosophila_melanogaster = "Drosophila melanogaster (RefSeq Release_6_plus_ISO1_MT)"
    # rhesus_macaque = "Macaca mulatta (RefSeq rheMac10/Mmul_10)"
    saccharomyces_cerevisiae = "Saccharomyces cerevisiae (RefSeq R64)"


class Trimmer(Enum):
    """
    Enumeration of supported trimming tools.

    Attributes:
        trimgalore: TrimGalore trimming tool.
        fastp: fastp trimming tool.
    """

    trimgalore = "trimgalore"
    fastp = "fastp"


class UMIToolsGrouping(Enum):
    """
    Enumeration of UMI-tools grouping methods.

    These methods are used for deduplicating reads based on UMIs.
    """

    directional = "directional"
    unique = "unique"
    cluster = "cluster"
    percentile = "percentile"
    adjacency = "adjacency"


class Aligner(Enum):
    """
    Enumeration of supported alignment tools.

    Attributes:
        star_salmon: STAR aligner with Salmon quantification.
        star_rsem: STAR aligner with RSEM quantification.
        hisat2: HISAT2 aligner.
    """

    star_salmon = "star_salmon"
    star_rsem = "star_rsem"
    hisat2 = "hisat2"


class SalmonQuantLibType(Enum):
    """
    Enumeration of Salmon quantification library types.

    These specify the type of library preparation for accurate quantification.
    """

    A = "A"
    IS = "IS"
    ISF = "ISF"
    ISR = "ISR"
    IU = "IU"
    MS = "MS"
    MSF = "MSF"
    MSR = "MSR"
    MU = "MU"
    OS = "OS"
    OSF = "OSF"
    OSR = "OSR"
    OU = "OU"
    SF = "SF"
    SR = "SR"
    U = "U"


class PseudoAligner(Enum):
    """
    Enumeration of supported pseudo-alignment tools.

    Attributes:
        salmon: Salmon pseudo-aligner.
        kallisto: Kallisto pseudo-aligner.
    """

    salmon = "salmon"
    kallisto = "kallisto"


def get_flag_defaults(name: str, val: Any, default_val: Optional[Any]):
    """
    Generate command-line flags for Nextflow based on parameter values.

    Args:
        name (str): The name of the parameter.
        val (Any): The current value of the parameter.
        default_val (Optional[Any]): The default value of the parameter.

    Returns:
        str: A string containing the appropriate command-line flag, or an empty string if the value is default or None.
    """
    if val == default_val or val is None:
        return ""
    else:
        return get_flag(name=name, val=val)


def custom_samplesheet_constructor(
    samples: List[SampleSheet], shared_dir: Path
) -> Path:
    """
    Construct a custom sample sheet CSV file from the provided samples.

    This function creates a CSV file containing information about each sample,
    including sample name, FASTQ file paths, and strandedness. It also handles
    compression of FASTQ files if they are not already gzipped.

    Args:
        samples (List[SampleSheet]): A list of SampleSheet objects containing sample information.
        shared_dir (Path): The shared directory path for storing compressed files.

    Returns:
        Path: The path to the created sample sheet CSV file.
    """
    samplesheet = Path("/root/samplesheet.csv")

    columns = ["sample", "fastq_1", "fastq_2", "strandedness"]

    with open(samplesheet, "w") as f:
        writer = csv.DictWriter(f, columns, delimiter=",")
        writer.writeheader()

        for sample in samples:
            # Check and compress fastq_1 if needed
            fastq_1_path = sample.fastq_1.remote_path
            if not sample.fastq_1.remote_path.endswith(".gz"):
                local_path = Path(sample.fastq_1.local_path)
                compressed_path = shared_dir / f"{local_path.name}.gz"
                print(f"Compressing to {compressed_path}")
                subprocess.run(
                    ["pigz", "-p", "8", "-c", local_path],
                    stdout=open(compressed_path, "wb"),
                    check=True,
                )
                fastq_1_path = compressed_path

            # Check and compress fastq_2 if it exists and needs compression
            if sample.fastq_2:
                fastq_2_path = sample.fastq_2.remote_path
                if not sample.fastq_2.remote_path.endswith(".gz"):
                    local_path = Path(sample.fastq_2.local_path)
                    compressed_path = shared_dir / f"{local_path.name}.gz"
                    print(f"Compressing to {compressed_path}")
                    subprocess.run(
                        ["pigz", "-p", "8", "-c", local_path],
                        stdout=open(compressed_path, "wb"),
                        check=True,
                    )
                    fastq_2_path = compressed_path

            row_data = {
                "sample": sample.sample,
                "fastq_1": fastq_1_path,
                "fastq_2": fastq_2_path,
                "strandedness": sample.strandedness
                if sample.strandedness is not None
                else "auto",
            }
            writer.writerow(row_data)

    return samplesheet


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    """
    Initialize the workflow by provisioning a shared storage volume.

    This function requests a shared storage volume from the Nextflow dispatcher service
    and returns the name of the provisioned volume.

    Returns:
        str: The name of the provisioned storage volume.

    Raises:
        RuntimeError: If the execution token is not available.
    """
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_expiration_hours": 1,
            "storage_gib": 100,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@nextflow_runtime_task(cpu=8, memory=32, storage_gib=2000)
def nextflow_runtime(
    pvc_name: str,
    input: typing.List[SampleSheet],
    run_name: Annotated[
        str,
        FlyteAnnotation(
            {
                "rules": [
                    {
                        "regex": r"^[a-zA-Z0-9_-]+$",
                        "message": "ID name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
                    }
                ],
            }
        ),
    ],
    outdir: LatchOutputDir,
    genome_source: str,
    fasta: typing.Optional[LatchFile],
    gtf: typing.Optional[LatchFile],
    gff: typing.Optional[LatchFile],
    gene_bed: typing.Optional[LatchFile],
    transcript_fasta: typing.Optional[LatchFile],
    additional_fasta: typing.Optional[LatchFile],
    splicesites: typing.Optional[LatchFile],
    star_index: typing.Optional[LatchFile],
    hisat2_index: typing.Optional[LatchFile],
    rsem_index: typing.Optional[LatchFile],
    salmon_index: typing.Optional[LatchFile],
    kallisto_index: typing.Optional[LatchFile],
    extra_trimgalore_args: typing.Optional[str],
    extra_fastp_args: typing.Optional[str],
    bbsplit_fasta_list: typing.Optional[LatchFile],
    bbsplit_index: typing.Optional[LatchFile],
    ribo_database_manifest: typing.Optional[LatchFile],
    umitools_bc_pattern: typing.Optional[str],
    umitools_bc_pattern2: typing.Optional[str],
    umi_discard_read: typing.Optional[int],
    umitools_umi_separator: typing.Optional[str],
    pseudo_aligner: typing.Optional[PseudoAligner],
    salmon_quant_libtype: typing.Optional[SalmonQuantLibType],
    seq_center: typing.Optional[str],
    extra_star_align_args: typing.Optional[str],
    extra_salmon_quant_args: typing.Optional[str],
    extra_kallisto_quant_args: typing.Optional[str],
    email: typing.Optional[str],
    multiqc_title: typing.Optional[str],
    latch_genome: Reference_Type,
    hisat2_build_memory: int,
    gencode: bool,
    gtf_extra_attributes: str,
    gtf_group_features: str,
    featurecounts_group_type: str,
    featurecounts_feature_type: str,
    trimmer: Trimmer,
    min_trimmed_reads: int,
    remove_ribo_rna: bool,
    with_umi: bool,
    umitools_extract_method: str,
    umitools_grouping_method: UMIToolsGrouping,
    umitools_dedup_stats: bool,
    aligner: Aligner,
    pseudo_aligner_kmer_size: int,
    bam_csi_index: bool,
    star_ignore_sjdbgtf: bool,
    min_mapped_reads: float,
    stringtie_ignore_gtf: bool,
    kallisto_quant_fraglen: int,
    kallisto_quant_fraglen_sd: int,
    save_merged_fastq: bool,
    save_umi_intermeds: bool,
    save_non_ribo_reads: bool,
    save_bbsplit_reads: bool,
    save_reference: bool,
    save_trimmed: bool,
    save_align_intermeds: bool,
    save_unaligned: bool,
    deseq2_vst: bool,
    rseqc_modules: str,
    skip_gtf_filter: bool,
    skip_gtf_transcript_filter: bool,
    skip_bbsplit: bool,
    skip_umi_extract: bool,
    skip_trimming: bool,
    skip_alignment: bool,
    skip_pseudo_alignment: bool,
    skip_markduplicates: bool,
    skip_bigwig: bool,
    skip_stringtie: bool,
    skip_fastqc: bool,
    skip_preseq: bool,
    skip_dupradar: bool,
    skip_qualimap: bool,
    skip_rseqc: bool,
    skip_biotype_qc: bool,
    skip_deseq2_qc: bool,
    skip_multiqc: bool,
    skip_qc: bool,
) -> str:
    """
    Execute the Nextflow RNA-seq pipeline with the provided parameters.

    This function sets up the Nextflow environment, prepares input files,
    constructs the Nextflow command with all specified parameters, and runs
    the pipeline. It also handles log file uploading after execution.

    Args:
        pvc_name (str): Name of the persistent volume claim.
        input (List[SampleSheet]): List of input samples.
        run_name (str): Name of the analysis run.
        outdir (LatchOutputDir): Output directory for results.
        genome_source (str): Source of the reference genome.
        ... (other parameters as described in the function signature)

    Returns:
        str: The name of the analysis run.

    Note:
        This function has a large number of parameters to control various
        aspects of the RNA-seq analysis pipeline. Refer to the nf-core/rnaseq
        documentation for detailed information on each parameter.
    """
    try:
        shared_dir = Path("/nf-workdir")
        rename_current_execution(str(run_name))

        # Create custom sample sheet
        input_samplesheet = custom_samplesheet_constructor(
            samples=input, shared_dir=shared_dir
        )

        # List of directories and files to ignore when copying
        ignore_list = [
            "latch",
            ".latch",
            "nextflow",
            ".nextflow",
            "work",
            "results",
            "miniconda",
            "anaconda3",
            "mambaforge",
        ]

        # Copy necessary files to the shared directory
        shutil.copytree(
            Path("/root"),
            shared_dir,
            ignore=lambda src, names: ignore_list,
            ignore_dangling_symlinks=True,
            dirs_exist_ok=True,
        )

        # Construct the Nextflow command
        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            "docker",
            "-resume",
            "-process.executor",
            "k8s",
            *get_flag(
                "input",
                input_samplesheet,
            ),
            *get_flag_defaults(
                "outdir", LatchOutputDir(f"{outdir.remote_path}/{run_name}"), None
            ),
        ]

        # Add genome-specific parameters if using Latch genome source
        if genome_source == "latch_genome_source":
            cmd += [
                "--fasta",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.fna",
                "--gtf",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.gtf",
                "--gene_bed",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.filtered.bed",
                "--transcript_fasta",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.transcripts.fa",
                "--splicesites",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/index/{latch_genome.name}.genomic.filtered.splice_sites.txt",
                "--star_index",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/index/star.tar.gz",
                "--hisat2_index",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/index/hisat2.tar.gz",
                "--rsem_index",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/index/rsem.tar.gz",
                "--salmon_index",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/index/salmon.tar.gz",
                "--kallisto_index",
                f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/index/kallisto",
            ]

        # Add all other parameters to the command
        cmd += [
            *get_flag_defaults("fasta", fasta, None),
            *get_flag_defaults("gtf", gtf, None),
            *get_flag_defaults("gff", gff, None),
            *get_flag_defaults("gene_bed", gene_bed, None),
            *get_flag_defaults("transcript_fasta", transcript_fasta, None),
            *get_flag_defaults("additional_fasta", additional_fasta, None),
            *get_flag_defaults("splicesites", splicesites, None),
            *get_flag_defaults("star_index", star_index, None),
            *get_flag_defaults("hisat2_index", hisat2_index, None),
            *get_flag_defaults("rsem_index", rsem_index, None),
            *get_flag_defaults("salmon_index", salmon_index, None),
            *get_flag_defaults("kallisto_index", kallisto_index, None),
            *get_flag_defaults("hisat2_build_memory", hisat2_build_memory, 200),
            *get_flag_defaults("gencode", gencode, False),
            *get_flag_defaults(
                "gtf_extra_attributes", gtf_extra_attributes, "gene_name"
            ),
            *get_flag_defaults("gtf_group_features", gtf_group_features, "gene_id"),
            *get_flag_defaults(
                "featurecounts_group_type", featurecounts_group_type, "gene_biotype"
            ),
            *get_flag_defaults(
                "featurecounts_feature_type", featurecounts_feature_type, "exon"
            ),
            *get_flag_defaults("trimmer", trimmer, Trimmer.trimgalore),
            *get_flag_defaults("extra_trimgalore_args", extra_trimgalore_args, None),
            *get_flag_defaults("extra_fastp_args", extra_fastp_args, None),
            *get_flag_defaults("min_trimmed_reads", min_trimmed_reads, 10000),
            *get_flag_defaults("bbsplit_fasta_list", bbsplit_fasta_list, None),
            *get_flag_defaults("bbsplit_index", bbsplit_index, None),
            *get_flag_defaults("remove_ribo_rna", remove_ribo_rna, False),
            *get_flag_defaults("ribo_database_manifest", ribo_database_manifest, None),
            *get_flag_defaults("with_umi", with_umi, False),
            *get_flag_defaults(
                "umitools_extract_method", umitools_extract_method, "string"
            ),
            *get_flag_defaults("umitools_bc_pattern", umitools_bc_pattern, None),
            *get_flag_defaults("umitools_bc_pattern2", umitools_bc_pattern2, None),
            *get_flag_defaults("umi_discard_read", umi_discard_read, None),
            *get_flag_defaults("umitools_umi_separator", umitools_umi_separator, None),
            *get_flag_defaults(
                "umitools_grouping_method",
                umitools_grouping_method,
                UMIToolsGrouping.directional,
            ),
            *get_flag_defaults("umitools_dedup_stats", umitools_dedup_stats, False),
            *get_flag_defaults("aligner", aligner, Aligner.star_salmon),
            *get_flag_defaults("pseudo_aligner", pseudo_aligner, None),
            *get_flag_defaults(
                "pseudo_aligner_kmer_size", pseudo_aligner_kmer_size, 31
            ),
            *get_flag_defaults("bam_csi_index", bam_csi_index, False),
            *get_flag_defaults("star_ignore_sjdbgtf", star_ignore_sjdbgtf, False),
            *get_flag_defaults("salmon_quant_libtype", salmon_quant_libtype, None),
            *get_flag_defaults("min_mapped_reads", min_mapped_reads, 5.0),
            *get_flag_defaults("seq_center", seq_center, None),
            *get_flag_defaults("stringtie_ignore_gtf", stringtie_ignore_gtf, False),
            *get_flag_defaults("extra_star_align_args", extra_star_align_args, None),
            *get_flag_defaults(
                "extra_salmon_quant_args", extra_salmon_quant_args, None
            ),
            *get_flag_defaults(
                "extra_kallisto_quant_args", extra_kallisto_quant_args, None
            ),
            *get_flag_defaults("kallisto_quant_fraglen", kallisto_quant_fraglen, 200),
            *get_flag_defaults(
                "kallisto_quant_fraglen_sd", kallisto_quant_fraglen_sd, 200
            ),
            *get_flag_defaults("save_merged_fastq", save_merged_fastq, False),
            *get_flag_defaults("save_umi_intermeds", save_umi_intermeds, False),
            *get_flag_defaults("save_non_ribo_reads", save_non_ribo_reads, False),
            *get_flag_defaults("save_bbsplit_reads", save_bbsplit_reads, False),
            *get_flag_defaults("save_reference", save_reference, False),
            *get_flag_defaults("save_trimmed", save_trimmed, False),
            *get_flag_defaults("save_align_intermeds", save_align_intermeds, False),
            *get_flag_defaults("save_unaligned", save_unaligned, False),
            *get_flag_defaults("deseq2_vst", deseq2_vst, True),
            *get_flag_defaults(
                "rseqc_modules",
                rseqc_modules,
                "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication",
            ),
            *get_flag_defaults("skip_gtf_filter", skip_gtf_filter, False),
            *get_flag_defaults(
                "skip_gtf_transcript_filter", skip_gtf_transcript_filter, False
            ),
            *get_flag_defaults("skip_bbsplit", skip_bbsplit, True),
            *get_flag_defaults("skip_umi_extract", skip_umi_extract, False),
            *get_flag_defaults("skip_trimming", skip_trimming, False),
            *get_flag_defaults("skip_alignment", skip_alignment, False),
            *get_flag_defaults("skip_pseudo_alignment", skip_pseudo_alignment, False),
            *get_flag_defaults("skip_markduplicates", skip_markduplicates, False),
            *get_flag_defaults("skip_bigwig", skip_bigwig, False),
            *get_flag_defaults("skip_stringtie", skip_stringtie, False),
            *get_flag_defaults("skip_fastqc", skip_fastqc, False),
            *get_flag_defaults("skip_preseq", skip_preseq, True),
            *get_flag_defaults("skip_dupradar", skip_dupradar, False),
            *get_flag_defaults("skip_qualimap", skip_qualimap, False),
            *get_flag_defaults("skip_rseqc", skip_rseqc, False),
            *get_flag_defaults("skip_biotype_qc", skip_biotype_qc, False),
            *get_flag_defaults("skip_deseq2_qc", skip_deseq2_qc, False),
            *get_flag_defaults("skip_multiqc", skip_multiqc, False),
            *get_flag_defaults("skip_qc", skip_qc, False),
            *get_flag_defaults("email", email, None),
            *get_flag_defaults("multiqc_title", multiqc_title, None),
        ]

        print("Launching Nextflow Runtime")
        print(" ".join(cmd))
        print(flush=True)

        # Set up environment variables for Nextflow
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }

        # Run the Nextflow command
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )

        # Copying multiqc report
        try:
            time.sleep(10)
            if skip_alignment is True:
                # Pseudo aligner case
                multiqc_src_path = LPath(
                    f"{outdir.remote_path}/{run_name}/multiqc/multiqc_report.html"
                )
            else:
                multiqc_src_path = LPath(
                    f"{outdir.remote_path}/{run_name}/multiqc/{aligner.value}/multiqc_report.html"
                )

            multiqc_dst_path = LPath(
                f"{outdir.remote_path}/{run_name}/{run_name}_report.html"
            )
            multiqc_src_path.copy_to(multiqc_dst_path)
        except Exception as e:
            print(e)
            print("Could not copy multiqc report to outer folder.")

        return run_name

    finally:
        print()

        # Upload Nextflow log file
        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(
                    urljoins(
                        "latch:///Bulk_RNAseq_logs/nfcore_rnaseq", name, "nextflow.log"
                    )
                )
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)

        #
        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")
