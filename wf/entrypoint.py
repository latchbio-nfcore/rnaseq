import csv
import os
import shutil
import subprocess
import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, List, Optional

import requests
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins
from wf.prep_dge import prep_dge

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)

with open("./README.md", "r") as readme_file:
    readme_contents = readme_file.read()


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    strandedness: str


class Reference_Type(Enum):
    homo_sapiens = "Homo sapiens (RefSeq GRCh38.p14)"
    mus_musculus = "Mus musculus (RefSeq GRCm39)"
    rattus_norvegicus = "Rattus norvegicus (RefSeq GRCr8)"
    # drosophila_melanogaster = "Drosophila melanogaster (RefSeq Release_6_plus_ISO1_MT)"
    # rhesus_macaque = "Macaca mulatta (RefSeq rheMac10/Mmul_10)"
    saccharomyces_cerevisiae = "Saccharomyces cerevisiae (RefSeq R64)"


class Trimmer(Enum):
    trimgalore = "trimgalore"
    fastp = "fastp"


class UMIToolsGrouping(Enum):
    directional = "directional"
    unique = "unique"
    cluster = "cluster"
    percentile = "percentile"
    adjacency = "adjacency"


class Aligner(Enum):
    star_salmon = "star_salmon"
    star_rsem = "star_rsem"
    hisat2 = "hisat2"


class SalmonQuantLibType(Enum):
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
    salmon = "salmon"
    kallisto = "kallisto"


class DifferentialGeneTool(Enum):
    deseq2 = "deseq2"
    # sleuth = "sleuth"
    # edgeR = "edgeR"


def get_flag_defaults(name: str, val: Any, default_val: Optional[Any]):
    if val == default_val or val is None:
        return ""
    else:
        return get_flag(name=name, val=val)


def custom_samplesheet_constructor(samples: List[SampleSheet]) -> Path:
    samplesheet = Path("/root/samplesheet.csv")

    columns = ["sample", "fastq_1", "fastq_2", "strandedness"]

    with open(samplesheet, "w") as f:
        writer = csv.DictWriter(f, columns, delimiter=",")
        writer.writeheader()

        for sample in samples:
            row_data = {
                "sample": sample.sample,
                "fastq_1": sample.fastq_1.remote_path,
                "fastq_2": "" if sample.fastq_2 is None else sample.fastq_2.remote_path,
                "strandedness": sample.strandedness,
            }
            writer.writerow(row_data)

    return samplesheet


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_gib": 100,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@nextflow_runtime_task(cpu=4, memory=16, storage_gib=1000)
def nextflow_runtime(
    pvc_name: str,
    input: typing.List[SampleSheet],
    run_name: str,
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
    try:
        shared_dir = Path("/nf-workdir")

        input_samplesheet = custom_samplesheet_constructor(samples=input)

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

        shutil.copytree(
            Path("/root"),
            shared_dir,
            ignore=lambda src, names: ignore_list,
            ignore_dangling_symlinks=True,
            dirs_exist_ok=True,
        )

        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            "docker",
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

        env = {
            **os.environ,
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms2048M -Xmx8G -XX:ActiveProcessorCount=4",
            "K8S_STORAGE_CLAIM_NAME": pvc_name,
            "NXF_DISABLE_CHECK_LATEST": "true",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )

        return run_name

    finally:
        print()

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


@workflow(metadata._nextflow_metadata)
def nf_nf_core_rnaseq(
    input: typing.List[SampleSheet],
    run_name: str,
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
    outdir: LatchOutputDir = LatchOutputDir("latch:///Bulk_RNAseq"),
    latch_genome: Reference_Type = Reference_Type.homo_sapiens,
    hisat2_build_memory: int = 200,
    gencode: bool = False,
    gtf_extra_attributes: str = "gene_name",
    gtf_group_features: str = "gene_id",
    featurecounts_group_type: str = "gene_biotype",
    featurecounts_feature_type: str = "exon",
    trimmer: Trimmer = Trimmer.trimgalore,
    min_trimmed_reads: int = 10000,
    remove_ribo_rna: bool = False,
    with_umi: bool = False,
    umitools_extract_method: str = "string",
    umitools_grouping_method: UMIToolsGrouping = UMIToolsGrouping.directional,
    umitools_dedup_stats: bool = False,
    aligner: Aligner = Aligner.star_salmon,
    pseudo_aligner_kmer_size: int = 31,
    bam_csi_index: bool = False,
    star_ignore_sjdbgtf: bool = False,
    min_mapped_reads: float = 5.0,
    stringtie_ignore_gtf: bool = False,
    kallisto_quant_fraglen: int = 200,
    kallisto_quant_fraglen_sd: int = 200,
    save_merged_fastq: bool = False,
    save_umi_intermeds: bool = False,
    save_non_ribo_reads: bool = False,
    save_bbsplit_reads: bool = False,
    save_reference: bool = False,
    save_trimmed: bool = False,
    save_align_intermeds: bool = False,
    save_unaligned: bool = False,
    deseq2_vst: bool = True,
    rseqc_modules: str = "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication",
    skip_gtf_filter: bool = False,
    skip_gtf_transcript_filter: bool = False,
    skip_bbsplit: bool = True,
    skip_umi_extract: bool = False,
    skip_trimming: bool = False,
    skip_alignment: bool = False,
    skip_pseudo_alignment: bool = False,
    skip_markduplicates: bool = False,
    skip_bigwig: bool = False,
    skip_stringtie: bool = False,
    skip_fastqc: bool = False,
    skip_preseq: bool = True,
    skip_dupradar: bool = False,
    skip_qualimap: bool = False,
    skip_rseqc: bool = False,
    skip_biotype_qc: bool = False,
    skip_deseq2_qc: bool = False,
    skip_multiqc: bool = False,
    skip_qc: bool = False,
) -> LatchOutputDir:
    f"""
    nf-core/rnaseq

    test test
    {readme_contents}
    """

    pvc_name: str = initialize()
    run_name = nextflow_runtime(
        pvc_name=pvc_name,
        input=input,
        run_name=run_name,
        outdir=outdir,
        genome_source=genome_source,
        latch_genome=latch_genome,
        fasta=fasta,
        gtf=gtf,
        gff=gff,
        gene_bed=gene_bed,
        transcript_fasta=transcript_fasta,
        additional_fasta=additional_fasta,
        splicesites=splicesites,
        star_index=star_index,
        hisat2_index=hisat2_index,
        rsem_index=rsem_index,
        salmon_index=salmon_index,
        kallisto_index=kallisto_index,
        hisat2_build_memory=hisat2_build_memory,
        gencode=gencode,
        gtf_extra_attributes=gtf_extra_attributes,
        gtf_group_features=gtf_group_features,
        featurecounts_group_type=featurecounts_group_type,
        featurecounts_feature_type=featurecounts_feature_type,
        trimmer=trimmer,
        extra_trimgalore_args=extra_trimgalore_args,
        extra_fastp_args=extra_fastp_args,
        min_trimmed_reads=min_trimmed_reads,
        bbsplit_fasta_list=bbsplit_fasta_list,
        bbsplit_index=bbsplit_index,
        remove_ribo_rna=remove_ribo_rna,
        ribo_database_manifest=ribo_database_manifest,
        with_umi=with_umi,
        umitools_extract_method=umitools_extract_method,
        umitools_bc_pattern=umitools_bc_pattern,
        umitools_bc_pattern2=umitools_bc_pattern2,
        umi_discard_read=umi_discard_read,
        umitools_umi_separator=umitools_umi_separator,
        umitools_grouping_method=umitools_grouping_method,
        umitools_dedup_stats=umitools_dedup_stats,
        aligner=aligner,
        pseudo_aligner=pseudo_aligner,
        pseudo_aligner_kmer_size=pseudo_aligner_kmer_size,
        bam_csi_index=bam_csi_index,
        star_ignore_sjdbgtf=star_ignore_sjdbgtf,
        salmon_quant_libtype=salmon_quant_libtype,
        min_mapped_reads=min_mapped_reads,
        seq_center=seq_center,
        stringtie_ignore_gtf=stringtie_ignore_gtf,
        extra_star_align_args=extra_star_align_args,
        extra_salmon_quant_args=extra_salmon_quant_args,
        extra_kallisto_quant_args=extra_kallisto_quant_args,
        kallisto_quant_fraglen=kallisto_quant_fraglen,
        kallisto_quant_fraglen_sd=kallisto_quant_fraglen_sd,
        save_merged_fastq=save_merged_fastq,
        save_umi_intermeds=save_umi_intermeds,
        save_non_ribo_reads=save_non_ribo_reads,
        save_bbsplit_reads=save_bbsplit_reads,
        save_reference=save_reference,
        save_trimmed=save_trimmed,
        save_align_intermeds=save_align_intermeds,
        save_unaligned=save_unaligned,
        deseq2_vst=deseq2_vst,
        rseqc_modules=rseqc_modules,
        skip_gtf_filter=skip_gtf_filter,
        skip_gtf_transcript_filter=skip_gtf_transcript_filter,
        skip_bbsplit=skip_bbsplit,
        skip_umi_extract=skip_umi_extract,
        skip_trimming=skip_trimming,
        skip_alignment=skip_alignment,
        skip_pseudo_alignment=skip_pseudo_alignment,
        skip_markduplicates=skip_markduplicates,
        skip_bigwig=skip_bigwig,
        skip_stringtie=skip_stringtie,
        skip_fastqc=skip_fastqc,
        skip_preseq=skip_preseq,
        skip_dupradar=skip_dupradar,
        skip_qualimap=skip_qualimap,
        skip_rseqc=skip_rseqc,
        skip_biotype_qc=skip_biotype_qc,
        skip_deseq2_qc=skip_deseq2_qc,
        skip_multiqc=skip_multiqc,
        skip_qc=skip_qc,
        email=email,
        multiqc_title=multiqc_title,
    )

    return prep_dge(
        run_name=run_name,
        latch_genome=latch_genome,
        outdir=outdir,
    )
