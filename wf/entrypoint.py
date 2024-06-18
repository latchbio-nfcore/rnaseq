from dataclasses import dataclass
from enum import Enum
import os
import subprocess
import requests
import shutil
from pathlib import Path
import typing
import typing_extensions

from latch.resources.workflow import workflow
from latch.resources.tasks import nextflow_runtime_task, custom_task
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir
from latch.ldata.path import LPath
from latch_cli.nextflow.workflow import get_flag
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.utils import urljoins
from latch.types import metadata
from flytekit.core.annotation import FlyteAnnotation

from latch_cli.services.register.utils import import_module_by_path

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata

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
        }
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]






@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(pvc_name: str, input: LatchFile, outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], gtf: typing.Optional[LatchFile], gff: typing.Optional[LatchFile], gene_bed: typing.Optional[LatchFile], transcript_fasta: typing.Optional[LatchFile], additional_fasta: typing.Optional[LatchFile], splicesites: typing.Optional[LatchFile], star_index: typing.Optional[str], hisat2_index: typing.Optional[str], rsem_index: typing.Optional[str], salmon_index: typing.Optional[str], kallisto_index: typing.Optional[str], gencode: typing.Optional[bool], extra_trimgalore_args: typing.Optional[str], extra_fastp_args: typing.Optional[str], bbsplit_fasta_list: typing.Optional[LatchFile], bbsplit_index: typing.Optional[str], remove_ribo_rna: typing.Optional[bool], ribo_database_manifest: typing.Optional[LatchFile], with_umi: typing.Optional[bool], umitools_bc_pattern: typing.Optional[str], umitools_bc_pattern2: typing.Optional[str], umi_discard_read: typing.Optional[int], umitools_umi_separator: typing.Optional[str], umitools_dedup_stats: typing.Optional[bool], pseudo_aligner: typing.Optional[str], bam_csi_index: typing.Optional[bool], star_ignore_sjdbgtf: typing.Optional[bool], salmon_quant_libtype: typing.Optional[str], seq_center: typing.Optional[str], stringtie_ignore_gtf: typing.Optional[bool], extra_star_align_args: typing.Optional[str], extra_salmon_quant_args: typing.Optional[str], extra_kallisto_quant_args: typing.Optional[str], save_merged_fastq: typing.Optional[bool], save_umi_intermeds: typing.Optional[bool], save_non_ribo_reads: typing.Optional[bool], save_bbsplit_reads: typing.Optional[bool], save_reference: typing.Optional[bool], save_trimmed: typing.Optional[bool], save_align_intermeds: typing.Optional[bool], save_unaligned: typing.Optional[bool], skip_gtf_filter: typing.Optional[bool], skip_gtf_transcript_filter: typing.Optional[bool], skip_umi_extract: typing.Optional[bool], skip_trimming: typing.Optional[bool], skip_alignment: typing.Optional[bool], skip_pseudo_alignment: typing.Optional[bool], skip_markduplicates: typing.Optional[bool], skip_bigwig: typing.Optional[bool], skip_stringtie: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_dupradar: typing.Optional[bool], skip_qualimap: typing.Optional[bool], skip_rseqc: typing.Optional[bool], skip_biotype_qc: typing.Optional[bool], skip_deseq2_qc: typing.Optional[bool], skip_multiqc: typing.Optional[bool], skip_qc: typing.Optional[bool], multiqc_methods_description: typing.Optional[LatchFile], hisat2_build_memory: typing.Optional[str], gtf_extra_attributes: typing.Optional[str], gtf_group_features: typing.Optional[str], featurecounts_group_type: typing.Optional[str], featurecounts_feature_type: typing.Optional[str], trimmer: typing.Optional[str], min_trimmed_reads: typing.Optional[int], umitools_extract_method: typing.Optional[str], umitools_grouping_method: typing.Optional[str], aligner: typing.Optional[str], pseudo_aligner_kmer_size: typing.Optional[int], min_mapped_reads: typing.Optional[float], kallisto_quant_fraglen: typing.Optional[int], kallisto_quant_fraglen_sd: typing.Optional[int], deseq2_vst: typing.Optional[bool], rseqc_modules: typing.Optional[str], skip_bbsplit: typing.Optional[bool], skip_preseq: typing.Optional[bool]) -> None:
    try:
        shared_dir = Path("/nf-workdir")



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
            "-c",
            "latch.config",
                *get_flag('input', input),
                *get_flag('outdir', outdir),
                *get_flag('email', email),
                *get_flag('multiqc_title', multiqc_title),
                *get_flag('genome', genome),
                *get_flag('fasta', fasta),
                *get_flag('gtf', gtf),
                *get_flag('gff', gff),
                *get_flag('gene_bed', gene_bed),
                *get_flag('transcript_fasta', transcript_fasta),
                *get_flag('additional_fasta', additional_fasta),
                *get_flag('splicesites', splicesites),
                *get_flag('star_index', star_index),
                *get_flag('hisat2_index', hisat2_index),
                *get_flag('rsem_index', rsem_index),
                *get_flag('salmon_index', salmon_index),
                *get_flag('kallisto_index', kallisto_index),
                *get_flag('hisat2_build_memory', hisat2_build_memory),
                *get_flag('gencode', gencode),
                *get_flag('gtf_extra_attributes', gtf_extra_attributes),
                *get_flag('gtf_group_features', gtf_group_features),
                *get_flag('featurecounts_group_type', featurecounts_group_type),
                *get_flag('featurecounts_feature_type', featurecounts_feature_type),
                *get_flag('trimmer', trimmer),
                *get_flag('extra_trimgalore_args', extra_trimgalore_args),
                *get_flag('extra_fastp_args', extra_fastp_args),
                *get_flag('min_trimmed_reads', min_trimmed_reads),
                *get_flag('bbsplit_fasta_list', bbsplit_fasta_list),
                *get_flag('bbsplit_index', bbsplit_index),
                *get_flag('remove_ribo_rna', remove_ribo_rna),
                *get_flag('ribo_database_manifest', ribo_database_manifest),
                *get_flag('with_umi', with_umi),
                *get_flag('umitools_extract_method', umitools_extract_method),
                *get_flag('umitools_bc_pattern', umitools_bc_pattern),
                *get_flag('umitools_bc_pattern2', umitools_bc_pattern2),
                *get_flag('umi_discard_read', umi_discard_read),
                *get_flag('umitools_umi_separator', umitools_umi_separator),
                *get_flag('umitools_grouping_method', umitools_grouping_method),
                *get_flag('umitools_dedup_stats', umitools_dedup_stats),
                *get_flag('aligner', aligner),
                *get_flag('pseudo_aligner', pseudo_aligner),
                *get_flag('pseudo_aligner_kmer_size', pseudo_aligner_kmer_size),
                *get_flag('bam_csi_index', bam_csi_index),
                *get_flag('star_ignore_sjdbgtf', star_ignore_sjdbgtf),
                *get_flag('salmon_quant_libtype', salmon_quant_libtype),
                *get_flag('min_mapped_reads', min_mapped_reads),
                *get_flag('seq_center', seq_center),
                *get_flag('stringtie_ignore_gtf', stringtie_ignore_gtf),
                *get_flag('extra_star_align_args', extra_star_align_args),
                *get_flag('extra_salmon_quant_args', extra_salmon_quant_args),
                *get_flag('extra_kallisto_quant_args', extra_kallisto_quant_args),
                *get_flag('kallisto_quant_fraglen', kallisto_quant_fraglen),
                *get_flag('kallisto_quant_fraglen_sd', kallisto_quant_fraglen_sd),
                *get_flag('save_merged_fastq', save_merged_fastq),
                *get_flag('save_umi_intermeds', save_umi_intermeds),
                *get_flag('save_non_ribo_reads', save_non_ribo_reads),
                *get_flag('save_bbsplit_reads', save_bbsplit_reads),
                *get_flag('save_reference', save_reference),
                *get_flag('save_trimmed', save_trimmed),
                *get_flag('save_align_intermeds', save_align_intermeds),
                *get_flag('save_unaligned', save_unaligned),
                *get_flag('deseq2_vst', deseq2_vst),
                *get_flag('rseqc_modules', rseqc_modules),
                *get_flag('skip_gtf_filter', skip_gtf_filter),
                *get_flag('skip_gtf_transcript_filter', skip_gtf_transcript_filter),
                *get_flag('skip_bbsplit', skip_bbsplit),
                *get_flag('skip_umi_extract', skip_umi_extract),
                *get_flag('skip_trimming', skip_trimming),
                *get_flag('skip_alignment', skip_alignment),
                *get_flag('skip_pseudo_alignment', skip_pseudo_alignment),
                *get_flag('skip_markduplicates', skip_markduplicates),
                *get_flag('skip_bigwig', skip_bigwig),
                *get_flag('skip_stringtie', skip_stringtie),
                *get_flag('skip_fastqc', skip_fastqc),
                *get_flag('skip_preseq', skip_preseq),
                *get_flag('skip_dupradar', skip_dupradar),
                *get_flag('skip_qualimap', skip_qualimap),
                *get_flag('skip_rseqc', skip_rseqc),
                *get_flag('skip_biotype_qc', skip_biotype_qc),
                *get_flag('skip_deseq2_qc', skip_deseq2_qc),
                *get_flag('skip_multiqc', skip_multiqc),
                *get_flag('skip_qc', skip_qc),
                *get_flag('multiqc_methods_description', multiqc_methods_description)
        ]

        print("Launching Nextflow Runtime")
        print(' '.join(cmd))
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
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(urljoins("latch:///your_log_dir/nf_nf_core_rnaseq", name, "nextflow.log"))
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)



@workflow(metadata._nextflow_metadata)
def nf_nf_core_rnaseq(input: LatchFile, outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})], email: typing.Optional[str], multiqc_title: typing.Optional[str], genome: typing.Optional[str], fasta: typing.Optional[LatchFile], gtf: typing.Optional[LatchFile], gff: typing.Optional[LatchFile], gene_bed: typing.Optional[LatchFile], transcript_fasta: typing.Optional[LatchFile], additional_fasta: typing.Optional[LatchFile], splicesites: typing.Optional[LatchFile], star_index: typing.Optional[str], hisat2_index: typing.Optional[str], rsem_index: typing.Optional[str], salmon_index: typing.Optional[str], kallisto_index: typing.Optional[str], gencode: typing.Optional[bool], extra_trimgalore_args: typing.Optional[str], extra_fastp_args: typing.Optional[str], bbsplit_fasta_list: typing.Optional[LatchFile], bbsplit_index: typing.Optional[str], remove_ribo_rna: typing.Optional[bool], ribo_database_manifest: typing.Optional[LatchFile], with_umi: typing.Optional[bool], umitools_bc_pattern: typing.Optional[str], umitools_bc_pattern2: typing.Optional[str], umi_discard_read: typing.Optional[int], umitools_umi_separator: typing.Optional[str], umitools_dedup_stats: typing.Optional[bool], pseudo_aligner: typing.Optional[str], bam_csi_index: typing.Optional[bool], star_ignore_sjdbgtf: typing.Optional[bool], salmon_quant_libtype: typing.Optional[str], seq_center: typing.Optional[str], stringtie_ignore_gtf: typing.Optional[bool], extra_star_align_args: typing.Optional[str], extra_salmon_quant_args: typing.Optional[str], extra_kallisto_quant_args: typing.Optional[str], save_merged_fastq: typing.Optional[bool], save_umi_intermeds: typing.Optional[bool], save_non_ribo_reads: typing.Optional[bool], save_bbsplit_reads: typing.Optional[bool], save_reference: typing.Optional[bool], save_trimmed: typing.Optional[bool], save_align_intermeds: typing.Optional[bool], save_unaligned: typing.Optional[bool], skip_gtf_filter: typing.Optional[bool], skip_gtf_transcript_filter: typing.Optional[bool], skip_umi_extract: typing.Optional[bool], skip_trimming: typing.Optional[bool], skip_alignment: typing.Optional[bool], skip_pseudo_alignment: typing.Optional[bool], skip_markduplicates: typing.Optional[bool], skip_bigwig: typing.Optional[bool], skip_stringtie: typing.Optional[bool], skip_fastqc: typing.Optional[bool], skip_dupradar: typing.Optional[bool], skip_qualimap: typing.Optional[bool], skip_rseqc: typing.Optional[bool], skip_biotype_qc: typing.Optional[bool], skip_deseq2_qc: typing.Optional[bool], skip_multiqc: typing.Optional[bool], skip_qc: typing.Optional[bool], multiqc_methods_description: typing.Optional[LatchFile], hisat2_build_memory: typing.Optional[str] = '200.GB', gtf_extra_attributes: typing.Optional[str] = 'gene_name', gtf_group_features: typing.Optional[str] = 'gene_id', featurecounts_group_type: typing.Optional[str] = 'gene_biotype', featurecounts_feature_type: typing.Optional[str] = 'exon', trimmer: typing.Optional[str] = 'trimgalore', min_trimmed_reads: typing.Optional[int] = 10000, umitools_extract_method: typing.Optional[str] = 'string', umitools_grouping_method: typing.Optional[str] = 'directional', aligner: typing.Optional[str] = 'star_salmon', pseudo_aligner_kmer_size: typing.Optional[int] = 31, min_mapped_reads: typing.Optional[float] = 5.0, kallisto_quant_fraglen: typing.Optional[int] = 200, kallisto_quant_fraglen_sd: typing.Optional[int] = 200, deseq2_vst: typing.Optional[bool] = True, rseqc_modules: typing.Optional[str] = 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication', skip_bbsplit: typing.Optional[bool] = True, skip_preseq: typing.Optional[bool] = True) -> None:
    """
    nf-core/rnaseq

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(pvc_name=pvc_name, input=input, outdir=outdir, email=email, multiqc_title=multiqc_title, genome=genome, fasta=fasta, gtf=gtf, gff=gff, gene_bed=gene_bed, transcript_fasta=transcript_fasta, additional_fasta=additional_fasta, splicesites=splicesites, star_index=star_index, hisat2_index=hisat2_index, rsem_index=rsem_index, salmon_index=salmon_index, kallisto_index=kallisto_index, hisat2_build_memory=hisat2_build_memory, gencode=gencode, gtf_extra_attributes=gtf_extra_attributes, gtf_group_features=gtf_group_features, featurecounts_group_type=featurecounts_group_type, featurecounts_feature_type=featurecounts_feature_type, trimmer=trimmer, extra_trimgalore_args=extra_trimgalore_args, extra_fastp_args=extra_fastp_args, min_trimmed_reads=min_trimmed_reads, bbsplit_fasta_list=bbsplit_fasta_list, bbsplit_index=bbsplit_index, remove_ribo_rna=remove_ribo_rna, ribo_database_manifest=ribo_database_manifest, with_umi=with_umi, umitools_extract_method=umitools_extract_method, umitools_bc_pattern=umitools_bc_pattern, umitools_bc_pattern2=umitools_bc_pattern2, umi_discard_read=umi_discard_read, umitools_umi_separator=umitools_umi_separator, umitools_grouping_method=umitools_grouping_method, umitools_dedup_stats=umitools_dedup_stats, aligner=aligner, pseudo_aligner=pseudo_aligner, pseudo_aligner_kmer_size=pseudo_aligner_kmer_size, bam_csi_index=bam_csi_index, star_ignore_sjdbgtf=star_ignore_sjdbgtf, salmon_quant_libtype=salmon_quant_libtype, min_mapped_reads=min_mapped_reads, seq_center=seq_center, stringtie_ignore_gtf=stringtie_ignore_gtf, extra_star_align_args=extra_star_align_args, extra_salmon_quant_args=extra_salmon_quant_args, extra_kallisto_quant_args=extra_kallisto_quant_args, kallisto_quant_fraglen=kallisto_quant_fraglen, kallisto_quant_fraglen_sd=kallisto_quant_fraglen_sd, save_merged_fastq=save_merged_fastq, save_umi_intermeds=save_umi_intermeds, save_non_ribo_reads=save_non_ribo_reads, save_bbsplit_reads=save_bbsplit_reads, save_reference=save_reference, save_trimmed=save_trimmed, save_align_intermeds=save_align_intermeds, save_unaligned=save_unaligned, deseq2_vst=deseq2_vst, rseqc_modules=rseqc_modules, skip_gtf_filter=skip_gtf_filter, skip_gtf_transcript_filter=skip_gtf_transcript_filter, skip_bbsplit=skip_bbsplit, skip_umi_extract=skip_umi_extract, skip_trimming=skip_trimming, skip_alignment=skip_alignment, skip_pseudo_alignment=skip_pseudo_alignment, skip_markduplicates=skip_markduplicates, skip_bigwig=skip_bigwig, skip_stringtie=skip_stringtie, skip_fastqc=skip_fastqc, skip_preseq=skip_preseq, skip_dupradar=skip_dupradar, skip_qualimap=skip_qualimap, skip_rseqc=skip_rseqc, skip_biotype_qc=skip_biotype_qc, skip_deseq2_qc=skip_deseq2_qc, skip_multiqc=skip_multiqc, skip_qc=skip_qc, multiqc_methods_description=multiqc_methods_description)

