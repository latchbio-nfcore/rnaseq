import csv
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchAuthor,
    NextflowMetadata,
    NextflowParameter,
    NextflowRuntimeResources,
    Params,
    Section,
    Spoiler,
    Text,
)

flow = [
    Section(
        "Samples",
        Params(
            "input",
        ),
    ),
    Section(
        "Reference Genome",
        Fork(
            "genome_source",
            "",
            latch_genome_source=ForkBranch(
                "Latch Verified Reference Genome",
                Params(
                    "latch_genome",
                ),
            ),
            custom=ForkBranch(
                "Custom Reference Genome",
                Params(
                    "fasta",
                    "gtf",
                ),
            ),
        ),
    ),
    Section(
        "Alignment",
        Fork(
            "alignment_method",
            "",
            traditional_alignment=ForkBranch(
                "Traditional Alignment",
                Params(
                    "aligner",
                ),
            ),
            pseudo_alignment=ForkBranch(
                "Pseudo Aignment",
                Text(
                    "A pseudoaligner is run in addition to the standard alignment workflow defined by --aligner for quantification."
                ),
                Text("The default is Salmon but this can be changed to Kallisto."),
                Text(
                    "Toggle on the '--skip_alignment' parameter to run Salmon or Kallisto in isolation."
                ),
                Params(
                    "skip_alignment",
                    "pseudo_aligner",
                ),
            ),
        ),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "Optional Arguments",
        Text("Additional optional arguments"),
        Section(
            "General Options",
            Params(
                "email",
                "multiqc_title",
            ),
        ),
        Section(
            "Reference Files",
            Params(
                "gff",
                "gene_bed",
                "transcript_fasta",
                "additional_fasta",
                "splicesites",
                "star_index",
                "hisat2_index",
                "rsem_index",
                "salmon_index",
                "kallisto_index",
                "hisat2_build_memory",
                "gencode",
                "gtf_extra_attributes",
                "gtf_group_features",
                "featurecounts_group_type",
                "featurecounts_feature_type",
            ),
        ),
        Section(
            "Trimming Options",
            Params(
                "trimmer",
                "extra_trimgalore_args",
                "extra_fastp_args",
                "min_trimmed_reads",
            ),
        ),
        Section(
            "Read Filtering",
            Params(
                "bbsplit_fasta_list",
                "bbsplit_index",
                "remove_ribo_rna",
                "ribo_database_manifest",
            ),
        ),
        Section(
            "UMI Options",
            Params(
                "with_umi",
                "umitools_extract_method",
                "umitools_bc_pattern",
                "umitools_bc_pattern2",
                "umi_discard_read",
                "umitools_umi_separator",
                "umitools_grouping_method",
                "umitools_dedup_stats",
            ),
        ),
        Section(
            "Alignment Options",
            Params(
                "pseudo_aligner_kmer_size",
                "bam_csi_index",
                "star_ignore_sjdbgtf",
                "salmon_quant_libtype",
                "min_mapped_reads",
                "seq_center",
                "stringtie_ignore_gtf",
                "extra_star_align_args",
                "extra_salmon_quant_args",
                "extra_kallisto_quant_args",
                "kallisto_quant_fraglen",
                "kallisto_quant_fraglen_sd",
            ),
        ),
        Section(
            "Optional Outputs",
            Params(
                "save_merged_fastq",
                "save_umi_intermeds",
                "save_non_ribo_reads",
                "save_bbsplit_reads",
                "save_reference",
                "save_trimmed",
                "save_align_intermeds",
                "save_unaligned",
            ),
        ),
        Section(
            "Quality Control",
            Params(
                "deseq2_vst",
                "rseqc_modules",
            ),
        ),
        Section(
            "Process Skipping",
            Params(
                "skip_gtf_filter",
                "skip_gtf_transcript_filter",
                "skip_bbsplit",
                "skip_umi_extract",
                "skip_trimming",
                "skip_pseudo_alignment",
                "skip_markduplicates",
                "skip_bigwig",
                "skip_stringtie",
                "skip_fastqc",
                "skip_preseq",
                "skip_dupradar",
                "skip_qualimap",
                "skip_rseqc",
                "skip_biotype_qc",
                "skip_deseq2_qc",
                "skip_multiqc",
                "skip_qc",
            ),
        ),
    ),
]


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


NextflowMetadata(
    display_name="nf-core/rnaseq",
    author=LatchAuthor(
        name="nf-core",
    ),
    parameters={
        "input": NextflowParameter(
            type=List[SampleSheet],
            display_name="Samplesheet",
            description="Information about the samples in the experiment.",
            batch_table_column=True,
            samplesheet_type=custom_samplesheet_constructor,
            samplesheet=True,
        ),
        "run_name": NextflowParameter(
            type=str,
            display_name="Run Name",
            description="Name of run",
            batch_table_column=True,
        ),
        "outdir": NextflowParameter(
            type=LatchOutputDir,
            display_name="Output Directory",
            description="The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
            batch_table_column=True,
            default=LatchOutputDir("latch:///Bulk_RNAseq"),
        ),
        "genome_source": NextflowParameter(
            type=str,
            display_name="Reference Genome",
            description="Choose Reference Genome",
        ),
        "latch_genome": NextflowParameter(
            type=Reference_Type,
            display_name="Latch Verfied Reference Genome",
            description="Name of Latch Verfied Reference Genome.",
            default=Reference_Type.homo_sapiens,
        ),
        "fasta": NextflowParameter(
            type=Optional[LatchFile],
            display_name="FASTA Genome File",
            description="Path to FASTA genome file.",
        ),
        "gtf": NextflowParameter(
            type=Optional[LatchFile],
            display_name="GTF Annotation File",
            description="Path to GTF annotation file. --gff must be used in the Optional parameters if this is not available.",
        ),
        "gff": NextflowParameter(
            type=Optional[LatchFile],
            display_name="GFF3 Annotation File",
            description="Path to GFF3 annotation file. This parameter must be specified if --gtf is not specified.",
        ),
        "gene_bed": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Gene BED File",
            description="Path to BED file containing gene intervals. This will be created from the GTF file if not specified.",
        ),
        "transcript_fasta": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Transcriptome FASTA File",
            description="Path to FASTA transcriptome file.",
        ),
        "additional_fasta": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Additional FASTA File",
            description="FASTA file to concatenate to genome FASTA file e.g. containing spike-in sequences. If provided, the sequences in this file will get concatenated to the existing genome FASTA file, a GTF file will be automatically created using the entire sequence as the gene, transcript, and exon features, and any alignment index will get created from the combined FASTA and GTF. It is recommended to save the reference with --save_reference to re-use the index for future runs so you do not need to create it again.",
        ),
        "splicesites": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Splice Sites File",
            description="Splice sites file required for HISAT2.",
        ),
        "star_index": NextflowParameter(
            type=Optional[LatchFile],
            display_name="STAR Index",
            description="Path to directory or tar.gz archive for pre-built STAR index.",
        ),
        "hisat2_index": NextflowParameter(
            type=Optional[LatchFile],
            display_name="HISAT2 Index",
            description="Path to directory or tar.gz archive for pre-built HISAT2 index.",
        ),
        "rsem_index": NextflowParameter(
            type=Optional[LatchFile],
            display_name="RSEM Index",
            description="Path to directory or tar.gz archive for pre-built RSEM index.",
        ),
        "salmon_index": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Salmon Index",
            description="Path to directory or tar.gz archive for pre-built Salmon index.",
        ),
        "kallisto_index": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Kallisto Index",
            description="Path to directory or tar.gz archive for pre-built Kallisto index.",
        ),
        "hisat2_build_memory": NextflowParameter(
            type=int,
            display_name="HISAT2 Build Memory (GB)",
            description="Minimum memory required to use splice sites and exons in the HiSAT2 index build process.",
            default=200,
        ),
        "gencode": NextflowParameter(
            type=bool,
            display_name="GENCODE Format",
            description="Specify if your GTF annotation is in GENCODE format. If your GTF file is in GENCODE format and you would like to run Salmon i.e. --pseudo_aligner salmon, you will need to provide this parameter in order to build the Salmon index appropriately.",
            default=False,
        ),
        "gtf_extra_attributes": NextflowParameter(
            type=str,
            display_name="GTF Extra Attributes",
            description="By default, the pipeline uses the gene_name field to obtain additional gene identifiers from the input GTF file when running Salmon. This behaviour can be modified by specifying --gtf_extra_attributes when running the pipeline. Note that you can also specify more than one desired value, separated by a comma e.g. --gtf_extra_attributes gene_id,....",
            default="gene_name",
        ),
        "gtf_group_features": NextflowParameter(
            type=str,
            display_name="GTF Group Features",
            description="Define the attribute type used to group features in the GTF file when running Salmon.",
            default="gene_id",
        ),
        "featurecounts_group_type": NextflowParameter(
            type=str,
            display_name="FeatureCounts Group Type",
            description="The attribute type used to group feature types in the GTF file when generating the biotype plot with featureCounts.",
            default="gene_biotype",
        ),
        "featurecounts_feature_type": NextflowParameter(
            type=str,
            display_name="FeatureCounts Feature Type",
            description="By default, the pipeline assigns reads based on the 'exon' attribute within the GTF file.",
            default="exon",
        ),
        "trimmer": NextflowParameter(
            type=Trimmer,
            display_name="Trimming Tool",
            description="Specifies the trimming tool to use - available options are 'trimgalore' and 'fastp'.",
            default=Trimmer.trimgalore,
        ),
        "extra_trimgalore_args": NextflowParameter(
            type=Optional[str],
            display_name="Extra Trim Galore! Args",
            description="Extra arguments to pass to Trim Galore! command in addition to defaults defined by the pipeline.",
        ),
        "extra_fastp_args": NextflowParameter(
            type=Optional[str],
            display_name="Extra fastp Args",
            description="Extra arguments to pass to fastp command in addition to defaults defined by the pipeline.",
        ),
        "min_trimmed_reads": NextflowParameter(
            type=int,
            display_name="Minimum Trimmed Reads",
            description="Minimum number of trimmed reads below which samples are removed from further processing. Some downstream steps in the pipeline will fail if this threshold is too low.",
            default=10000,
        ),
        "bbsplit_fasta_list": NextflowParameter(
            type=Optional[LatchFile],
            display_name="BBSplit FASTA List",
            description="Comma-separated file containing a list of reference genomes to filter reads against with BBSplit. You have to also explicitly set --skip_bbsplit false if you want to use BBSplit. The file should contain 2 columns: short name and full path to reference genome(s) e.g.",
        ),
        "bbsplit_index": NextflowParameter(
            type=Optional[LatchFile],
            display_name="BBSplit Index",
            description="Path to directory or tar.gz archive for pre-built BBSplit index.",
        ),
        "remove_ribo_rna": NextflowParameter(
            type=bool,
            display_name="Remove Ribosomal RNA",
            description="Enable the removal of reads derived from ribosomal RNA using SortMeRNA.",
            default=False,
        ),
        "ribo_database_manifest": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Ribo Database Manifest",
            description="Text file containing paths to fasta files (one per line) that will be used to create the database for SortMeRNA.",
        ),
        "with_umi": NextflowParameter(
            type=bool,
            display_name="Enable UMI",
            description="Enable UMI-based read deduplication.",
            default=False,
        ),
        "umitools_extract_method": NextflowParameter(
            type=str,
            display_name="UMI Extract Method",
            description="UMI pattern to use. Can be either 'string' (default) or 'regex'.",
            default="string",
        ),
        "umitools_bc_pattern": NextflowParameter(
            type=Optional[str],
            display_name="UMI Barcode Pattern",
            description="The UMI barcode pattern to use e.g. 'NNNNNN' indicates that the first 6 nucleotides of the read are from the UMI.",
        ),
        "umitools_bc_pattern2": NextflowParameter(
            type=Optional[str],
            display_name="UMI Barcode Pattern 2",
            description="The UMI barcode pattern to use if the UMI is located in read 2.",
        ),
        "umi_discard_read": NextflowParameter(
            type=Optional[int],
            display_name="UMI Discard Read",
            description="After UMI barcode extraction discard either R1 or R2 by setting this parameter to 1 or 2, respectively.",
        ),
        "umitools_umi_separator": NextflowParameter(
            type=Optional[str],
            display_name="UMI Separator",
            description="The character that separates the UMI in the read name. Most likely a colon if you skipped the extraction with UMI-tools and used other software.",
        ),
        "umitools_grouping_method": NextflowParameter(
            type=UMIToolsGrouping,
            display_name="UMI Grouping Method",
            description="Method to use to determine read groups by subsuming those with similar UMIs. All methods start by identifying the reads with the same mapping position, but treat similar yet nonidentical UMIs differently.",
            default=UMIToolsGrouping.directional,
        ),
        "umitools_dedup_stats": NextflowParameter(
            type=bool,
            display_name="UMI Deduplication Stats",
            description="Generate output stats when running 'umi_tools dedup'.",
            default=False,
        ),
        "aligner": NextflowParameter(
            type=Aligner,
            display_name="Alignment Algorithm",
            description="Specifies the alignment algorithm to use - available options are 'star_salmon', 'star_rsem' and 'hisat2'.",
            default=Aligner.star_salmon,
        ),
        "pseudo_aligner": NextflowParameter(
            type=Optional[PseudoAligner],
            display_name="Pseudo Aligner",
            description="Specifies the pseudo aligner to use - available options are 'salmon'. Runs in addition to '--aligner'. Set 'skip_aignment' to be true if only pseudo alignemnt needed.",
        ),
        "pseudo_aligner_kmer_size": NextflowParameter(
            type=int,
            display_name="Pseudo Aligner Kmer Size",
            description="Kmer length passed to indexing step of pseudoaligners.",
            default=31,
        ),
        "bam_csi_index": NextflowParameter(
            type=bool,
            display_name="Create CSI Index",
            description="Create a CSI index for BAM files instead of the traditional BAI index. This will be required for genomes with larger chromosome sizes.",
            default=False,
        ),
        "star_ignore_sjdbgtf": NextflowParameter(
            type=bool,
            display_name="STAR Ignore SJDB GTF",
            description="When using pre-built STAR indices do not re-extract and use splice junctions from the GTF file.",
            default=False,
        ),
        "salmon_quant_libtype": NextflowParameter(
            type=Optional[SalmonQuantLibType],
            display_name="Salmon Library Type",
            description="Override Salmon library type inferred based on strandedness defined in meta object.",
        ),
        "min_mapped_reads": NextflowParameter(
            type=float,
            display_name="Minimum Mapped Reads",
            description="Minimum percentage of uniquely mapped reads below which samples are removed from further processing.",
            default=5.0,
        ),
        "seq_center": NextflowParameter(
            type=Optional[str],
            display_name="Sequencing Center",
            description="Sequencing center information to be added to read group of BAM files.",
        ),
        "stringtie_ignore_gtf": NextflowParameter(
            type=bool,
            display_name="StringTie Ignore GTF",
            description="Perform reference-guided de novo assembly of transcripts using StringTie i.e. dont restrict to those in GTF file.",
            default=False,
        ),
        "extra_star_align_args": NextflowParameter(
            type=Optional[str],
            display_name="Extra STAR Align Args",
            description="Extra arguments to pass to STAR alignment command in addition to defaults defined by the pipeline. Only available for the STAR-Salmon route.",
        ),
        "extra_salmon_quant_args": NextflowParameter(
            type=Optional[str],
            display_name="Extra Salmon Quant Args",
            description="Extra arguments to pass to Salmon quant command in addition to defaults defined by the pipeline.",
        ),
        "extra_kallisto_quant_args": NextflowParameter(
            type=Optional[str],
            display_name="Extra Kallisto Quant Args",
            description="Extra arguments to pass to Kallisto quant command in addition to defaults defined by the pipeline.",
        ),
        "kallisto_quant_fraglen": NextflowParameter(
            type=int,
            display_name="Kallisto Quant Fragment Length",
            description="In single-end mode Kallisto requires an estimated fragment length. Specify a default value for that here. TODO: use existing RSeQC results to do this dynamically.",
            default=200,
        ),
        "kallisto_quant_fraglen_sd": NextflowParameter(
            type=int,
            display_name="Kallisto Quant Fragment Length SD",
            description="In single-end mode, Kallisto requires an estimated standard error for fragment length. Specify a default value for that here. TODO: use existing RSeQC results to do this dynamically.",
            default=200,
        ),
        "save_merged_fastq": NextflowParameter(
            type=bool,
            display_name="Save Merged FastQ",
            description="Save FastQ files after merging re-sequenced libraries in the results directory.",
            default=False,
        ),
        "save_umi_intermeds": NextflowParameter(
            type=bool,
            display_name="Save UMI Intermediates",
            description="If this option is specified, intermediate FastQ and BAM files produced by UMI-tools are also saved in the results directory.",
            default=False,
        ),
        "save_non_ribo_reads": NextflowParameter(
            type=bool,
            display_name="Save Non-rRNA Reads",
            description="If this option is specified, intermediate FastQ files containing non-rRNA reads will be saved in the results directory.",
            default=False,
        ),
        "save_bbsplit_reads": NextflowParameter(
            type=bool,
            display_name="Save BBSplit Reads",
            description="If this option is specified, FastQ files split by reference will be saved in the results directory.",
            default=False,
        ),
        "save_reference": NextflowParameter(
            type=bool,
            display_name="Save Reference",
            description="If generated by the pipeline save the STAR index in the results directory.",
            default=False,
        ),
        "save_trimmed": NextflowParameter(
            type=bool,
            display_name="Save Trimmed",
            description="Save the trimmed FastQ files in the results directory.",
            default=False,
        ),
        "save_align_intermeds": NextflowParameter(
            type=bool,
            display_name="Save Alignment Intermediates",
            description="Save the intermediate BAM files from the alignment step.",
            default=False,
        ),
        "save_unaligned": NextflowParameter(
            type=bool,
            display_name="Save Unaligned Reads",
            description="Where possible, save unaligned reads from either STAR, HISAT2 or Salmon to the results directory.",
            default=False,
        ),
        "deseq2_vst": NextflowParameter(
            type=bool,
            display_name="DESeq2 VST",
            description="Use vst transformation instead of rlog with DESeq2.",
            default=True,
        ),
        "rseqc_modules": NextflowParameter(
            type=str,
            display_name="RSeQC Modules",
            description="Specify the RSeQC modules to run.",
            default="bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication",
        ),
        "skip_gtf_filter": NextflowParameter(
            type=bool,
            display_name="Skip GTF Filter",
            description="Skip filtering of GTF for valid scaffolds and/ or transcript IDs.",
            default=False,
        ),
        "skip_gtf_transcript_filter": NextflowParameter(
            type=bool,
            display_name="Skip GTF Transcript Filter",
            description="Skip the 'transcript_id' checking component of the GTF filtering script used in the pipeline.",
            default=False,
        ),
        "skip_bbsplit": NextflowParameter(
            type=bool,
            display_name="Skip BBSplit",
            description="Skip BBSplit for removal of non-reference genome reads.",
            default=True,
        ),
        "skip_umi_extract": NextflowParameter(
            type=bool,
            display_name="Skip UMI Extraction",
            description="Skip the UMI extraction from the read in case the UMIs have been moved to the headers in advance of the pipeline run.",
            default=False,
        ),
        "skip_trimming": NextflowParameter(
            type=bool,
            display_name="Skip Trimming",
            description="Skip the adapter trimming step.",
            default=False,
        ),
        "skip_alignment": NextflowParameter(
            type=bool,
            display_name="Skip Alignment",
            description="Skip all of the alignment-based processes within the pipeline.",
            default=False,
        ),
        "skip_pseudo_alignment": NextflowParameter(
            type=bool,
            display_name="Skip Pseudoalignment",
            description="Skip all of the pseudoalignment-based processes within the pipeline.",
            default=False,
        ),
        "skip_markduplicates": NextflowParameter(
            type=bool,
            display_name="Skip MarkDuplicates",
            description="Skip picard MarkDuplicates step.",
            default=False,
        ),
        "skip_bigwig": NextflowParameter(
            type=bool,
            display_name="Skip bigWig Creation",
            description="Skip bigWig file creation.",
            default=False,
        ),
        "skip_stringtie": NextflowParameter(
            type=bool,
            display_name="Skip StringTie",
            description="Skip StringTie.",
            default=False,
        ),
        "skip_fastqc": NextflowParameter(
            type=bool,
            display_name="Skip FastQC",
            description="Skip FastQC.",
            default=False,
        ),
        "skip_preseq": NextflowParameter(
            type=bool,
            display_name="Skip Preseq",
            description="Skip Preseq.",
            default=True,
        ),
        "skip_dupradar": NextflowParameter(
            type=bool,
            display_name="Skip dupRadar",
            description="Skip dupRadar.",
            default=False,
        ),
        "skip_qualimap": NextflowParameter(
            type=bool,
            display_name="Skip Qualimap",
            description="Skip Qualimap.",
            default=False,
        ),
        "skip_rseqc": NextflowParameter(
            type=bool,
            display_name="Skip RSeQC",
            description="Skip RSeQC.",
            default=False,
        ),
        "skip_biotype_qc": NextflowParameter(
            type=bool,
            display_name="Skip Biotype QC",
            description="Skip additional featureCounts process for biotype QC.",
            default=False,
        ),
        "skip_deseq2_qc": NextflowParameter(
            type=bool,
            display_name="Skip DESeq2 QC",
            description="Skip DESeq2 PCA and heatmap plotting.",
            default=False,
        ),
        "skip_multiqc": NextflowParameter(
            type=bool,
            display_name="Skip MultiQC",
            description="Skip MultiQC.",
            default=False,
        ),
        "skip_qc": NextflowParameter(
            type=bool,
            display_name="Skip QC",
            description="Skip all QC steps except for MultiQC.",
            default=False,
        ),
        "email": NextflowParameter(
            type=Optional[str],
            display_name="Email Address",
            description="Email address for completion summary.",
        ),
        "multiqc_title": NextflowParameter(
            type=Optional[str],
            display_name="MultiQC Report Title",
            description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
        ),
    },
    runtime_resources=NextflowRuntimeResources(cpus=4, memory=16, storage_gib=100),
    log_dir=LatchOutputDir("latch:///Bulk_RNAseq_logs"),
    flow=flow,
)
