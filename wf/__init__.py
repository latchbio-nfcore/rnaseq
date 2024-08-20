import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
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

from wf.entrypoint import custom_samplesheet_constructor, initialize, nextflow_runtime
from wf.prep_dge import prep_dge

# Define the structure of the workflow UI
flow = [
    Section(
        "Samples",
        Params(
            "input",
        ),
        Text(
            "Sample identifier and FASTQ files should not contain spaces in file names or full directory locations."
        ),
        Text(
            "Strandedness can be set to 'auto', 'reverse', 'forward'. If left untoggled, it will default to 'auto'."
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
                "tx2gene_file",
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


# Define Nextflow metadata for the workflow
NextflowMetadata(
    display_name="nf-core/rnaseq",
    documentation="https://wiki.latch.bio/workflows/bulk-rna-seq#bulk-rnaseq-quantification",
    wiki_url="https://wiki.latch.bio/workflows/bulk-rna-seq#bulk-rnaseq-quantification",
    author=LatchAuthor(
        name="nf-core",
        github="https://github.com/latchbio-nfcore/rnaseq",
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
        "tx2gene_file": NextflowParameter(
            type=Optional[LatchFile],
            display_name="Tx2gene TSV File",
            description="Tx2Gene file for tximport grouping before Deseq file generation for RSEM runs.",
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
    about_page_path=Path("./README.md"),
)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_rnaseq(
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
    tx2gene_file: Optional[LatchFile],
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
    """nf-core/rnaseq is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">
    Produce transcript/count matrices from sequencing reads.

    nf-core/rnaseq is an open-source bioinformatics workflow that processes raw sequencing reads, aligns them to genes, and performs quality control checks. The pipeline uses gold-standard tools, is maintained by the growing nf-core community and can be modified or extended as needed to adapt to bespoke biology.

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

    The managed computing infrastructure scales to hundreds of samples, with clear logging and error-reporting. Data provenance links versioned and containerized workflow code to input and output files.

    </p>

    </html>

    # nf-core/rnaseq

    [![GitHub Actions CI Status](https://github.com/nf-core/rnaseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/rnaseq/actions?query=workflow%3A%22nf-core+CI%22)
    [![GitHub Actions Linting Status](https://github.com/nf-core/rnaseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/rnaseq/actions?query=workflow%3A%22nf-core+linting%22)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/rnaseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.1400710-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.1400710)

    [![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)

    [![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rnaseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/rnaseq)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

    ## Introduction

    **nf-core/rnaseq** is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.

    ![nf-core/rnaseq metro map](docs/images/nf-core-rnaseq_metro_map_grey.png)

    1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
    2. Sub-sample FastQ files and auto-infer strandedness ([`fq`](https://github.com/stjude-rust-labs/fq), [`Salmon`](https://combine-lab.github.io/salmon/))
    3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    4. UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
    5. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
    6. Removal of genome contaminants ([`BBSplit`](http://seqanswers.com/forums/showthread.php?t=41288))
    7. Removal of ribosomal RNA ([`SortMeRNA`](https://github.com/biocore/sortmerna))
    8. Choice of multiple alignment and quantification routes:
    1. [`STAR`](https://github.com/alexdobin/STAR) -> [`Salmon`](https://combine-lab.github.io/salmon/)
    2. [`STAR`](https://github.com/alexdobin/STAR) -> [`RSEM`](https://github.com/deweylab/RSEM)
    3. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) -> **NO QUANTIFICATION**
    9. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    10. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
    11. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
    12. Transcript assembly and quantification ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
    13. Create bigWig coverage files ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    14. Extensive quality control:
        1. [`RSeQC`](http://rseqc.sourceforge.net/)
        2. [`Qualimap`](http://qualimap.bioinfo.cipf.es/)
        3. [`dupRadar`](https://bioconductor.org/packages/release/bioc/html/dupRadar.html)
        4. [`Preseq`](http://smithlabresearch.org/software/preseq/)
        5. [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
    15. Pseudoalignment and quantification ([`Salmon`](https://combine-lab.github.io/salmon/) or ['Kallisto'](https://pachterlab.github.io/kallisto/); _optional_)
    16. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

    > **Note**
    > The SRA download functionality has been removed from the pipeline (`>=3.2`) and ported to an independent workflow called [nf-core/fetchngs](https://nf-co.re/fetchngs). You can provide `--nf_core_pipeline rnaseq` when running nf-core/fetchngs to download and auto-create a samplesheet containing publicly available samples that can be accepted directly as input by this pipeline.

    > **Warning**
    > Quantification isn't performed if using `--aligner hisat2` due to the lack of an appropriate option to calculate accurate expression estimates from HISAT2 derived genomic alignments. However, you can use this route if you have a preference for the alignment, QC and other types of downstream analysis compatible with the output of HISAT2.

    ## Usage

    First, prepare a samplesheet with your input data that looks as follows:

    **samplesheet.csv**:

    sample,fastq_1,fastq_2,strandedness
    CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto
    CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,auto
    CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,auto

    Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Rows with the same sample identifier are considered technical replicates and merged automatically. The strandedness refers to the library preparation and will be automatically inferred if set to `auto`.

    > **Warning:**
    > Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
    > provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
    > see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

    For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/rnaseq/usage) and the [parameter documentation](https://nf-co.re/rnaseq/parameters).

    ## Pipeline output

    To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/rnaseq/results) tab on the nf-core website pipeline page.
    For more details about the output files and reports, please refer to the
    [output documentation](https://nf-co.re/rnaseq/output).

    This pipeline quantifies RNA-sequenced reads relative to genes/transcripts in the genome and normalizes the resulting data. It does not compare the samples statistically in order to assign significance in the form of FDR or P-values. For downstream analyses, the output files from this pipeline can be analysed directly in statistical environments like [R](https://www.r-project.org/), [Julia](https://julialang.org/) or via the [nf-core/differentialabundance](https://github.com/nf-core/differentialabundance/) pipeline.

    ## Online videos

    A short talk about the history, current status and functionality on offer in this pipeline was given by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) on [8th February 2022](https://nf-co.re/events/2022/bytesize-32-nf-core-rnaseq) as part of the nf-core/bytesize series.

    You can find numerous talks on the [nf-core events page](https://nf-co.re/events) from various topics including writing pipelines/modules in Nextflow DSL2, using nf-core tooling, running nf-core pipelines as well as more generic content like contributing to Github. Please check them out!

    ## Credits

    These scripts were originally written for use at the [National Genomics Infrastructure](https://ngisweden.scilifelab.se), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

    The pipeline was re-written in Nextflow DSL2 and is primarily maintained by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [Seqera Labs, Spain](https://seqera.io/).

    The pipeline workflow diagram was initially designed by Sarah Guinchard ([@G-Sarah](https://github.com/G-Sarah)) and James Fellows Yates ([@jfy133](https://github.com/jfy133)), further modifications where made by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) and Maxime Garcia ([@maxulysse](https://github.com/maxulysse)).

    Many thanks to other who have helped out along the way too, including (but not limited to):

    - [Alex Peltzer](https://github.com/apeltzer)
    - [Colin Davenport](https://github.com/colindaven)
    - [Denis Moreno](https://github.com/Galithil)
    - [Edmund Miller](https://github.com/Emiller88)
    - [Gregor Sturm](https://github.com/grst)
    - [Jacki Buros Novik](https://github.com/jburos)
    - [Lorena Pantano](https://github.com/lpantano)
    - [Matthias Zepper](https://github.com/MatthiasZepper)
    - [Maxime Garcia](https://github.com/maxulysse)
    - [Olga Botvinnik](https://github.com/olgabot)
    - [@orzechoj](https://github.com/orzechoj)
    - [Paolo Di Tommaso](https://github.com/pditommaso)
    - [Rob Syme](https://github.com/robsyme)

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on the [Slack `#rnaseq` channel](https://nfcore.slack.com/channels/rnaseq) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    If you use nf-core/rnaseq for your analysis, please cite it using the following doi: [10.5281/zenodo.1400710](https://doi.org/10.5281/zenodo.1400710)

    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

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
        tx2gene_file=tx2gene_file,
        outdir=outdir,
    )


LaunchPlan(
    nf_nf_core_rnaseq,
    "Test Data",
    {
        "input": [
            SampleSheet(
                sample="WT_REP1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/rnaseq/nfcore_test/SRR6357070_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/rnaseq/nfcore_test/SRR6357070_2.fastq.gz"
                ),
                strandedness="auto",
            ),
            SampleSheet(
                sample="WT_REP2",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/rnaseq/nfcore_test/SRR6357072_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/rnaseq/nfcore_test/SRR6357072_2.fastq.gz"
                ),
                strandedness="reverse",
            ),
            SampleSheet(
                sample="RAP1_UNINDUCED_REP1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/rnaseq/nfcore_test/SRR6357073_1.fastq.gz"
                ),
                fastq_2=None,
                strandedness="reverse",
            ),
            SampleSheet(
                sample="RAP1_IAA_30M_REP1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/rnaseq/nfcore_test/SRR6357076_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/rnaseq/nfcore_test/SRR6357076_2.fastq.gz"
                ),
                strandedness="reverse",
            ),
        ],
        "genome_source": "custom",
        "fasta": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/genome.fasta",
        "gtf": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/genes_with_empty_tid.gtf.gz",
        "gff": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/genes.gff.gz",
        "transcript_fasta": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/transcriptome.fasta",
        "additional_fasta": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/gfp.fa.gz",
        "hisat2_index": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/hisat2.tar.gz",
        "rsem_index": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/rsem.tar.gz",
        "salmon_index": "s3://latch-public/nf-core/rnaseq/nfcore_test/reference/salmon.tar.gz",
        "run_name": "Test",
    },
)
