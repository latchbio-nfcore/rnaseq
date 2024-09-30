from dataclasses import dataclass
from enum import Enum
from typing import Optional

from latch.types.file import LatchFile


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
    rhesus_macaque = "Macaca mulatta (RefSeq rheMac10/Mmul_10)"
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
