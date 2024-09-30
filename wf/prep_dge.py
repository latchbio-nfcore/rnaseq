import subprocess
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
from latch.functions.messages import message
from latch.ldata.path import LPath
from latch.resources.tasks import small_task
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile

from wf.dataclasses import (
    Reference_Type,
)

# Configure stdout for line buffering
sys.stdout.reconfigure(line_buffering=True)


def detect_method_type(selected_run_path: LPath):
    """
    Detects the quantification method used in the RNA-seq analysis and locates the gene count file.

    Args:
        selected_run_path (LPath): The path to the analysis run directory.

    Returns:
        tuple: A tuple containing the path to the gene count file and the detected method.
    """
    run_contents = [str(dir.path.split("/")[-1]) for dir in selected_run_path.iterdir()]

    print(run_contents)

    if "salmon" in run_contents:
        gene_count_object = LPath(
            f"{selected_run_path.path}/salmon/salmon.merged.gene_counts_length_scaled.tsv"
        ).download()
        method = "salmon"

    elif "kallisto" in run_contents:
        gene_count_object = LPath(
            f"{selected_run_path.path}/kallisto/kallisto.merged.gene_counts_length_scaled.tsv"
        ).download()
        method = "kallisto"

    elif "star_salmon" in run_contents:
        gene_count_object = LPath(
            f"{selected_run_path.path}/star_salmon/salmon.merged.gene_counts_length_scaled.tsv"
        ).download()
        method = "star_salmon"

    elif "star_rsem" in run_contents:
        rsem_dir = LPath(f"{selected_run_path.path}/star_rsem")
        rsem_file_list = [
            file.download()
            for file in rsem_dir.iterdir()
            if file.path.endswith(".isoforms.results")
        ]
        method = "rsem"

        gene_count_object = "/root/rsem_file_list.txt"
        with open(gene_count_object, "w") as file_list:
            for file_path in rsem_file_list:
                file_list.write(f"{file_path}\n")

    elif "hisat2" in run_contents:
        print(
            "HISAT2 flow has no quantification outputs - https://nf-co.re/rnaseq/3.14.0"
        )
        gene_count_object = None
        method = None

    return gene_count_object, method


@small_task
def prep_dge(
    run_name: str,
    latch_genome: Reference_Type,
    tx2gene_file: Optional[LatchFile],
    outdir: LatchOutputDir,
) -> LatchOutputDir:
    """
    Prepares differential gene expression (DGE) analysis input files.

    This function processes the gene count data from various quantification methods
    and prepares it for use in differential expression analysis tools like DESeq2.

    Args:
        run_name (str): The name of the analysis run.
        latch_genome (Reference_Type): The reference genome used in the analysis.
        tx2gene_file (Optional[LatchFile]): A file mapping transcripts to genes (optional).
        outdir (LatchOutputDir): The output directory for the analysis results.

    Returns:
        LatchOutputDir: The directory containing the prepared DGE input files.
    """
    print("Setting up local directories")
    local_output_directory = Path(f"/root/output/{run_name}/deseq2_counts")
    local_output_directory.mkdir(parents=True, exist_ok=True)

    gene_count_object, method = detect_method_type(
        selected_run_path=LPath(f"{outdir.remote_path}/{run_name}"),
    )

    if method == "rsem":
        try:
            if tx2gene_file:
                tx2gene_file_p = Path(tx2gene_file)
            else:
                tx2gene_file_p = Path(
                    f"/root/assets/latch_tx2gene/tx2gene_{latch_genome.name}.tsv"
                )

            prep_cmd = [
                "Rscript",
                "/root/latch_scripts/prep_deseq2_rsem.R",
                str(gene_count_object),
                str(tx2gene_file_p),
                str(local_output_directory),
            ]
            subprocess.run(prep_cmd, check=True, cwd=local_output_directory)
        except Exception as e:
            message(
                "error",
                {
                    "title": "Could not make Deseq2 input file for RSEM run",
                    "body": f"Error: {str(e)}",
                },
            )

    elif method is not None and gene_count_object is not None:
        print(f"Processing {method} gene counts")
        df = pd.read_csv(gene_count_object, sep="\t")

        if "gene_name" in df.columns:
            df = df.drop(columns=["gene_name"])

        if df.columns[0] != "gene_id":
            df = df.rename(columns={df.columns[0]: "gene_id"})

        output_file = (
            local_output_directory / f"deseq2_{method}_length_scaled_counts.tsv"
        )
        df.to_csv(output_file, index=False, sep="\t")

        print(f"Wrote DESeq2 input file to {output_file}")

    print("Uploading results")
    return LatchOutputDir(str("/root/output"), outdir.remote_path)
