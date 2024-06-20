import csv
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

from latch.ldata.path import LPath
from latch.resources.tasks import small_task
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile

sys.stdout.reconfigure(line_buffering=True)


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    strandedness: str
    differential_condition: Optional[str]


class DifferentialGeneTool(Enum):
    deseq2 = "deseq2"
    # sleuth = "sleuth"
    # edgeR = "edgeR"


def differential_sheet_constructor(samples: List[SampleSheet]) -> Path:
    samplesheet = Path("/root/dge_condition.csv")

    columns = ["sample", "condition"]

    with open(samplesheet, "w") as f:
        writer = csv.DictWriter(f, columns, delimiter=",")
        writer.writeheader()

        for sample in samples:
            row_data = {
                "sample": sample.sample,
                "condition": "Ignore"
                if sample.differential_condition is None
                else str(sample.differential_condition).lower(),
            }
            writer.writerow(row_data)

    return samplesheet


def detect_run_type(selected_run_path: LPath, tool: str):
    run_contents = [str(dir.path.split("/")[-1]) for dir in selected_run_path.iterdir()]

    print(run_contents)

    if "salmon" in run_contents:
        gene_count_object = LPath(
            f"{selected_run_path.path}/salmon/salmon.merged.gene_counts_length_scaled.rds"
        ).download()
        r_script = Path(f"/root/scripts/{tool}_salmon_kallisto.R")

    elif "kallisto" in run_contents:
        gene_count_object = LPath(
            f"{selected_run_path.path}/kallisto/kallisto.merged.gene_counts_length_scaled.rds"
        ).download()
        r_script = Path(f"/root/scripts/{tool}_salmon_kallisto.R")

    elif "star_salmon" in run_contents:
        gene_count_object = LPath(
            f"{selected_run_path.path}/star_salmon/salmon.merged.gene_counts_length_scaled.rds"
        ).download()
        r_script = Path(f"/root/scripts/{tool}_salmon_kallisto.R")

    elif "star_rsem" in run_contents:
        rsem_dir = LPath(f"{selected_run_path.path}/star_rsem")
        rsem_file_list = [
            file.download()
            for file in rsem_dir.iterdir()
            if file.path.endswith(".genes.results")
        ]
        r_script = Path(f"/root/scripts/{tool}_rsem.R")

        gene_count_object = "/root/rsem_file_list.txt"
        with open(gene_count_object, "w") as file_list:
            for file_path in rsem_file_list:
                file_list.write(f"{file_path}\n")

    elif "hisat2" in run_contents:
        print(
            "HISAT2 flow has no quantification outputs - https://nf-co.re/rnaseq/3.14.0"
        )
        gene_count_object = None
        r_script = None

    return gene_count_object, r_script


@small_task()
def dge(
    input: List[SampleSheet],
    run_name: str,
    run_latch_dge: Optional[DifferentialGeneTool],
    outdir: LatchOutputDir,
) -> Optional[LatchOutputDir]:
    if run_latch_dge:
        print("Setting up local directories")
        local_output_directory = Path(f"/root/output/{run_name}/{run_latch_dge.name}")
        local_output_directory.mkdir(parents=True, exist_ok=True)

        gene_count_object, r_script = detect_run_type(
            selected_run_path=LPath(f"{outdir.remote_path}/{run_name}"),
            tool=run_latch_dge.name,
        )
        condition_sheet = differential_sheet_constructor(samples=input)

        if gene_count_object and r_script:
            print("Running Deseq2")
            deseq_cmd = [
                "Rscript",
                str(r_script),
                str(gene_count_object),
                str(condition_sheet),
                str(local_output_directory),
            ]
            subprocess.run(deseq_cmd, check=True, cwd=local_output_directory)

            print("Uploading results")
            return LatchOutputDir(str("/root/output"), outdir.remote_path)

    return None
