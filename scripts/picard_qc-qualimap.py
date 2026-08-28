from pathlib import Path
import re


def parse_genome_results(file):

    metrics = {}
    section = None
   
    lines = Path(file).read_text().splitlines()


    for line in lines:

        line = line.strip()

        # zjisti, ve které sekci jsme
        if line.startswith(">>>>>>> Globals inside"):
            section = "globals_inside"
            continue

        elif line.startswith(">>>>>>> Globals"):
            section = "globals"
            continue


        if "number of mapped reads" in line:

            # číslo v závorce (99.79%)
            m = re.search(r"\(([\d.]+)%\)", line)

            if m:
                
                if section =="globals":
                    metrics["mapped_reads_percent"] = float(m.group(1))

                elif section == "globals_inside":
                    metrics["mapped_reads_inside_percent"] = float(m.group(1))

        # mean mapping quality = 58.1376
        if "mean mapping quality" in line:
            m = re.search(r"=\s*([\d.]+)", line)
            if m:
                metrics["mean_mapping_quality"] = float(m.group(1))
        
                # mean mapping quality = 58.1376
        if "mean coverageData" in line:
            m = re.search(r"=\s*([\d.]+)", line)
            if m:
                metrics["mean_coverage"] = float(m.group(1))

        # coverage_20x
        if "reference with a coverageData >= 20X" in line:
            m = re.search(r"(\d+\.\d+)%", line)
            if m:
                metrics["coverage_20x_percent"] = float(m.group(1))

    return metrics


# projdi všechny genome_metrics soubory

for folder in Path(".").glob("*.qualimap"):

	genome_file = folder / "genome_results.txt"
 
	if genome_file.exists():

		metrics = parse_genome_results(genome_file)

		print(folder.name)
		print(f"Mapped reads (all):   {metrics['mapped_reads_percent']} %")
		print(f"Mapped reads (inside):   {metrics['mapped_reads_inside_percent']} %")
		print(f"Mapping quality:    {metrics['mean_mapping_quality' ]} ")
		print(f"Mean coverage:    {metrics['mean_coverage' ]} ")
		print(f"Coverage 20x: {metrics['coverage_20x_percent']} %")
		print()

def evaluate(metrics):

    return {
        "mapped_reads_percent_qc": float(metrics["mapped_reads_percent"]) >= 98,
        "mapped_reads_inside_percent_qc": float(metrics["mapped_reads_inside_percent"]) >= 70,
        "mean_mapping_quality_qc": float(metrics["mean_mapping_quality"]) >= 50,
        "coverage_20x_percent_qc": float(metrics["coverage_20x_percent"]) >= 95,
        "mean_coverage_qc": float(metrics["mean_coverage"]) >= 100,
    }


results = evaluate(metrics)

for metric, ok in results.items():
    status = "PASS" if ok else "FAIL"
    print(f"{metric:15} {status}")

### zapise do csv
from pathlib import Path
import csv

with open("qualimap_qc_summary.tsv", "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "sample",
        "mapped_reads_percent",
        "mapped_reads_inside_percent",
        "mean_mapping_quality",
        "coverage_20x_percent",
        "mean_coverage",
        "mapped_reads_percent_qc",
        "mapped_reads_inside_percent_qc",
        "mean_mapping_quality_qc",
        "coverage_20x_percent_qc",
        "mean_coverage_qc",

    ])

    for folder in Path(".").glob("*.qualimap"):

        genome_file = folder / "genome_results.txt"

        if not genome_file.exists():

            continue

        metrics = parse_genome_results(genome_file)

        results = evaluate(metrics)

        writer.writerow([
            folder.stem.replace(".qualimap", ""),
            metrics["mapped_reads_percent"],
            metrics["mapped_reads_inside_percent"],
            metrics["mean_mapping_quality"],
            metrics["coverage_20x_percent"],
            metrics["mean_coverage"],
            "PASS" if results["mapped_reads_percent_qc"] else "FAIL",
            "PASS" if results["mapped_reads_inside_percent_qc"] else "FAIL",
            "PASS" if results["mean_mapping_quality_qc"] else "FAIL",
            "PASS" if results["coverage_20x_percent_qc"] else "FAIL",
            "PASS" if results["mean_coverage_qc"] else "FAIL",
        ])
