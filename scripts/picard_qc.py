from pathlib import Path
import csv


def parse_hs_metrics(file):
    lines = Path(file).read_text().splitlines()

    header = None
    values = None

    for i, line in enumerate(lines):
        if line.startswith("BAIT_SET"):
            header = line.split("\t")
            values = lines[i + 1].split("\t")
            break

    if header is None:
        raise ValueError(f"Header nebyl nalezen v {file}")

    return dict(zip(header, values))


def evaluate(metrics):

    return {
        "fold80_qc": float(metrics["FOLD_80_BASE_PENALTY"]) <= 2.0,
        "fold_enrichment_qc": float(metrics["FOLD_ENRICHMENT"]) >= 40,
        "at_dropout_qc": float(metrics["AT_DROPOUT"]) <= 0.50,
        "gc_dropout_qc": float(metrics["GC_DROPOUT"]) <= 10,
        "duplicates_qc": float(metrics["PCT_EXC_DUPE"]) <= 0.10,
    }


with open("picard_qc_summary.tsv", "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "sample",
        "fold80",
        "fold_enrichment",
        "at_dropout",
        "gc_dropout",
        "duplicates",
        "fold80_qc",
        "fold_enrichment_qc",
        "at_dropout_qc",
        "gc_dropout_qc",
        "duplicates_qc",
    ])

    for file in sorted(Path(".").glob("*hs_metrics.txt")):

        metrics = parse_hs_metrics(file)
        results = evaluate(metrics)

        # odstraní ".hs_metrics" z názvu
        sample = file.name.replace(".hs_metrics.txt", "")

        # pokud by soubory měly jinou koncovku
        sample = sample.replace(".hs_metrics", "")

        writer.writerow([
            sample,
            metrics["FOLD_80_BASE_PENALTY"],
            metrics["FOLD_ENRICHMENT"],
            metrics["AT_DROPOUT"],
            metrics["GC_DROPOUT"],
            metrics["PCT_EXC_DUPE"],
            "PASS" if results["fold80_qc"] else "FAIL",
            "PASS" if results["fold_enrichment_qc"] else "FAIL",
            "PASS" if results["at_dropout_qc"] else "FAIL",
            "PASS" if results["gc_dropout_qc"] else "FAIL",
            "PASS" if results["duplicates_qc"] else "FAIL",
        ])
