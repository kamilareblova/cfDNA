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
        "on_bait_bases_qc": float(metrics["ON_BAIT_BASES"]) <= 2.0,
        "off_bait_bases_qc": float(metrics["OFF_BAIT_BASES"]) >= 40,
        "pct_off_bait_qc": float(metrics["PCT_OFF_BAIT"]) <= 0.50,
        "pf_reads_qc": float(metrics["PF_READS"]) <= 10,
        "mean_target_coverage_qc": float(metrics["MEAN_TARGET_COVERAGE"]) <= 0.10,
        "median_target_coverage_qc": float(metrics["MEDIAN_TARGET_COVERAGE"]) <= 100,
        "zero_cvg_targets_pct_qc": float(metrics["ZERO_CVG_TARGETS_PCT"]) <= 100,
        "pct_target_bases_10X_qc": float(metrics["PCT_TARGET_BASES_10X"]) <= 100,
        "pct_target_bases_100X_qc": float(metrics["PCT_TARGET_BASES_100X"]) <= 100,
        "pct_target_bases_1000X_qc": float(metrics["PCT_TARGET_BASES_1000X"]) <= 100,
        "pct_target_bases_2500X_qc": float(metrics["PCT_TARGET_BASES_2500X"]) <= 100,  
        "pct_target_bases_5000X_qc": float(metrics["PCT_TARGET_BASES_5000X"]) <= 100
    }


with open("picard_qc_summary.tsv", "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "sample",
        "on_bait_bases",
        "off_bait_bases",
        "pct_off_bait",
        "pf_reads",
        "mean_target_coverage",
        "median_target_coverage",
        "zero_cvg_targets_pct",
        "pct_target_bases_10X",
        "pct_target_bases_100X",
        "pct_target_bases_1000X",
        "pct_target_bases_2500X",
        "pct_target_bases_5000X",
        "on_bait_bases_qc",
        "off_bait_bases_qc",
        "pct_off_bait_qc",
        "pf_reads_qc",
        "mean_target_coverage_qc",
        "median_target_coverage_qc",
        "zero_cvg_targets_pct_qc",
        "pct_target_bases_10X_qc",
        "pct_target_bases_100X_qc",
        "pct_target_bases_1000X_qc",
        "pct_target_bases_2500X_qc",
        "pct_target_bases_5000X_qc"
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
            metrics["ON_BAIT_BASES"],
            metrics["OFF_BAIT_BASES"],
            metrics["PCT_OFF_BAIT"],
            metrics["PF_READS"],
            metrics["MEAN_TARGET_COVERAGE"],
            metrics["MEDIAN_TARGET_COVERAGE"],
            metrics["ZERO_CVG_TARGETS_PCT"],
            metrics["PCT_TARGET_BASES_10X"],
            metrics["PCT_TARGET_BASES_100X"],
            metrics["PCT_TARGET_BASES_1000X"],
            metrics["PCT_TARGET_BASES_2500X"],
            metrics["PCT_TARGET_BASES_5000X"],
            "PASS" if results["on_bait_bases_qc"] else "FAIL",
            "PASS" if results["off_bait_bases_qc"] else "FAIL",
            "PASS" if results["pct_off_bait_qc"] else "FAIL",
            "PASS" if results["pf_reads_qc"] else "FAIL",
            "PASS" if results["mean_target_coverage_qc"] else "FAIL",
            "PASS" if results["median_target_coverage_qc"] else "FAIL",
            "PASS" if results["zero_cvg_targets_pct_qc"] else "FAIL",
            "PASS" if results["pct_target_bases_10X_qc"] else "FAIL",
            "PASS" if results["pct_target_bases_100X_qc"] else "FAIL",
            "PASS" if results["pct_target_bases_1000X_qc"] else "FAIL",
            "PASS" if results["pct_target_bases_2500X_qc"] else "FAIL",
            "PASS" if results["pct_target_bases_5000X_qc"] else "FAIL"
        ])
