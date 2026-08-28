import csv
import sys


def read_tsv(filename):
    data = {}
    columns = []

    with open(filename) as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = [c for c in reader.fieldnames if c != "sample"]

        for row in reader:
            sample = row.pop("sample")
            data[sample] = row

    return data, columns


qualimap, qualimap_cols = read_tsv(sys.argv[1])
picard, picard_cols = read_tsv(sys.argv[2])
fastqc, fastqc_cols = read_tsv(sys.argv[3])

samples = sorted(
    set(qualimap) |
    set(picard) |
    set(fastqc)
)

fieldnames = [
    "sample",
    "overall_qc",

    # PASS / FAIL
    "per_base_sequence_quality",
    "fold80_qc",
    "fold_enrichment_qc",
    "at_dropout_qc",
    "gc_dropout_qc",
    "duplicates_qc",
    "mapped_reads_percent_qc",
    "mapped_reads_inside_percent_qc",
    "mean_mapping_quality_qc",
    "coverage_20x_percent_qc",
    "mean_coverage_qc",

    # hodnoty
    "fold80",
    "fold_enrichment",
    "at_dropout",
    "gc_dropout",
    "duplicates",
    "mapped_reads_percent",
    "mapped_reads_inside_percent",
    "mean_mapping_quality",
    "coverage_20x_percent",
    "mean_coverage",
]

with open(sys.argv[4], "w", newline="") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fieldnames,
        delimiter="\t",
        extrasaction="ignore"
    )

    writer.writeheader()

    qc_columns = [
        "per_base_sequence_quality",
        "fold80_qc",
        "fold_enrichment_qc",
        "at_dropout_qc",
        "gc_dropout_qc",
        "duplicates_qc",
        "mapped_reads_percent_qc",
        "mapped_reads_inside_percent_qc",
        "mean_mapping_quality_qc",
        "coverage_20x_percent_qc",
        "mean_coverage_qc",
]

    for sample in samples:

        row = {"sample": sample}

        row.update(fastqc.get(sample, {}))
        row.update(picard.get(sample, {}))
        row.update(qualimap.get(sample, {}))

        qc_columns = [
            key for key, value in row.items()
            if key.endswith("_qc") or key == "per_base_sequence_quality"
]

        row["overall_qc"] = (
            "PASS"
            if all(value == "PASS" for key, value in row.items()
            if key.endswith("_qc") or key == "per_base_sequence_quality")
            else "FAIL"
        )

        writer.writerow(row)
