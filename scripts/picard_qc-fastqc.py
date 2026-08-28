from pathlib import Path
import zipfile
import csv


def parse_fastqc(zip_file):
    with zipfile.ZipFile(zip_file) as z:

        summary = next(x for x in z.namelist() if x.endswith("summary.txt"))

        with z.open(summary) as f:

            for line in f:
                status, module, sample = line.decode().strip().split("\t")

                if module == "Per base sequence quality":
                    return status

    raise ValueError(f"'Per base sequence quality' nebyl nalezen v {zip_file}")


# uloží statusy R1/R2 pro každý vzorek
samples = {}

for zip_file in Path(".").glob("*_fastqc.zip"):

    name = zip_file.stem.replace("_fastqc", "")

    if name.endswith("_R1"):
        sample = name[:-3]
        read = "R1"
    elif name.endswith("_R2"):
        sample = name[:-3]
        read = "R2"
    else:
        continue

    samples.setdefault(sample, {})[read] = parse_fastqc(zip_file)


with open("fastqc_summary.tsv", "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "sample",
        "per_base_sequence_quality",
    ])

    for sample in sorted(samples):

        r1 = samples[sample].get("R1")
        r2 = samples[sample].get("R2")

        if r1 == "PASS" and r2 == "PASS":
            qc = "PASS"
        else:
            qc = "FAIL"

        writer.writerow([sample, qc])
