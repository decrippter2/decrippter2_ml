import csv


def parse_csv(filepath):
    access_list = []
    with open(filepath, newline="") as smash_csv:
        reader = csv.DictReader(smash_csv, delimiter="\t")
        for row in reader:
            info = [row["NCBI accession"], row["From"], row["To"]]
            access_list.append(info)
    return access_list
