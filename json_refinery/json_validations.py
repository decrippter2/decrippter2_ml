import json
import os


def import_dict(json_file):
    with open(json_file) as in_json:
        json_dict = json.load(in_json)
    return json_dict


def check_sequence_info(json_dict, filename):
    ids = json_dict["mibig_acc"]["protein_ids"]
    for id in ids:
        try:
            protein_dict = json_dict["mibig_acc"][id]
            concatenation = (
                protein_dict["leader_seq"]
                + protein_dict["core_seq"]
                + protein_dict["follower_seq"]
            )
            if protein_dict["complete_seq"] == concatenation:
                continue
            else:
                print(filename, " contains WRONG CONCATENATION")
        except KeyError:
            print("missing protein entry")


def check_cores(json_dict, filename):
    ids = json_dict["mibig_acc"]["protein_ids"]
    cores = json_dict["mibig_acc"]["core_sequences"]
    extracted_cores = []
    for id in ids:
        try:
            id_core = json_dict["mibig_acc"][id]["core_seq"]
            extracted_cores.append(id_core)
        except KeyError:
            print("missing protein entry")
    for core in cores:
        if core in extracted_cores:
            continue
        else:
            print(filename, " is missing a proper core in the breakdown")
    for core in extracted_cores:
        if core in cores:
            continue
        else:
            print(filename, " is missing a core in the corelist")


def check_prot_ids(json_dict, filename):
    ids = json_dict["mibig_acc"]["protein_ids"]
    for id in ids:
        if id in json_dict["mibig_acc"]:
            continue
        else:
            print(filename, " is missing a proper protein id")


def check_accession(json_dict, filename):
    acc = json_dict["mibig_acc"]["accession"]
    if acc.startswith("WP"):
        print(filename, " does not have proper genome accession")


def check_emptiness(json_dict, filename):
    ids = json_dict["mibig_acc"]["protein_ids"]
    for id in ids:
        try:
            if "TBA" in json_dict["mibig_acc"][id].values():
                print(filename, " has TBA")
            else:
                continue
        except KeyError:
            print(KeyError)
    if json_dict["mibig_acc"]["ref"] == "":
        print(filename, " is missing reference")


def check_class(json_dict, filename):
    class_info = json_dict["mibig_acc"]["ripp_class"]
    if class_info == "TBA" or class_info == "":
        print(filename, " is missing RiPP class")


def find_class_members(json_dict, filename, class_name):
    if json_dict["mibig_acc"]["ripp_class"] == class_name:
        print(filename, " presents class: ", class_name)


def main():
    count_core = 0
    count_prots = 0
    for _dirpath, _dirnames, filenames in os.walk(
        "C:/Users/rsanz/PycharmProjects/decrippter2_ml/json_data"
    ):
        for filename in filenames:
            print("------")
            print("file: " + filename)
            file = (
                "C:/Users/rsanz/PycharmProjects/decrippter2_ml/json_data"
                + "/"
                + filename
            )
            dictionary = import_dict(file)
            check_emptiness(dictionary, file)
            check_accession(dictionary, file)
            check_prot_ids(dictionary, file)
            check_sequence_info(dictionary, file)
            check_class(dictionary, file)
            check_cores(dictionary, file)
            count_core += len(dictionary["mibig_acc"]["core_sequences"])
            count_prots += len(dictionary["mibig_acc"]["protein_ids"])

    print(count_core, " core sequences")
    print(count_prots, " proteins")


main()
