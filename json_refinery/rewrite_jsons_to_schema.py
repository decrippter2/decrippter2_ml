import json
import os


def write_json_dict(infolder, outfolder):
    for _dirpath, _dirnames, filenames in os.walk(infolder):
        file_count = 1
        for filename in filenames:
            out_dict = {
                "identifiers": {"ncbi_acc": "", "mibig_acc": ""},
                "compound_name": [],
                "ripp_class": "",
                "ref": [],
                "entries": [],
            }

            with open(infolder + "/" + filename) as json_file:
                json_dict = json.load(json_file)
                json_name_split = filename.split("/")
                json_name = json_name_split[-1]
                out_dict["identifiers"]["ncbi_acc"] = json_dict["mibig_acc"][
                    "accession"
                ]
                out_dict["identifiers"]["mibig_acc"] = json_name[:-5]
                for compound in json_dict["mibig_acc"]["compound_name"]:
                    out_dict["compound_name"].append(compound)
                out_dict["ripp_class"] = json_dict["mibig_acc"]["ripp_class"]
                out_dict["ref"].append(json_dict["mibig_acc"]["ref"])

                for protein_id in json_dict["mibig_acc"]["protein_ids"]:
                    entry_dict = {
                        "protein_ids": {"genpept": "", "uniprot": ""},
                        "complete": "",
                        "leader": [],
                        "core": [],
                        "follower": [],
                    }
                    entry_dict["protein_ids"]["genpept"] = protein_id
                    entry_dict["complete"] = json_dict["mibig_acc"][protein_id][
                        "complete_seq"
                    ]
                    # check for entry duplicates
                    """
                    exists_dupl=False
                    for key in json_dict["mibig_acc"]:
                        if protein_id[:-2] in key:
                            if "_" in key or ".2" in key:
                                exists_dupl=True
                        else:
                            continue
                    if exists_dupl==True:
                        lead_list=[]
                        fol_list=[]
                        for key in json_dict["mibig_acc"]:
                            if protein_id[:-2] in key:
                                lead_list.append(json_dict["mibig_acc"][key]["leader_seq"])
                        out_dict["entries"].append(entry_dict)
                    else:
                    """
                    lead_list = []
                    cor_list = []
                    fol_list = []
                    for key in json_dict["mibig_acc"]:
                        if protein_id[:-2] in key:
                            lead_list.append(json_dict["mibig_acc"][key]["leader_seq"])
                            fol_list.append(json_dict["mibig_acc"][key]["follower_seq"])
                            if json_dict["mibig_acc"][key]["core_seq"] not in cor_list:
                                cor_list.append(json_dict["mibig_acc"][key]["core_seq"])
                    shortest_lead = min(lead_list, key=len)
                    shortest_follower = min(fol_list, key=len)

                    entry_dict["leader"].append(shortest_lead)
                    entry_dict["follower"].append(shortest_follower)
                    entry_dict["core"].extend(cor_list)
                    out_dict["entries"].append(entry_dict)
            outfile_name = "ripp" + "0" * (7 - len(str(file_count))) + str(file_count)
            with open(outfolder + "/" + outfile_name + ".json", "w") as outfile:
                json.dump(out_dict, outfile, indent=6)
            file_count += 1


write_json_dict(
    "C:/Users/rsanz/PycharmProjects/decrippter2_ml/json_data",
    "C:/Users/rsanz/PycharmProjects/decrippter2_ml/json_data",
)
