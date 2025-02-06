import json
import os


def import_dict(json_file):
    with open(json_file) as in_json:
        json_dict = json.load(in_json)
    return json_dict


for _dirpath, _dirnames, filenames in os.walk(
        "/json_data"
):
    for filename in filenames:
        file = (
            "C:/Users/rsanz/PycharmProjects/decrippter2_ml/json_data"
            + "/"
            + filename
        )
        dictionary = import_dict(file)
        bgc_number = dictionary["identifiers"]["mibig_acc"]
        for entry in dictionary["entries"]:
            if entry["leader"] == [] and entry["follower"] == []:
                print("missing leader/follower in ", filename)
            for core in entry["core"]:
                if core in entry["leader"][0] and core != "TBA":
                    print(
                        "Core repetion in leader sequence of ",
                        filename[:-5],
                        " from BGC ",
                        bgc_number,
                    )
                elif core in entry["follower"][0] and core != "TBA":
                    print(
                        "Core repetion in follower sequence of ",
                        filename[:-5],
                        " from BGC ",
                        bgc_number,
                    )
