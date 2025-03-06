import json
import os

import matplotlib.pyplot as plt

abbrev_dict = {
    "AIP": "AIP",
    "Atropopeptide": "Atro",
    "Bacteriocin": "Btcn",
    "Biarylitide": "Bryl",
    "Borosin": "Boro",
    "Bottromycin": "Bttr",
    "Cittilin": "Citt",
    "Crocagin": "Croc",
    "Cyanobactin": "Cyan",
    "Cyptide": "Cypt",
    "Dikaritin": "Dikr",
    "Glycocin": "Glyc",
    "Graspetide": "Gras",
    "Guanidinotide": "Guan",
    "Head-to-tail cyclized peptide": "HtTC",
    "Lanthipeptide": "Lant",
    "Lasso": "Lass",
    "Linaridin": "Linr",
    "Linear azole-containing peptide": "LAP",
    "Methanobactin": "Meth",
    "Microcin": "Micr",
    "Microviridin": "Micv",
    "Mycofactocin": "Mycf",
    "N-Formylated TBA": "NFor",
    "Pantocin": "Pant",
    "Pearlin": "Prln",
    "Proteusin": "Prou",
    "Ranthipeptide": "Rant",
    "Sactipeptide": "Sact",
    "Selidamide": "Seli",
    "Streptide": "Strp",
    "Sulfatyrotide": "Sulf",
    "TBA": "TBA",
    "Thioamitide": "Thmd",
    "Thiopeptide": "Thpt",
    "Ustiloxin": "Uslx",
    "other": "Othr",
    "Amatoxin": "Amtx",
    "Cyclotide": "Cycl",
    "Lyciumin": "Lycm",
    "Rotapeptide": "Rtpp",
    "Ryptide": "Rypt",
    "Epipeptide": "Epip",
}


def show_class_distribution(class_dict):
    key_list = []
    value_list = []
    x_coord = []
    class_dict = dict(
        sorted(class_dict.items(), key=lambda item: item[1], reverse=True)
    )
    for key in class_dict:
        key_list.append(key)
        value_list.append(class_dict[key])
    for i in range(len(key_list)):
        x_coord.append(i)
    print(sorted(key_list))
    print(sum(value_list))

    fig, ax = plt.subplots()
    ax.barh(x_coord, value_list, align="center")
    ax.set_yticks(x_coord, labels=key_list)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel("N. proteins")
    ax.set_title("Family distribution")
    fig.set_figheight(10)
    for index, value in enumerate(value_list):
        plt.text(value, index, str(value))
    plt.savefig("graph.png")


def import_dict(json_file):
    with open(json_file) as in_json:
        json_dict = json.load(in_json)
    return json_dict


no_amplified_class_dictionary = {}
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
        class_abrev = abbrev_dict[dictionary["ripp_class"]]
        if class_abrev not in no_amplified_class_dictionary:
            no_amplified_class_dictionary[class_abrev] = len(dictionary["entries"])
        elif class_abrev in no_amplified_class_dictionary:
            no_amplified_class_dictionary[class_abrev] += len(dictionary["entries"])


show_class_distribution(noo_amplified_class_dictionary)
