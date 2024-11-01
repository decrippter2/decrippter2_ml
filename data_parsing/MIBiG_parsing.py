import json
import os
import re


def detect_gene_and_core(bgc_dict, pos):
    core_seqs = []
    prot_id_pat = r"[a-z]+"
    protein_ids = []
    for gene in bgc_dict["biosynthesis"]["classes"][pos]["precursors"]:
        if type(gene["core_sequence"]) == list:
            core_seqs.append(gene["core_sequence"][0])
        elif type(gene["core_sequence"]) == str:
            core_seqs.append(gene["core_sequence"].upper())
        if not re.match(prot_id_pat, gene["gene"]):
            protein_ids.append(gene["gene"])
        else:
            protein_ids = []
    return protein_ids, core_seqs


def detect_ripp_core_id(mibig_dict):
    for i in range(len(mibig_dict["biosynthesis"]["classes"])):
        if mibig_dict["biosynthesis"]["classes"][i]["class"] == "ribosomal":
            flag = True
            protein_ids, core_seqs = detect_gene_and_core(mibig_dict, i)
    try:
        print(flag)
    except NameError:
        flag = False
        core_seqs = ["None"]
    return flag, core_seqs, protein_ids


def json_ripp_access(bgc_dict):
    access = bgc_dict["loci"][0]["accession"]
    n_comp = len(bgc_dict["compounds"])
    try:
        start = bgc_dict["loci"][0]["location"]["from"]
        end = bgc_dict["loci"][0]["location"]["to"]
    except KeyError:
        start = None
        end = None
    return access, n_comp, start, end


def json_parse_ripps(json_folder):
    ripp_dict = {}
    for _dirpath, _dirnames, filenames in os.walk(json_folder):
        for file in filenames:
            filename = str(json_folder) + "/" + file
            with open(filename) as raw_json:
                bgc_json = raw_json.read()
            bgc_dict = json.loads(bgc_json)
            try:
                for elem in bgc_dict["biosynthesis"]["classes"][0]["precursors"]:
                    print("accession: ", bgc_dict["loci"][0]["accession"])
                    print("gene entry: ", elem["gene"])
            except KeyError:
                continue
            flag, core_seqs, prec_id = detect_ripp_core_id(bgc_dict)
            if flag is True:
                access, ncomp, start, end = json_ripp_access(bgc_dict)
                ripp_dict[file] = {}
                ripp_dict[file]["accession"] = access
                ripp_dict[file]["ncomp"] = ncomp
                ripp_dict[file]["start"] = start
                ripp_dict[file]["end"] = end
                ripp_dict[file]["core_sequences"] = core_seqs
                ripp_dict[file]["protein_ids"] = []
                if prec_id != []:
                    ripp_dict[file]["protein_ids"] = prec_id
    print(ripp_dict)
    return ripp_dict
