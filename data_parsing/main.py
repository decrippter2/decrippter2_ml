import json
from pathlib import Path

import parse_fastas
import parse_gbk

from data_parsing import MIBiG_parsing, access_ncbi_sequence


def main():
    folder = Path("/mibig_json_4.0rc1")
    mibig_dict = MIBiG_parsing.json_parse_ripps(folder)
    for key in mibig_dict:
        acc = mibig_dict[key]["accession"]
        start = mibig_dict[key]["start"]
        end = mibig_dict[key]["end"]
        parsed_prot_ids = mibig_dict[key]["protein_ids"]
        out_folder_gbk = "ncbi_gbk"
        print(acc, start, end)
        print("parsed protein ids: ", parsed_prot_ids)
        access_ncbi_sequence.download_gbk(acc, start, end, out_folder_gbk)
        if parsed_prot_ids == []:
            try:
                prot_lens, prot_ids = parse_gbk.parse_lines(
                    str(out_folder_gbk + "/" + acc + ".gbk"), mibig_dict[key]["ncomp"]
                )
                mibig_dict[key]["protein_ids"] = prot_ids
            except FileNotFoundError:
                mibig_dict[key]["protein_ids"] = []
                print(acc, "NOT FOUND ON NCBI")

            print("------------------")
        else:
            prot_ids = parsed_prot_ids
        for i in range(len(prot_ids)):
            try:
                if prot_ids[i] != None:
                    access_ncbi_sequence.download_prot(prot_ids[i], "prot_fastas")
                    mibig_dict[key][prot_ids[i]] = {}
                    mibig_dict[key][prot_ids[i]]["complete_seq"] = (
                        parse_fastas.parse_prot_seq(
                            f"prot_fastas/" + prot_ids[i] + ".fa"
                        )
                    )
                    try:
                        sequence, core_seq = (
                            mibig_dict[key][prot_ids[i]]["complete_seq"],
                            mibig_dict[key]["core_sequences"],
                        )
                    except IndexError:
                        sequence, core_seq = (
                            mibig_dict[key][prot_ids[i]]["complete_seq"],
                            "TBA",
                        )

                    lead, cor, fol = parse_fastas.classify_aa_seq(sequence, core_seq)
                    mibig_dict[key][prot_ids[i]]["leader_seq"] = lead
                    mibig_dict[key][prot_ids[i]]["core_seq"] = cor
                    mibig_dict[key][prot_ids[i]]["follower_seq"] = fol
                else:
                    prov_key = f"Unidentified protein{i}"
                    mibig_dict[key][prov_key] = {}
                    mibig_dict[key][prov_key]["complete_seq"] = "TBA"
                    mibig_dict[key][prov_key]["leader_seq"] = "TBA"
                    mibig_dict[key][prov_key]["core_seq"] = "TBA"
                    mibig_dict[key][prov_key]["follower_seq"] = "TBA"
            except FileNotFoundError:
                prov_key = f"Unidentified protein{i}"
                mibig_dict[key][prov_key] = {}
                mibig_dict[key][prov_key]["complete_seq"] = "TBA"
                mibig_dict[key][prov_key]["leader_seq"] = "TBA"
                mibig_dict[key][prov_key]["core_seq"] = "TBA"
                mibig_dict[key][prov_key]["follower_seq"] = "TBA"
    with open("repo.json", "w") as outfile:
        json.dump(mibig_dict, outfile, indent=6)


main()
