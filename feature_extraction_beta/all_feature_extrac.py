import contextlib
import json
from sys import argv

import peptides
from pyPept import conformer
from numpy import log2
from feature_extraction import RiPP

def alpha_helical_propensity(sequence):
    predictor = conformer.SecStructPredictor()
    sec_str = predictor.predict_active_ss(sequence)
    alpha_helix = (sec_str.count('H')) / len(sec_str)
    return alpha_helix

def write_dataset(fasta_file, outfile_name, is_ripp):
    with open(
        "/lustre/BIF/nobackup/sanz006/positive_training/entry_category.json"
    ) as json_cat:
        json_category_dict = json.load(json_cat)
    with open(outfile_name, "w") as outfile:
        header_list = [
            "protein_id",
            "sequence",
            "RiPP",
            "Class",
            "validation",
            "AF1",
            "AF2",
            "AF3",
            "AF4",
            "AF5",
            "BLOSUM1",
            "BLOSUM2",
            "BLOSUM3",
            "BLOSUM4",
            "BLOSUM5",
            "BLOSUM6",
            "BLOSUM7",
            "BLOSUM8",
            "BLOSUM9",
            "BLOSUM10",
            "PP1",
            "PP2",
            "PP3",
            "F1",
            "F2",
            "F3",
            "F4",
            "F5",
            "F6",
            "KF1",
            "KF2",
            "KF3",
            "KF4",
            "KF5",
            "KF6",
            "KF7",
            "KF8",
            "KF9",
            "KF10",
            "MSWHIM1",
            "MSWHIM2",
            "MSWHIM3",
            "E1",
            "E2",
            "E3",
            "E4",
            "E5",
            "PD1",
            "PD2",
            "PRIN1",
            "PRIN2",
            "PRIN3",
            "ProtFP1",
            "ProtFP2",
            "ProtFP3",
            "ProtFP4",
            "ProtFP5",
            "ProtFP6",
            "ProtFP7",
            "ProtFP8",
            "SV1",
            "SV2",
            "SV3",
            "SV4",
            "ST1",
            "ST2",
            "ST3",
            "ST4",
            "ST5",
            "ST6",
            "ST7",
            "ST8",
            "SVGER1",
            "SVGER2",
            "SVGER3",
            "SVGER4",
            "SVGER5",
            "SVGER6",
            "SVGER7",
            "SVGER8",
            "SVGER9",
            "SVGER10",
            "SVGER11",
            "T1",
            "T2",
            "T3",
            "T4",
            "T5",
            "VHSE1",
            "VHSE2",
            "VHSE3",
            "VHSE4",
            "VHSE5",
            "VHSE6",
            "VHSE7",
            "VHSE8",
            "VSTPV1",
            "VSTPV2",
            "VSTPV3",
            "VSTPV4",
            "VSTPV5",
            "VSTPV6",
            "Z1",
            "Z2",
            "Z3",
            "Z4",
            "Z5",
            "cys30",
            "cys20",
            "cys_ser30",
            "cys_ser20",
            "charge",
            "avgcharge",
            "avghydrop",
            "length",
            "entropy",
            "entropyratio",
            "boman",
            "instability",
            "aliphatic",
            "A",
            "R",
            "N",
            "D",
            "C",
            "E",
            "Q",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
            "RHK",
            "DE",
            "STNQ",
            "CGP",
            "AVIL",
            "MFYW"
        ]
        header = ",".join(header_list) + "\n"
        outfile.write(header)
        with open(fasta_file) as fasta_text:
            for line in fasta_text:
                if line.startswith(">"):
                    protein_id = line[1:-1]
                    descriptors = {}
                    descriptors["protein_id"]=protein_id
                    descriptors["validation"]="yes"
                if not line.startswith(">"):
                    sequence = line[:-1]
                    peptides_features = peptides.Peptide(sequence).descriptors()
                    descriptors.update(peptides_features)
                    decript_obj = RiPP(sequence)
                    decript_obj.calculate_features()
                    decript_feat = decript_obj.__dict__
                    aa_freq = decript_feat["aafreq"]
                    clr_freq = decript_feat["clfreq"]
                    decript_feat.pop("aafreq")
                    decript_feat.pop("clfreq")
                    decript_feat.update(aa_freq)
                    decript_feat.update(clr_freq)
                    descriptors.update(peptides_features)
                    descriptors.update(decript_feat)
                    if is_ripp == "True":
                        descriptors["RiPP"] = "RiPP"
                        descriptors["Class"] = json_category_dict[protein_id]
                    else:
                        descriptors["RiPP"] = "No_RiPP"
                        descriptors["Class"] = "No_RiPP"
                    descriptor_line = []
                    for descriptor in header_list:
                        descriptor_line.append(str(descriptors[descriptor]))
                    outline = ",".join(descriptor_line) + "\n"
                    print(header)
                    print(outline)
                    outfile.write(outline)


if __name__ == "__main__":
    #usage: python3 /home/sanz006/thesis/feature_extraction_beta/all_feature_extrac.py /lustre/BIF/nobackup/sanz006/positive_training/positive_sequences.fasta /lustre/BIF/nobackup/sanz006/positive_training/complete_positive_v2.csv True
    # usage: python3 /home/sanz006/thesis/feature_extraction_beta/all_feature_extrac.py /lustre/BIF/nobackup/sanz006/negative_training/RiPP_negative_training_set.fasta /lustre/BIF/nobackup/sanz006/positive_training/complete_negative_v2.csv False
    write_dataset(argv[1], argv[2], argv[3])
