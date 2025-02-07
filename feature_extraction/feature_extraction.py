
import os
import json
import peptides
from ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd

class FeatureExtractor(BaseModel):
    """Pydantic_based class to extract all features from folder

    Attributes:
            folder_path=path to JSON with RiPP entries
    """

    folder_path: str | None = None

    def write_multifasta(self: Self, multifasta: str) -> str:
        """Write multifasta file with all ripp sequences

        Args:
            multifasta: name of multifasta file

        Returns:
            The name of multifasta file

        """
        for _dirpath, _dirnames, filenames in os.walk(
                self.folder_path
        ):
            with open(multifasta,'w') as outfile:
                for filename in filenames:
                    file = (
                            self.folder_path
                            + filename
                    )
                    #we open each of the ripp jsons
                    with open(file) as ripp_entry:
                        ripp_dict = json.load(ripp_entry)
                        for entry in ripp_dict["entries"]:
                            header='>'+str(entry["protein_ids"]["genpept"])+'\n'
                            sequence=entry["complete"]+'\n'
                            outfile.write(header)
                            outfile.write(sequence)
        return multifasta
    def record_ripp_category(self: Self, fasta_file:str) -> dict:
        """Extracts a record of ripp subclasses

        Args:
            fasta_file: a str containing the path to the multifasta file

        Returns:
            A dictionary containing the genpept accessions as keys and the
            corresponding ripp category as values
        """

        entry_cat_dict = {}
        with open(fasta_file) as positives_fasta:
            for line in positives_fasta:
                if line.startswith('>'):
                    prot_query = line[1:-1]
                    for _dirpath, _dirnames, filenames in os.walk(self.folder_path):
                        for filename in filenames:
                            with open(self.folder_path + filename) as json_file:
                                json_dict = json.load(json_file)
                                for entry in json_dict["entries"]:
                                    if prot_query in entry["protein_ids"].values():
                                        entry_cat_dict[prot_query] = json_dict["ripp_class"]
        return entry_cat_dict

    def write_dataset(self: Self, fasta_file:str, outfile_name:str,ripp_classes:dict, is_ripp:bool) -> str:
        """Writes the dataset extracting the features from sequence

        Args:
            fasta_file: str containing path to multifasta file
            outfile_name: str stating path of generated csv file
            ripp_classes: dictionary containing genpept accessions and
            corresponding ripp families, generated with record_ripp_category
            is_ripp: True is multifasta corresponds to RiPP sequences, False if not

        Returns:
            A string with the path of the saved .csv file

        """
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
                        descriptors["protein_id"] = protein_id
                        descriptors["validation"] = "yes"
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
                            descriptors["Class"] = ripp_classes[protein_id]
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
        return outfile_name