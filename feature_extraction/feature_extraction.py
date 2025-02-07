
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
    ripp_categories: dict | None= None


    def read_file(self: Self, file_path:str):
        with open(file_path) as infile:
            output=infile
        return output

    def read_json(self: Self, json_path:str) -> dict:
        with open(json_path) as json_file:
            json_dict=json.load(json_file)
        return json_dict

    def write_multifasta(self: Self, multifasta: str) -> str:
        """Write multifasta file with all ripp sequences

        Args:
            multifasta: name of multifasta file

        Returns:
            The name of multifasta file

        """
        for __, __, filenames in os.walk(
                self.folder_path
        ):
            with open(multifasta,'w') as outfile:
                for filename in filenames:
                    file = (
                            self.folder_path
                            + filename
                    )
                    #we open each of the ripp jsons
                    ripp_dict=self.read_json(file)
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
        positives_fasta=self.read_file(fasta_file)
        for line in positives_fasta:
            if line.startswith('>'):
                prot_query = line[1:-1]
                for __, __, filenames in os.walk(self.folder_path):
                    for filename in filenames:
                        json_path=self.folder_path+filename
                        json_dict=self.read_json(json_path)
                        for entry in json_dict["entries"]:
                            if prot_query in entry["protein_ids"].values():
                                entry_cat_dict[prot_query] = json_dict["ripp_class"]
        self.ripp_categories=entry_cat_dict

    def write_dataset(self: Self, fasta_file:str, outfile_name:str, is_ripp:bool) -> str:
        """Writes the dataset extracting the features from sequence

        Args:
            fasta_file: str containing path to multifasta file
            outfile_name: str stating path of generated csv file
            is_ripp: True is multifasta corresponds to RiPP sequences, False if not

        Returns:
            A string with the path of the saved .csv file

        """
        header_list=self.read_json('feature_extraction/header_list.json')
        header_list=header_list["header_list"]
        with open(outfile_name, "w") as outfile:
            header = ",".join(header_list) + "\n"
            outfile.write(header)
            fasta_text=self.read_file(fasta_file)
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
                        descriptors["Class"] = self.ripp_categories[protein_id]
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

class SequenceManager(FeatureExtractor):

    aa_seq: str | None = None
    def
