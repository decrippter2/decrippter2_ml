
import os
import json
import peptides
from ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd
from sequence_manager import SequenceManager

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

    def write_df_row(self:Self,row_dict: dict):
        index_map = {v: i for i, v in enumerate(self.header_list)}
        sorted(row_dict.items(), key=lambda pair: index_map[pair[0]])
        df_row=pd.DataFrame(row_dict.values(),columns=self.header_list)
        return df_row

    def write_dataset(self: Self, fasta_file:str, is_ripp:bool):
        """Writes the dataset extracting the features from sequence

        Args:
            fasta_file: str containing path to multifasta file
            is_ripp: True is multifasta corresponds to RiPP sequences, False if not

        Returns:
            A string with the path of the saved .csv file

        """
        header_list=self.read_json('feature_extraction/header_list.json')
        header_list=header_list["header_list"]
        fasta_text=self.read_file(fasta_file)
        ripp_dataframe=pd.DataFrame(columns=header_list)
        for line in fasta_text:
            if line.startswith(">"):
                protein_id = line[1:-1]
                descriptors = {}
                descriptors["protein_id"] = protein_id
                descriptors["validation"] = "yes"
            if not line.startswith(">"):
                sequence = line[:-1]

                feature_extractor=SequenceManager(aa_seq=sequence)
                extracted_features=feature_extractor.calculate_features()
                descriptors.update(extracted_features)

                if is_ripp:
                    descriptors["RiPP"] = "RiPP"
                    descriptors["Class"] = self.ripp_categories[protein_id]
                else:
                    descriptors["RiPP"] = "No_RiPP"
                    descriptors["Class"] = "No_RiPP"
                new_row=self.write_df_row(descriptors)
                ripp_dataframe= pd.concat([ripp_dataframe,new_row], ignore_index=True)

        return ripp_dataframe

