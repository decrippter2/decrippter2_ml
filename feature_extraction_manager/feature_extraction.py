from pathlib import Path
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

    json_folder: Path = Path(__file__).parent.joinpath("json_data")
    prot_id_subclass: dict ={}
    entry_subclass: dict= {}


    def read_file(self: Self, file_path:str):
        """Reads text file

        Args:
            file_path: str, path to text file

        Returns:
            Text file content read into a variable

        """
        with open(file_path) as infile:
            output=infile
        return output

    def read_json(self: Self, json_path:str) -> dict:
        """Reads .json file content

        Args:
            json_path: str, path to json file

        Returns:
            json file content as a dictionary

        """
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
        with open(multifasta, 'w') as outfile:
            for ripp_json in self.folder_path.iterdir():
                    #we open each of the ripp jsons
                    ripp_dict=self.read_json(ripp_json)
                    for entry in ripp_dict["entries"]:
                        header='>'+str(entry["protein_ids"]["genpept"])+'\n'
                        sequence=entry["complete"]+'\n'
                        outfile.write(header)
                        outfile.write(sequence)

    def record_ripp_category(self: Self, fasta_file:str):
        """Extracts a record of ripp subclasses

        Args:
            fasta_file: a str containing the path to the multifasta file

        Returns:
            A dictionary containing the genpept accessions as keys and the
            corresponding ripp category as values
        """

        id_subclass = {}
        entry_subclass = {}

        positives_fasta=self.read_file(fasta_file)
        for line in positives_fasta:
            if line.startswith('>'):
                prot_query = line[1:-1]
                for json_file in self.json_folder.iterdir():
                        json_dict=self.read_json(json_file)
                        entry_subclass[str(json_file)]=json_dict["ripp_class"]
                        for entry in json_dict["entries"]:
                            if prot_query in entry["protein_ids"].values():
                                id_subclass[prot_query] = json_dict["ripp_class"]
        self.prot_id_subclass=id_subclass
        self.entry_subclass=entry_subclass

    def write_df_row(self:Self,row_dict: dict):
        """Writes a row to Pandas Dataframe format

        Args:
            row_dict: dict containing new row's content

        Returns:
            row_dicts content as a single-row Dataframe

        """
        index_map = {v: i for i, v in enumerate(self.header_list)}
        sorted(row_dict.items(), key=lambda pair: index_map[pair[0]])
        df_row=pd.DataFrame(row_dict.values(),columns=self.header_list)
        return df_row

    def write_dataset(self: Self, fasta_file:str, is_ripp:bool,is_validated:bool):
        """Writes the dataset extracting the features from sequence

        Args:
            fasta_file: str containing path to multifasta file
            is_ripp: True is multifasta corresponds to RiPP sequences, False if not

        Returns:
            A string with the path of the saved .csv file

        """
        header_list=self.read_json('feature_extraction_manager/header_list.json')
        header_list=header_list["header_list"]
        fasta_text=self.read_file(fasta_file)
        ripp_dataframe=pd.DataFrame(columns=header_list)
        for line in fasta_text:
            if line.startswith(">"):
                protein_id = line[1:-1]
                descriptors = {}
                descriptors["protein_id"] = protein_id
                if is_ripp and is_validated:
                    descriptors["validation"] = "yes"
                    descriptors["RiPP"] = "RiPP"
                    descriptors["Class"] = self.prot_id_subclass[protein_id]
                elif is_ripp and not is_validated:
                    descriptors["validation"] = "no"
                    descriptors["RiPP"] = "RiPP"
                    descriptors["Class"] = self.entry_subclass[fasta_file[:5]]#needs testing
                elif not is_ripp and is_validated:
                    descriptors["validation"] = "yes"
                    descriptors["RiPP"] = "No_RiPP"
                    descriptors["Class"] = "No_RiPP"
            if not line.startswith(">"):
                sequence = line[:-1]

                feature_extractor=SequenceManager(aa_seq=sequence)
                extracted_features=feature_extractor.calculate_features()
                descriptors.update(extracted_features)

                new_row=self.write_df_row(descriptors)
                ripp_dataframe= pd.concat([ripp_dataframe,new_row], ignore_index=True)

        return ripp_dataframe

    def build_dataset(self,fasta_file,is_ripp,is_validated):
        self.write_multifasta(fasta_file)
        self.record_ripp_category(fasta_file)
        dataframe_obj = self.write_dataset(fasta_file, is_ripp, is_validated)


        return dataframe_obj