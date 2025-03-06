from pathlib import Path
import json
from typing import Self
from pydantic import BaseModel
import pandas as pd
from .sequence_manager import SequenceManager

class FeatureExtractor(BaseModel):
    """Pydantic_based class to extract all features from folder

    Attributes:
            folder_path=path to JSON with RiPP entries
    """

    json_folder: Path = Path(__file__).parent.parent.joinpath("json_data")
    header_path: Path = Path(__file__).parent.joinpath("header_list.json")
    header_list: list = []
    feature_list: Path = Path(__file__).parent.joinpath("feature_list.json")
    prot_id_subclass: dict ={}
    entry_subclass: dict= {}


    def read_file(self: Self, file_path:str):
        """Reads text file

        Args:
            file_path: str, path to text file

        Returns:
            Text file content read into a variable

        """
        if type(file_path)==list:
            output=file_path
        else:
            with open(file_path) as infile:
                output=infile.readlines()
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

    def record_ripp(self:Self, multifasta: str):
        """Records RiPP subclass of each RiPP entry and protein id

        Args:
            multifasta: str, path to multifasta file to be written
        """
        with open(multifasta,'w') as outfile:
            for ripp_json in self.json_folder.iterdir():
                ripp_dict=self.read_json(ripp_json)
                subclass=ripp_dict["ripp_class"]
                self.entry_subclass[ripp_json] = subclass
                for entry in ripp_dict["entries"]:
                    protein_id=entry["protein_ids"]["genpept"]
                    header='>'+protein_id+'\n'
                    sequence=entry["complete"]+'\n'
                    outfile.write(header)
                    outfile.write(sequence)
                    self.prot_id_subclass[protein_id]=subclass


    def write_df_row(self:Self,row_dict: dict):
        """Writes a row to Pandas Dataframe format

        Args:
            row_dict: dict containing new row's content

        Returns:
            row_dicts content as a single-row Dataframe

        """
        index_map = {v: i for i, v in enumerate(self.header_list)}
        sorted(row_dict.items(), key=lambda pair: index_map[pair[0]])
        df_row=pd.DataFrame([row_dict],columns=self.header_list)
        return df_row

    def write_dataset(self: Self, fasta_file:str, is_ripp:bool,is_validated:bool):
        """Writes the dataset extracting the features from sequence

        Args:
            fasta_file: str containing path to multifasta file
            is_ripp: True if multifasta corresponds to RiPP sequences, False if not

        Returns:
            A string with the path of the saved .csv file

        """
        header_list=self.read_json(self.header_path)
        self.header_list=header_list["header_list"]
        fasta_text=self.read_file(fasta_file)
        ripp_dataframe=pd.DataFrame(columns=self.header_list)
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
        """Builds a Pandas Dataframe object when new data comes out

        Args:
            fasta_file: str containing path to multifasta file
            is_ripp: bool, True is multifasta corresponds to RiPP sequences, False if not
            is_validated: bool, True is multifasta corresponds to validated RiPP sequences, False if not

        Returns:
            A Pandas dataframe object
        """
        #self.write_multifasta(fasta_file)
        self.record_ripp(fasta_file)
        dataframe_obj = self.write_dataset(fasta_file, is_ripp, is_validated)


        return dataframe_obj