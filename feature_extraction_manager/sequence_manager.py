
import os
from pathlib import Path
import json
import peptides
from ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
from nsp3_manager import NSP3Manager
import pandas as pd

class SequenceManager(BaseModel):
    """Pydantic_based class to extract all features from aminoacid sequence

    Attributes:
            aa_seq: str with aminoacid sequence
            descriptors_dict: dict containing all extracted features
            header_list: list containing pandas Dataframe colnames in order

    """
    aa_seq: str
    descriptors_dict: dict = {}
    header_list : list
    sec_str: str
    def read_json(self: Self, json_path:str) -> dict:
        """

        """
        with open(json_path) as json_file:
            json_dict=json.load(json_file)
        return json_dict

    def read_header(self:Self):
        header_list = self.read_json('feature_extraction_manager/header_list.json')
        self.header_list=header_list["header_list"]

    def calculate_alpha(self:Self):
        ahc_leader = (self.sec_str[:50].count('H')) / len(self.sec_str[:50])
        ahc50_total = (self.sec_str[:50].count('H')) / len(self.sec_str)
        ahc_total = (self.sec_str.count('H')) / len(self.sec_str)
        return ahc_leader, ahc50_total, ahc_total

    def calculate_beta(self:Self):
        bsc_leader = (self.sec_str[:50].count('E')) / len(self.sec_str[:50])
        bsc50_total = (self.sec_str[:50].count('E')) / len(self.sec_str)
        bsc_total = (self.sec_str.count('E')) / len(self.sec_str)
        return bsc_leader, bsc50_total, bsc_total

    def calculate_coil(self:Self):
        cc_leader = (self.sec_str[:50].count('C')) / len(self.sec_str[:50])
        cc50_total = (self.sec_str[:50].count('C')) / len(self.sec_str)
        cc_total = (self.sec_str.count('C')) / len(self.sec_str)
        return cc_leader, cc50_total, cc_total

    def secondary_structure(self:Self):
        """
        Generates secondary structure information, to be improved while I
        get NSP3 to work
        """

        path_to_fasta=Path('temp/sect_str.txt')
        output_folder=Path('temp')
        output_folder.mkdir()
        with open(path_to_fasta) as temp_fasta:
            temp_fasta.write('>protein001\n')
            temp_fasta.write(f'{self.aa_seq}')
        if len(self.aa_seq)>=10:
            nsp3_predictor=NSP3Manager(input_fasta=path_to_fasta, output_folder=output_folder)
            self.sec_str=nsp3_predictor.run()
        else:
            self.sec_str='N'
        alpha_leader,alpha_50_total,alpha_total=self.calculate_alpha()
        beta_leader, beta_50_total, beta_total = self.calculate_beta()
        coil_leader, coil_50_total, coil_total = self.calculate_coil()
        secondary_structure_dict = {
            "Alpha50_leader": alpha_leader,
            "Alpha50_total": alpha_50_total,
            "Total_alpha":alpha_total,
            "Beta50_leader":beta_leader,
            "Beta50_total":beta_50_total,
            "Total_beta":beta_total,
            "Coil50_leader":coil_leader,
            "Coil50_total":coil_50_total,
            "Total_coil":coil_total
        }
        return secondary_structure_dict
    def calculate_features(self: Self):
        peptides_features=peptides.Peptide(self.aa_seq).descriptors()
        self.descriptors_dict.update(peptides_features)
        decrippt=RiPP(self.aa_seq)
        decrippt.calculate_features()
        decript_feat = decrippt.__dict__
        aa_freq = decript_feat["aafreq"]
        clr_freq = decript_feat["clfreq"]
        decript_feat.pop("aafreq")
        decript_feat.pop("clfreq")
        decript_feat.update(aa_freq)
        decript_feat.update(clr_freq)
        self.descriptors_dict.update(peptides_features)
        self.descriptors_dict.update(decript_feat)
        secondary_str=self.secondary_structure()
        self.descriptors_dict.update(secondary_str)
        return self.descriptors_dict

