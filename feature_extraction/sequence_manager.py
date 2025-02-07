
import os
import json
import peptides
from ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd

class SequenceManager(BaseModel):
    """

    """
    aa_seq: str | None = None
    descriptors_dict: {}
    header_list = list
    def read_json(self: Self, json_path:str) -> dict:
        with open(json_path) as json_file:
            json_dict=json.load(json_file)
        return json_dict

    def read_header(self:Self):
        header_list = self.read_json('feature_extraction/header_list.json')
        self.header_list=header_list["header_list"]
    def secondary_structure(self:Self):
        alpha_leader,alpha_50_total,alpha_total=self.aa_seq,self.aa_seq,self.aa_seq
        beta_leader, beta_50_total, beta_total = self.aa_seq, self.aa_seq, self.aa_seq
        coil_leader, coil_50_total, coil_total = self.aa_seq, self.aa_seq, self.aa_seq
        structure_cols=[alpha_leader,alpha_50_total,alpha_total,beta_leader, beta_50_total,
                        beta_total,coil_leader, coil_50_total, coil_total]
        secondary_structure_dict = {
            "Alpha50_leader": alpha_leader,
            "Alpha50_total": alpha_50_total,
            "Total_ssh dodnalpha":alpha_total,
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

