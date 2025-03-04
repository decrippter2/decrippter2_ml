

import json
import peptides
from decrippter2_ml.feature_extraction_manager.ripp_class import RiPP
from typing import Self
from pydantic import BaseModel


class SequenceManager(BaseModel):
    """Pydantic_based class to extract all features from aminoacid sequence

    Attributes:
            aa_seq: str with aminoacid sequence
            descriptors_dict: dict containing all extracted features
            header_list: list containing pandas Dataframe colnames in order

    """
    aa_seq: str
    descriptors_dict: dict = {}
    header_list : list =[]
    def read_json(self: Self, json_path:str) -> dict:
        """

        """
        with open(json_path) as json_file:
            json_dict=json.load(json_file)
        return json_dict

    def read_header(self:Self):
        header_list = self.read_json('feature_extraction_manager/header_list.json')
        self.header_list=header_list["header_list"]

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
        return self.descriptors_dict

