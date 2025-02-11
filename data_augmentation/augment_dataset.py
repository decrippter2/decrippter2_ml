import os
import json
import peptides
from feature_extraction_manager.ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd
from feature_extraction_manager.sequence_manager import SequenceManager
from feature_extraction_manager.feature_extraction import FeatureExtractor
from feature_extraction_manager.blast_manager import BlastManager
from pathlib import Path

class DataAugmentationManager(FeatureExtractor):


    def expansion_multifasta(self):
        blast_obj=BlastManager()
        blast_obj.jsons_to_fastas()
        blast_obj.run()
        for xml in blast_obj.ncbi_results.iterdir():
            ripp_id=xml[:11]



    def blast_subclasses(self):

    def expand_dataset(self):
        #build validated dataset
        positive_dataset=self.build_validated_dataset('zenodo_ripp.fa',True,True)
        negative_dataset=self.build_validated_dataset('negative.fa',False,True)
        #record subclasses for blast results

        #build augmented fraction of dataset

        #merge all into big dataframe


