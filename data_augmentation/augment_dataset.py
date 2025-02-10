import os
import json
import peptides
from feature_extraction.ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd
from feature_extraction.sequence_manager import SequenceManager
from feature_extraction.feature_extraction import FeatureExtractor
from feature_extraction.blast_manager import BlastManager

class DataAugmentationManager(FeatureExtractor):

    blast_obj=BlastManager()
    blast_obj.jsons_to_fastas()
    blast_obj.run()


    def blast_subclasses(self):

    def expand_dataset(self):
        #build validated dataset
        positive_dataset=self.build_validated_dataset('zenodo_ripp.fa',True,True)
        negative_dataset=self.build_validated_dataset('negative.fa',False,True)
        #record subclasses for blast results

        #build augmented fraction of dataset

        #merge all into big dataframe


