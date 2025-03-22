###STILL NEEDS TO BE TESTED
from decrippter2_ml_src import FeatureExtractor, BlastManager

from typing import Self
import pandas as pd
class DataAugmentationManager(FeatureExtractor):
    """Pydantic_based class to perform data augmentation on validated data
    """

    def expansion_multifasta(self:Self):
        """Calls BlastManager to perform dataset expansion
        """
        blast_obj=BlastManager()
        blast_obj.jsons_to_fastas()
        blast_obj.run()
        for xml in blast_obj.ncbi_results.iterdir():
            blast_obj.extract_xml(xml[:-2]) #all xmls file results are converted to multifasta files

    def erase_duplicates(self:Self,dataframe):
        """Function to erase duplicates from dataset after augmentation,
        prioritizing 'validated' tags

        Args:
            dataframe: pd dataframe containing full dataset with extracted features

        Returns:
            pd dataframe with erased duplicates
            """
        sequences = []
        drop_indices = []
        for index, row in dataframe[dataframe.validation == "yes"].iterrows():
            print(index)
            sequences.append(row["sequence"])
        for index, row in dataframe[dataframe.validation == "no"].iterrows():
            print(index)
            if row["sequence"] not in sequences:
                sequences.append(row["sequence"])
            else:
                drop_indices.append(index)
        dataframe = dataframe.drop(drop_indices)
        return dataframe

    def expand_dataset(self:Self):
        """Function to carry on dataset augmentation

        Returns:
            Expanded DataFrame as Pandas object"""
        self.expansion_multifasta()
        #build validated dataset
        positive_dataset=self.build_dataset('zenodo_ripp.fa',True,True)
        negative_dataset=self.build_dataset('negative.fa',False,True)
        blast_obj=BlastManager()
        augmented_data=pd.DataFrame()
        for fasta_file in blast_obj.ncbi_results_fasta.iterdir():
            result_dataframe=self.build_dataset(fasta_file,True,False)
            augmented_data=pd.concat([augmented_data,result_dataframe])
        expanded_dataset=pd.concat([positive_dataset,negative_dataset,augmented_data])#merge all into big dataframe
        return expanded_dataset #returns df as pd object



