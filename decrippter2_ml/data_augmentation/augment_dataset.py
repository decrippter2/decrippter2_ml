import pandas as pd
from decrippter2_ml.feature_extraction_manager.feature_extraction import FeatureExtractor
from decrippter2_ml.feature_extraction_manager.blast_manager import BlastManager


class DataAugmentationManager(FeatureExtractor):


    def expansion_multifasta(self):
        blast_obj=BlastManager()
        blast_obj.jsons_to_fastas()
        blast_obj.run()
        for xml in blast_obj.ncbi_results.iterdir():
            blast_obj.extract_xml(xml[:-2]) #all xmls file results are converted to multifasta files

    def erase_duplicates(self,dataframe):
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
    def expand_dataset(self):
        self.expansion_multifasta()
        #build validated dataset
        positive_dataset=self.build_dataset('zenodo_ripp.fa',True,True)
        negative_dataset=self.build_dataset('negative.fa',False,True)
        #record subclasses for blast results

        #build augmented fraction of dataset
        blast_obj=BlastManager()
        augmented_data=pd.DataFrame()
        for fasta_file in blast_obj.ncbi_results_fasta.iterdir():
            result_dataframe=self.build_dataset(fasta_file,True,False)
            augmented_data=pd.concat([augmented_data,result_dataframe])
        expanded_dataset=pd.concat([positive_dataset,negative_dataset,augmented_data])#merge all into big dataframe
        return expanded_dataset #returns df as pd object



