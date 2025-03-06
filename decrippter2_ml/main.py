import json
import logging
import sys
from importlib import metadata
import argparse

import pandas as pd

from download_manager.downloader import DownloadManager
from feature_extraction_manager.feature_extraction import FeatureExtractor
from classifiers_scripts.nn import NN_Classifier
from classifiers_scripts.nn_svm import NN_SVM_Classifier
from classifiers_scripts.svm import SVM_Classifier
from data_augmentation.augment_dataset import DataAugmentationManager

def config_logger() -> logging.Logger:
    """Set up a named logger with nice formatting

    :return
        A Logger object
    """
    logger = logging.getLogger("decrippter2_ml")
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler(sys.stdout)
    """
    console_handler.setFormatter(
        coloredlogs.ColoredFormatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
    )"""
    logger.addHandler(console_handler)
    return logger

def read_fasta(file):
    with open(file) as fasta:
            fasta_lines=fasta.readlines()
    return fasta_lines

def retrieve_data():
    #download=DownloadManager()
    #download.download_data()
    #download.organize_data()
    data_augmenter=DataAugmentationManager()
    training_set=data_augmenter.expand_dataset()
    return training_set
def retrain_models():
    training_set=retrieve_data()
    nn_obj=NN_Classifier(dataset=training_set)
    nn_obj.retrain_model()
    nn_svm_poly3=NN_SVM_Classifier(model_type='poly3',dataset=training_set)
    nn_svm_poly3.retrain_model()
    nn_svm_rbf=NN_SVM_Classifier(model_type='rbf',dataset=training_set)
    nn_svm_rbf.retrain_model()
    svm_poly3=SVM_Classifier(model_type='poly3',dataset=training_set)
    svm_poly3.retrain_model()
    svm_rbf = SVM_Classifier(model_type='rbf', dataset=training_set)
    svm_rbf.retrain_model()
def run_prediction(fasta_lines,model_type,output_name,logger):
    if model_type=='nn':
        model=NN_Classifier()
    elif model_type=='nn_svm_poly3':
        model=NN_SVM_Classifier(model_type='poly3')
    elif model_type=='nn_svm_rbf':
        model=NN_SVM_Classifier(model_type='rbf')
    elif model_type=='svm_poly3':
        model=NN_SVM_Classifier(model_type='poly3')
    elif model_type=='svm_rbf':
        model=NN_SVM_Classifier(model_type='rbf')
    else:
        logger.fatal('Error: invalid kernel')
        raise RuntimeError
    results_dict={}
    feature_extractor=FeatureExtractor()
    feature_pd=feature_extractor.write_dataset(fasta_lines,False,False)
    output=model.predict(feature_pd)
    output=pd.DataFrame({'output':output })
    feature_pd=pd.concat([feature_pd,output],axis=1)
    for line in feature_pd.iterrows():
        results_dict[line[1].protein_id]=line[1].output
    with open(output_name,'w') as outfile:
        json.dump(results_dict,outfile,indent=6)



def parse_arguments():
    """Function to parse all arguments"""
    settings = {}
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',
                        help='Input: either a fasta file or a single peptide sequence to analyze (specify which with -t)')
    parser.add_argument('-t', '--input_type',
                        help='Type of input given with the -i flag. Either fasta or seq (default)',
                        choices=['fasta', 'seq'], default='seq')
    parser.add_argument('-o', '--outputname', help='Path in which the results will be written to (add .json extension)', default='.')
    parser.add_argument('--redo', help='Flag to retrain the existing models with updated data',action=argparse.BooleanOptionalAction)
    parser.add_argument('--models_path',
                        help='Point to the folder containing the pretrained models (downloaded with decRiPPter)',
                        default='pretrained_models/')
    parser.add_argument('-m', '--model_type',
                        help='Type of model used from RiPP prediction',
                        choices=['nn','nn_svm_poly3','nn_svm_rbf','svm_poly3','svm_rbf'])

    args = parser.parse_args()
    for key, value in args.__dict__.items():
        # Update only if the args value is not None
        if value != None:
            settings[key] = value

    return args, settings

def main() -> None:
    """Function to execute main body of code"""
    logger = config_logger()
    logger.debug("Starting decRiPPter2...")

    args,settings = parse_arguments()
    if settings['redo'] == True:
        retrain_models()
    else:
        if settings['input_type'] == 'seq':
            sequence_lines=['>seq0000001\n',settings['input']]
        else:
            sequence_lines=read_fasta(settings['input'])
        run_prediction(sequence_lines,settings['model_type'],settings['outputname'],logger)




if __name__ == "__main__":
    main()
