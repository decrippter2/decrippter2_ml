"""Entry point of decrippter2_ml

Copyright (c) 2024 Mitja M. Zdouc and individual contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import json
import logging
import sys
from importlib import metadata
import argparse
#import coloredlogs
from download_manager.download_manager import DownloadManager
from feature_extraction_manager.feature_extraction import FeatureExtractor
from classifiers_scripts.nn import NN_Classifier
from classifiers_scripts.nn_svm import NN_SVM_Classifier
from classifiers_scripts.svm import SVM_Classifier

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
    fasta_dict={}
    with open(file) as fasta:
        for line in fasta:
            if line.startswith('>'):
                prot_id=line[1:-1]
            else:
                seq=line[:-1]
            fasta_dict[prot_id]=seq
    return fasta_dict

def retrain_models():
    download=DownloadManager()
    download.download_data()
    download.organize_data()
    feature_extractor=FeatureExtractor()
    training_set=feature_extractor.build_dataset()
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
def run_prediction(seq_dict,model_type,output_name,logger):
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
    for id in seq_dict:
        output=model.predict(seq_dict[id])
        results_dict[id]=output
    with open(output_name,'w') as outfile:
        json.dump(results_dict,outfile,indent=6)



def parse_arguments():

    settings = {}
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input',
                        help='Input: either a fasta file or a single peptide sequence to analyze (specify which with -t)',
                        required=True)
    parser.add_argument('-t', '--input_type',
                        help='Type of input given with the -i flag. Either fasta or seq (default)',
                        choices=['fasta', 'seq'], default='seq')
    #parser.add_argument('-c', '--cutoff', help='The cutoff of the SVM score to use (0-1).', default=0.9, type=float)
    parser.add_argument('-o', '--outputname', help='Path in which the results will be written to (json format)', default='.')
    parser.add_argument('--redo', help='Flag to retrain the existing models with updated data')
    #parser.add_argument('--output_type', choices=['simple', 'detailed'], default='simple',
                        #help='Choose between a simple output (only headers with final score) or a detailed one (all protein features and individual scores of the three SVMs')
    #    parser.add_argument('--keep_negatives', help='Also show the precursors in the output that did not make the cutoff', default=False, action='store_true')
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
    logger.info(f"Started decrippter2_ml CLI v{metadata.version('decrippter2_ml')}.")
    logger.debug("Hello, world!")

    args,settings = parse_arguments()
    sequence_dict={}
    if settings['--redo'] == True:
        retrain_models()
    if settings['input_type'] == 'seq':
        sequence_dict['seq0000001']=settings['input']
    else:
        sequence_dict=read_fasta(settings['input'])
    run_prediction(sequence_dict,settings['model_type'],settings['outputname'],logger)




if __name__ == "__main__":
    main()
