import os
import json
import peptides
from ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd
from subprocess import run
from pathlib import Path
"""
DRAFT CODE:

nsp_path='/lustre/BIF/nobackup/sanz006/nsp3/nsp3.py'
model_path='/lustre/BIF/nobackup/sanz006/nsp3/models/nsp3.pth'
input='/lustre/BIF/nobackup/sanz006/positive_training/short_multifasta.txt'
output_folder='/lustre/BIF/nobackup/sanz006/positive_training/nsp3_test'
worker_id='01'
run(f'conda activate nsp3')
run(f'python3 {nsp_path} -m {model_path} -i {input} -o {output_folder} -w {worker_id}')
with open(f'{output_folder}/{worker_id}/{worker_id}.json') as results_file:
    results_dict=json.load(results_file)
structure_dict={}
for result in results_dict:
    structure_dict[result["desc"]]=result["q8"]

"""


class NSP3Manager(BaseModel):
    """

    """
    nsp_path: Path=Path('/lustre/BIF/nobackup/sanz006/nsp3/nsp3.py') #I DONT KNOW WHERE THIS WILL BE STORED IN THE FINAL PROGRAM
    model_path: Path=Path('/lustre/BIF/nobackup/sanz006/nsp3/models/nsp3.pth')
    input_fasta: Path
    output_folder:Path=Path('/lustre/BIF/nobackup/sanz006/positive_training/nsp3_test')
    worker_id:str='01'

    def run(self):
        run(f'conda activate nsp3') #NEEDS INSTALLATION BEFOREHAND OF NSP3 IN ENVIRONMENT
        run(f'python3 {self.nsp_path} -m {self.model_path} -i {self.input_fasta} -o {self.output_folder} -w {self.worker_id}')
        with open(f'{self.output_folder}/{self.worker_id}/{self.worker_id}.json') as results_file:
            results_dict=json.load(results_file)
        sec_str=results_dict[0]["q8"]
        return sec_str
