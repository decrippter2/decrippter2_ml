import json
import peptides
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd
from pathlib import Path
import os
"""
from decrippter2_ml import BlastManager

blast_obj=BlastManager()
blast_obj.extract_xml('MITE0000001')
"""
for file in Path(__file__).parent.joinpath("blast_nr_results").iterdir():
    print(file)
    if str(file)[len(str(file))-3:] == ".fa":
        os.remove(file)