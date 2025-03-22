import json
import peptides
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd
from blast_manager import BlastManager

blast_obj=BlastManager()
blast_obj.extract_xml2('MITE0000001')
