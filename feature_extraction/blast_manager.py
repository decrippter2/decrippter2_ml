
import os
import json
import peptides
from ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import time
class BlastManager(BaseModel):
    """
    
    """
    def run(self):
        self.ncbi_results.mkdir(exist_ok=True)

        for fasta in self.fasta.iterdir():
            if self.ncbi_results.joinpath(f"{fasta.stem}.xml").exists():
                continue

            record = SeqIO.read(fasta, "fasta")
            query = f""">{record.id}
            {record.seq}
            """

            result_handle = NCBIWWW.qblast(
                program="blastp",
                database="nr",
                sequence=query,
                expect=1e-5,
                hitlist_size=5000,  # Maximum results
                format_type="XML",
            )
            raw_blast_output = result_handle.read()
            with open(self.ncbi_results.joinpath(f"{record.id}.xml"), "w") as file:
                file.write(raw_blast_output)

            time.sleep(10)  # limits rate to prevent IP ban by NCBI