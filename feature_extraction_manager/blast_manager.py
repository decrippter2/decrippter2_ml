
import os
import json
import peptides
from ripp_class import RiPP
import logging
from typing import Any, Self
from pydantic import BaseModel
from feature_extraction.feature_extraction import FeatureExtractor
import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import time
from pathlib import Path
class BlastManager(FeatureExtractor):
    """
    
    """

    ncbi_results: Path = Path(__file__).parent.joinpath("blast_nr_results")
    fasta: Path = Path(__file__).parent.joinpath("data/fasta")

    def write_fasta_file(self,filename:str,sequence:str):
        with open(self.fasta.joinpath(f"{filename}.fa"), "w") as file:
            file.write(f">{filename}\n")
            file.write(f"{sequence}")
    def jsons_to_fastas(self):
        for __, __, filenames in os.walk(
                self.folder_path
        ):
            for filename in filenames:
                file = (
                        self.folder_path
                        + filename
                )
                file_dict=self.read_json(file)
                n_entries=len(file_dict["entries"])
                for i in range(n_entries):
                    subfilename=filename[:-5]+'_'+str(i+1)
                    sequence=file_dict["entries"][i]["complete"]
                    self.write_fasta_file(subfilename,sequence)
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

    def extract_xml(self, acc: str):
        """Extracts metadata for SSN from BLAST XML file

        Counts matches >= 95% similarity and adds to efi-est metadata
        Also stores matches >= 95% similarity for dumping as csv

        Arguments:
            acc: a MITE accession
        """
        with open(self.ncbi_results.joinpath(f"{acc}.xml")) as xml_file:
            blast_record = NCBIXML.read(xml_file)

        counter = 0
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sim_perc = round((hsp.positives / hsp.align_length) * 100, 2)
                id_perc = round((hsp.identities / hsp.align_length) * 100, 2)

                if id_perc >= 50:
                    self.nr_blast_matches["mite_acc"].append(acc)
                    self.nr_blast_matches["accession"].append(alignment.accession)
                    self.nr_blast_matches["length"].append(alignment.length)
                    self.nr_blast_matches["e_value"].append(hsp.expect)
                    self.nr_blast_matches["score"].append(hsp.score)
                    self.nr_blast_matches["bitscore"].append(hsp.bits)
                    self.nr_blast_matches["percent_sim"].append(sim_perc)
                    self.nr_blast_matches["percent_id"].append(id_perc)

                    counter += 1

        self.metadata_efi_est["ncbi_nr_matches"].append(counter)

    def extract_xml2(self, acc: str):
        """Extracts accessions and sequences from XML file
        """
        with open(self.ncbi_results.joinpath(f"{acc}.xml")) as xml_file:
            blast_record = NCBIXML.read(xml_file)
        results={}
        counter = 0
        print(blast_record)
        """
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sim_perc = round((hsp.positives / hsp.align_length) * 100, 2)
                id_perc = round((hsp.identities / hsp.align_length) * 100, 2)

                if id_perc >= 50:
                    self.nr_blast_matches["mite_acc"].append(acc)
                    self.nr_blast_matches["accession"].append(alignment.accession)
                    self.nr_blast_matches["length"].append(alignment.length)
                    self.nr_blast_matches["e_value"].append(hsp.expect)
                    self.nr_blast_matches["score"].append(hsp.score)
                    self.nr_blast_matches["bitscore"].append(hsp.bits)
                    self.nr_blast_matches["percent_sim"].append(sim_perc)
                    self.nr_blast_matches["percent_id"].append(id_perc)

                    results[alignment.accession]=blast_record.
"""