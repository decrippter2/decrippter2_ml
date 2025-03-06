### STILL NEEDS TO BE TESTED
import os
from decrippter2_ml.feature_extraction_manager.feature_extraction import FeatureExtractor
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import time
from pathlib import Path
from typing import Self
class BlastManager(FeatureExtractor):
    """Pydantic_based class to perform a BLAST search over validated data

    Attributes:
        ncbi_results: path to directory where BLAST results will be stored
        fasta: path to directory containing all RiPP sequence entries in separated FASTA files
        ncbi_results_fasta: path to
    
    """

    ncbi_results: Path = Path(__file__).parent.joinpath("blast_nr_results")
    fasta: Path = Path(__file__).parent.joinpath("data/fasta")
    ncbi_results_fasta: Path = Path(__file__).parent.joinpath("blast_nr_results_fasta")

    def write_fasta_file(self,filename:str,sequence:str):
        """Function for fasta writing

        Args:
            filename: str, name of fasta file
            sequence: str, AA sequence of RiPP entry
        """
        with open(self.fasta.joinpath(f"{filename}.fa"), "w") as file:
            file.write(f">{filename}\n")
            file.write(f"{sequence}")
    def jsons_to_fastas(self):
        """Turns all protein entries into multiple fasta files"""
        for __, __, filenames in os.walk(
                self.folder_path
        ):
            for filename in filenames: #loops through each one of the json files
                file = (
                        self.folder_path
                        + filename #only json without path
                )
                file_dict=self.read_json(file)
                n_entries=len(file_dict["entries"])
                for i in range(n_entries):
                    subfilename=filename[:-5]+'_'+str(i+1) #fasta file name eg: ripp0000001_1.fa
                    sequence=file_dict["entries"][i]["complete"]
                    self.write_fasta_file(subfilename,sequence)
    def run(self:Self):
        """Runs blast search"""
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
                hitlist_size=200,  # Maximum results
                format_type="XML",
            )
            raw_blast_output = result_handle.read()
            with open(self.ncbi_results.joinpath(f"{record.id}.xml"), "w") as file:
                file.write(raw_blast_output)

            time.sleep(10)  # limits rate to prevent IP ban by NCBI

    def extract_xml(self, acc: str):
        """Extracts accessions and sequences from XML file
        """
        with open(self.ncbi_results.joinpath(f"{acc}.xml")) as xml_file:
            blast_record = NCBIXML.read(xml_file)
        results={}
        with open(self.ncbi_results_fasta.joinpath(f"{acc}.fa"),'w') as out_fasta:
            for alignment in blast_record.alignments:
                print(acc)
                print(alignment.accession)
                accession=alignment.accession
                for hsp in alignment.hsps:

                    id_perc = round((hsp.identities / hsp.align_length) * 100, 2)
                    if id_perc >=50:
                        out_fasta.write(f">{accession}\n")
                        out_fasta.write(f"{sequence}\n")
                        print(hsp.sbjct)
                        sequence=hsp.sbjct
                        results[accession]=sequence
