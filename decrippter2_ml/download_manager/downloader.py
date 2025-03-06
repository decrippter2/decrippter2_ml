###STILL TO BE TESTED WITH ACTUAL ZENODO RECORD
from pydantic import BaseModel
from pathlib import Path
import logging
import requests
import json
import os
import shutil

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)
class DownloadManager(BaseModel):
    """Pydantic_based class to download files from Zenodo record and unpack

    Attributes:
        record: the record to download
        location: the location to download data to
        record: path to the record file
        record_unzip: path to unzipped record file
        version: path to file containing the version of mite_data used
    """

    record: str
    location: Path = Path(__file__).parent.joinpath("data")
    record: Path = Path(__file__).parent.joinpath("data/record.zip")
    record_unzip: Path = Path(__file__).parent.joinpath("data/record")
    version: Path = Path(__file__).parent.joinpath("version.json")

    def run(self) -> None:
        """Call methods for downloading and moving data"""

        if self.location.exists():
            logger.warning(
                "RiPP data folder already present - skip download. Remove the folder manually if a new version needs to be downloaded."
            )
            return

        self.location.mkdir(parents=True)
        self.download_data()
        self.organize_data()

    def download_data(self) -> None:
        """Download data from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        response_metadata = requests.get(
            f"https://zenodo.org/api/records/{self.record}"
        )
        if response_metadata.status_code != 200:
            logger.fatal(
                f"Error fetching 'ripp_data' record metadata: {response_metadata.status_code}"
            )
            raise RuntimeError

        record_metadata = response_metadata.json()
        version = record_metadata["metadata"]["version"]
        files_url = record_metadata["files"][0]["links"]["self"]

        response_data = requests.get(files_url)
        if response_data.status_code != 200:
            logger.fatal(
                f"Error downloading 'ripp_data' record: {response_data.status_code}"
            )
            raise RuntimeError

        with open(self.version, "w") as f:
            f.write(json.dumps({"version_decriptter2_data_used": f"{version}"}))

        with open(self.record, "wb") as f:
            f.write(response_data.content)

    def organize_data(self) -> None:
        """Unpacks data, moves to convenient location, cleans up

        Raises:
            NotADirectoryError: directory not unzipped in expected location
            RuntimeError: Could not determine data location in downloaded folder
        """
        shutil.unpack_archive(
            filename=self.record, extract_dir=self.record_unzip, format="zip"
        )
        if not self.record_unzip.exists():
            logger.fatal(f"Could not find the unzipped directory {self.record_unzip}.")
            raise NotADirectoryError

        matching_dirs = list(self.record_unzip.glob("decrippter2_ml_json_data"))
        if not matching_dirs:
            logger.fatal(
                f"Could not determine data storage location in downloaded directory."
            )
            raise RuntimeError

        subdir = matching_dirs[0]

        shutil.move(
            src=self.record_unzip.joinpath(subdir).joinpath("mite_data/data").resolve(),
            dst=self.location.resolve(),
        )

        shutil.move(
            src=self.record_unzip.joinpath(subdir)
            .joinpath("mite_data/fasta")
            .resolve(),
            dst=self.location.resolve(),
        )

        os.remove(self.record)
        shutil.rmtree(self.record_unzip)