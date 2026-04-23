import re
import io
import requests
import pandas as pd
import numpy as np

from gzip import decompress
from pathlib import Path

from st2common.runners.base_action import Action


class RegulatoryData(Action):

    def run(
        self, regulatoryfileurl: str, outputfileregulatory: str
    ) -> tuple[bool, str]:

        if Path(outputfileregulatory).exists():
            return (True, "Regulator file already exists!")

        succeded, results = self.get_data(regulatoryfileurl)
        if not succeded:
            return (False, results)
        df = self.filter_regulatoryfactor_data(results)

        df.to_json(outputfileregulatory, orient="records")

        return (True, f"Saved regulatory data to {outputfileregulatory}")

    def get_data(self, url: str) -> tuple[bool, pd.DataFrame | str]:
        """
        Function for retrieving data from request
        """
        try:
            response = requests.get(url)
            response.raise_for_status()
            try:
                response = decompress(response.content)
            except:
                return (False, f"File {url} is not a gzip-compressed file")
        except:
            return (False, f"Problem with request from {url}")

        df = pd.read_csv(
            io.StringIO(response.decode()),
            sep="\t",
            header=None,
            comment="#",
            dtype=str,
        )

        df.columns = columns = [
            "chromosome",
            "source",
            "Type",
            "start",
            "end",
            "col6",
            "col7",
            "col8",
            "col9",
        ]
        df = self.get_attributes(
            df=df, attributes=["ID", "gene_id"], last_column="col9"
        )

        return (True, df)

    def get_attributes(
        self, df: pd.DataFrame, attributes: list[str], last_column: str
    ) -> pd.DataFrame:
        """
        Function for extracting key value data
        """

        for attribute in attributes:
            df[attribute] = df[last_column].apply(
                lambda x: (
                    re.findall(rf"{attribute}=([^;]*)", x)[0]
                    if f"{attribute}=" in x
                    else np.nan
                )
            )

        df.drop(last_column, axis=1, inplace=True)

        return df

    def filter_regulatoryfactor_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for filtering the data into the correct format
        """

        df["start"] = df["start"].astype("int")
        df["end"] = df["end"].astype("int")

        df = df[df["Type"] != "open_chromatin_region"]

        return df
