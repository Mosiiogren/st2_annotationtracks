import re
import io
import datetime
import requests

import pandas as pd
import numpy as np

from gzip import decompress
from pathlib import Path
from st2common.runners.base_action import Action


class RegulatoryData(Action):

    def run(
        self,
        regulatoryfileurl: list[str],
        outputfileregulatory: str,
        columns: list[str],
        attributes: list[str],
        outputfileregulatorytracks: str,
        filename: str,
    ) -> tuple[bool, str]:

        df = pd.DataFrame()
        for url in regulatoryfileurl:
            succeded, results = self.get_data(url, columns, attributes)
            if not succeded:
                return (False, results)

            df = pd.concat([df, results])

        df = self.add_comments(df)
        self.save_to_file(df, outputfileregulatory)
        filename = self.create_annotationtrack_files(
            df, outputfileregulatorytracks, filename
        )

        return (True, filename)

    def get_data(
        self, url: str, columns: list[str], attributes: list[str]
    ) -> tuple[bool, pd.DataFrame | str]:
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

        df.columns = columns

        df = self.get_attributes(
            df=df,
            attributes=attributes,
            last_column=columns[-1],
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
                    else "-"
                )
            )
            df[attribute] = df[attribute].astype("str")
            df.loc[df[attribute].isna(), attribute] = "-"

        df.drop(last_column, axis=1, inplace=True)

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """
        df = df.astype(str)

        df["comments"] = (
            "Ensembl"
            + ";"
            + "Regulatory element: "
            + df["Type"]
            + ";"
            + "ID: "
            + df["ID"]
            + ";"
            + "Extended start: "
            + df["extended_start"]
            + ";"
            + "Extended_end: "
            + df["extended_end"]
            + ";"
            + "Gene ID: "
            + df["gene_id"]
            + ";"
            + "Gene name: "
            + df["gene_name"]
            + ";"
            + "Gene biotype: "
            + df["gene_biotype"]
            + ";"
            + "Source: "
            + df["source"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def save_to_file(self, df: pd.DataFrame, outputfile):
        """
        Function for saving data to json file
        """
        df = df[df["Type"] != "open_chromatin_region"]
        df.to_json(outputfile, orient="records")

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str, filename: str
    ) -> list[str]:

        df = df.drop(
            df.columns.difference(
                [
                    "chromosome",
                    "start",
                    "end",
                    "comments",
                ]
            ),
            axis=1,
        )

        df["start"] = df["start"].astype("int")
        df["end"] = df["end"].astype("int")

        df.to_csv(
            (outputfolder + filename),
            sep="\t",
            index=False,
        )

        return [filename]
