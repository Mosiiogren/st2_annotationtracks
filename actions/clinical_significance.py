import re
import io
import requests
import datetime

import pandas as pd
import numpy as np
from gzip import decompress
from pathlib import Path

from st2common.runners.base_action import Action


class ClinicalSignificance(Action):

    def run(
        self,
        clinicalsignificancefileurl: list,
        outputclinicalsignificance: str,
        columns: list,
        attributes: list,
    ) -> tuple[bool, str]:

        finaldf = pd.DataFrame()
        for url in clinicalsignificancefileurl:
            succeded, results = self.get_data(url, columns, attributes)
            if not succeded:
                return (False, results)
            if results.empty:
                continue

            df = self.filter_data(results, attributes)
            df = self.add_comments(df, attributes)

            finaldf = pd.concat([finaldf, df])

        filename = self.create_annotationtrack_files(
            finaldf, outputclinicalsignificance
        )

        return (True, filename)

    def get_data(
        self, url: str, columns: list, attributes: list
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

        # Check how many rows after comments are over, if zero then skip
        try:
            df = pd.read_csv(
                io.StringIO(response.decode()),
                sep="\t",
                header=None,
                comment="#",
                dtype=str,
            )
        except pd.errors.EmptyDataError:
            # Return an empty dataframe
            return (True, pd.DataFrame())
        except:
            return (False, f"Problem with creating dataframe of {url}")

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
                    else np.nan
                )
            )

        df.drop(last_column, axis=1, inplace=True)

        return df

    def filter_data(self, df: pd.DataFrame, attributes: list) -> pd.DataFrame:
        """
        Function for filtering the dataframe
        """
        # Filter the chromosome column -> Looks like: NC_000019.10 (We want 19)
        df = df[df["chromosome"].str.contains("NC") == True].copy()
        df.loc[:, "chromosome"] = df["chromosome"].str.split(".").str[0]
        df.loc[:, "chromosome"] = df["chromosome"].str.split("_").str[1]
        df["chromosome"] = df["chromosome"].astype("int")

        # Remove any additional chromosomes
        df = df[(df["chromosome"] < 25)]

        # Change chromosome number 23 and 24 to X and Y respectively
        df["chromosome"] = df["chromosome"].astype("str")
        df.loc[df["chromosome"] == "23", "chromosome"] = "X"
        df.loc[df["chromosome"] == "24", "chromosome"] = "Y"

        for attribute in attributes:
            df[attribute] = df[attribute].astype("str")
            df.loc[df[attribute].isna(), attribute] = "-"

        return df

    def add_comments(self, df: pd.DataFrame, attributes: list) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """
        # df["comments"] = "SV TYPE: " + df["Name"].astype(str) + ";"
        # for attribute in attributes:
        #     df["comments"].map(
        #         lambda lst: lst.append(attribute + df[attribute].astype(str) + ";")
        #     )
        df["comments"] = (
            "Clinical significance: "
            + df["clinical_int"].astype(str)
            + ";"
            + "SV TYPE: "
            + df["Name"].astype(str)
            + ";"
            + "Phenotype: "
            + df["phenotype"].astype(str)
            + ";"
            + "Phenotype ID: "
            + df["phenotype_id"].astype(str)
            + ";"
            + "Consequences: "
            + df["consequence"].astype(str)
            + ";"
            + "Validated: "
            + df["validated"].astype(str)
            + ";"
            + "Copy number: "
            + df["copy_number"].astype(str)
            + ";"
            + "Outer and/or inner range of start position: "
            + df["Start_range"].astype(str)
            + ";"
            + "Outer and/or inner range of end position: "
            + df["End_range"].astype(str)
            + ";"
            + "Zygosity: "
            + df["zygosity"].astype(str)
            + ";"
            + "Web link to the variant: "
            + df["Dbxref"].astype(str)
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str
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

        filename = "ClinicalSignificance"

        df.to_csv(
            (outputfolder + filename),
            sep="\t",
            index=False,
        )

        return [filename]


# https://ftp.ensembl.org/pub/release-113/variation/MaveDB/
