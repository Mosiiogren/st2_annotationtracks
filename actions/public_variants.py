import io
import requests
import datetime
import time

import pandas as pd
from pathlib import Path
from gzip import decompress

from st2common.runners.base_action import Action


class PublicStructuralVariantsData(Action):

    def run(
        self,
        urls: list[str],
        outputfolder: str,
        columns: list[str],
        filename: str,
    ) -> tuple[bool, str]:

        finaldf = pd.DataFrame()
        for url in urls:
            succeded, results = self.get_data(url, columns)
            if not succeded:
                return (False, results)
            if results.empty:
                continue

            df = self.filter_data(results)
            df = self.add_clinical_significance(df, url)

            finaldf = pd.concat([finaldf, df])

            time.sleep(5)

        df = self.add_comments(finaldf)
        filename = self.create_annotationtrack_files(df, outputfolder, filename)

        return (True, filename)

    def get_data(self, url: str, columns: list[str]) -> tuple[bool, pd.DataFrame | str]:
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

        try:
            df = pd.read_csv(
                io.StringIO(response.decode()),
                sep="\t",
                header=None,
                comment="#",
                dtype=str,
            )
        except pd.errors.EmptyDataError:
            return (True, pd.DataFrame())
        except:
            return (False, f"Problem with creating dataframe of {url}")

        df.columns = columns

        return (True, df)

    def filter_data(self, df: pd.DataFrame) -> pd.DataFrame:

        df.loc[:, "chromosome"] = df["chromosome"].str.split("chr").str[1]
        df["Name"] = df["info"].str.split("_").str[-1]

        return df

    def add_clinical_significance(self, df: pd.DataFrame, url: str) -> pd.DataFrame:

        df["clincal_significance"] = (
            url.rsplit("/")[-1].rsplit(".bed")[0].rsplit(".")[-1]
        )

        return df

    def add_comments(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a column consisting of comments
        Comments are visually separate in Gens by ;
        """

        df = df.astype(str)

        df["comments"] = (
            "Structural variants recieved from dbVar"
            + ";"
            + "SV TYPE: "
            + df["Name"]
            + ";"
            + "Clinical significance: "
            + df["clincal_significance"]
            + ";"
            + "Track created at: "
            + datetime.datetime.now().strftime("%c")
        )

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str, filename: str
    ) -> list[str]:

        filenames = []
        for SVtype in df["Name"].unique():
            df_copy = df[(df["Name"] == SVtype)]

            if df_copy.empty:
                continue

            df_copy = df_copy.drop(
                df_copy.columns.difference(
                    [
                        "chromosome",
                        "start",
                        "end",
                        "comments",
                    ]
                ),
                axis=1,
            )
            df_copy.to_csv(
                (outputfolder + f"{filename}_{SVtype}.tsv"),
                sep="\t",
                index=False,
            )
            filenames.append(f"{filename}_{SVtype}.tsv")

        return filenames


# Contains all information -> Should be better although there will be al lot of SVs displayed
# https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/deletions/GRCh38.nr_deletions.tsv.gz

# Add ACMG genes
# https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/annotation/GRCh38/ACMG_genes.gff.gz
