import re
import io
import requests
import pandas as pd
import numpy as np

from gzip import decompress
from pathlib import Path
from st2common.runners.base_action import Action


class GeneData(Action):
    """
    Class for selecting MANE status or Ensemble tagged transcripts/genes and their corresponding exons
    """

    def run(
        self,
        MANEfileurl: str,
        genefileurl: str,
        outputfilegene: str,
        outputfileexon: str,
    ) -> tuple[bool, str]:

        if (Path(outputfilegene).exists()) & (Path(outputfileexon).exists()):
            return (True, "Gene file and Exon file already exists!")

        df_MANE_data, df_ENSEMBL_data = self.get_data(MANEfileurl, genefileurl)
        df_MANE_data = self.filter_MANE_data(df_MANE_data)
        df_ENSEMBL_data = self.filter_gene_data(df_ENSEMBL_data)

        df_ENSEMBL_genes = df_ENSEMBL_data[(df_ENSEMBL_data["feature"] == "gene")]
        df_transcripts = df_ENSEMBL_data[(df_ENSEMBL_data["feature"] == "transcript")]
        df_exons = df_ENSEMBL_data[(df_ENSEMBL_data["feature"] == "exon")]

        # Create genedataframe consisting of MANE status and Ensembl_canonical genes
        df_genes = self.create_MANE_genes_dataframe(
            df_transcripts, df_MANE_data, df_ENSEMBL_genes
        )

        # Create genedataframe consisting on only exons from MANE status transcripts
        df_exons = self.crete_MANE_exons_dataframe(df_exons, df_genes)

        df_genes.to_json(outputfilegene, orient="records")
        df_exons.to_json(outputfileexon, orient="records")

        return (
            True,
            f"Saved gene data to {outputfilegene} and exon data to {outputfileexon}",
        )

    def get_data(self, MANEurl: str, geneurl: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Function for retrieving data from request
        """
        response_MANE = self.get_response(MANEurl)
        df_MANE = pd.read_csv(
            io.StringIO(response_MANE.decode()), sep="\t", na_values="\\N", dtype=str
        )

        response_gene = self.get_response(geneurl)
        df_gene = pd.read_csv(
            io.StringIO(response_gene.decode()),
            sep="\t",
            header=None,
            comment="#",
            dtype=str,
        )

        df_gene.columns = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
        df_gene = self.get_attributes(
            df=df_gene,
            attributes=[
                "gene_id",
                "gene_name",
                "gene_biotype",
                "exon_id",
                "exon_number",
                "transcript_id",
                "tag",
            ],
            last_column="attribute",
        )

        return df_MANE, df_gene

    def get_response(self, url: str) -> str:
        """
        Function for retrieving the response from a request
        """

        response = requests.get(url)
        response = decompress(response.content)

        return response

    def get_attributes(
        self, df: pd.DataFrame, attributes: list, last_column
    ) -> pd.DataFrame:
        """
        Function for extracting key value data
        """

        for attribute in attributes:
            df[attribute] = df[last_column].apply(
                lambda x: (
                    re.findall(rf'{attribute} "([^"]*)"', x)[0]
                    if rf'{attribute} "' in x
                    else np.nan
                )
            )

        df = df.drop(last_column, axis=1)

        return df

    def filter_MANE_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for filtering the data into the correct format
        """

        # Filter the chromosome column -> Looks like: NC_000019.10 (We want 19)
        df = df[df["GRCh38_chr"].str.contains("NT|NW") == False].copy()
        df.loc[:, "GRCh38_chr"] = df["GRCh38_chr"].str.split(".").str[0]
        df.loc[:, "GRCh38_chr"] = df["GRCh38_chr"].str.split("_").str[1]
        # df["GRCh38_chr"] = df["GRCh38_chr"].astype("int").astype("str")
        df.loc[df["GRCh38_chr"] == "23", "GRCh38_chr"] = "X"
        df.loc[df["GRCh38_chr"] == "24", "GRCh38_chr"] = "Y"

        # Remove the version from Ensemble ID and Ensembl_nuc
        df.loc[:, "Ensembl_Gene"] = df["Ensembl_Gene"].str.split(".").str[0]
        df.loc[:, "Ensembl_nuc"] = df["Ensembl_nuc"].str.split(".").str[0]
        df = df.rename(
            columns={"Ensembl_nuc": "transcript_id", "Ensembl_Gene": "gene_id"}
        )

        df = df.drop(
            df.columns.difference(
                [
                    "gene_id",
                    "symbol",
                    "transcript_id",
                    "MANE_status",
                    "GRCh38_chr",
                    "chr_start",
                    "chr_end",
                ]
            ),
            axis=1,
        )

        df["chr_start"] = df["chr_start"].astype("int")
        df["chr_end"] = df["chr_end"].astype("int")

        return df

    def filter_gene_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for filtering the data into the correct format
        """

        # Remove version from chromosome
        df.loc[:, "seqname"] = df["seqname"].str.split(".").str[0]
        df.rename(columns={"seqname": "chromosome"}, inplace=True)

        df["start"] = df["start"].astype("int")
        df["end"] = df["end"].astype("int")

        df = df[df["chromosome"].str.contains("Un|EBV|random|M|KI|GL") == False]
        df = df.drop(["score", "frame", "source"], axis=1)

        return df

    def create_MANE_genes_dataframe(
        self,
        df_transcripts: pd.DataFrame,
        df_MANE: pd.DataFrame,
        df_genes: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Function for slecting MANE status and Ensemble tagged transcripts
        """

        # Merge MANE and ENSEMBLE dataframes based on transcript id
        df_transcripts_MANE_Ensemble = pd.merge(
            df_transcripts, df_MANE, on="transcript_id", how="outer"
        )

        # If there is a MANE status on the trasncript, change the start and end positons to the ones in the MANE dataframe
        df_transcripts_MANE_Ensemble.loc[
            df_transcripts_MANE_Ensemble["MANE_status"].notna(), "start"
        ] = df_transcripts_MANE_Ensemble["chr_start"]
        df_transcripts_MANE_Ensemble.loc[
            df_transcripts_MANE_Ensemble["MANE_status"].notna(), "end"
        ] = df_transcripts_MANE_Ensemble["chr_end"]

        # If MANE status is nan but there is a tag corresponding to the transcript id, add the tag as MANE status
        df_transcripts_MANE_Ensemble.loc[
            (
                (df_transcripts_MANE_Ensemble["tag"] == "Ensembl_canonical")
                | (df_transcripts_MANE_Ensemble["tag"] == "gencode_primary")
                | (df_transcripts_MANE_Ensemble["tag"] == "gencode_basic")
            )
            & (df_transcripts_MANE_Ensemble["MANE_status"].isna()),
            "MANE_status",
        ] = df_transcripts_MANE_Ensemble["tag"]

        # If the chromosme or gene name is gone for the trancript, add the corresponding chromosme from the MANE dataframe instead
        df_transcripts_MANE_Ensemble.loc[
            (df_transcripts_MANE_Ensemble["chromosome"].isna()), "chromosome"
        ] = df_transcripts_MANE_Ensemble["GRCh38_chr"]
        df_transcripts_MANE_Ensemble.loc[
            (df_transcripts_MANE_Ensemble["symbol"].isna())
            & (df_transcripts_MANE_Ensemble["gene_name"].notna()),
            "symbol",
        ] = df_transcripts_MANE_Ensemble["gene_name"]

        # Remove all transcripts that don't have a MANE status
        df_transcripts_MANE_Ensemble = df_transcripts_MANE_Ensemble.drop(
            df_transcripts_MANE_Ensemble[
                df_transcripts_MANE_Ensemble["MANE_status"].isna()
            ].index
        )

        # There are duplicates where different transcripts_id has been taken out for the same gene, these needs to be removed
        df_transcripts_MANE_Ensemble["MANE_status"] = pd.Categorical(
            df_transcripts_MANE_Ensemble["MANE_status"],
            categories=[
                "MANE Plus Clinical",
                "MANE Select",
                "gencode_primary",
                "Ensembl_canonical",
                "gencode_basic",
                "-",
            ],
            ordered=True,
        )
        df_transcripts_MANE_Ensemble.sort_values("MANE_status", inplace=True)
        df_transcripts_MANE_Ensemble.drop_duplicates(
            subset=["gene_id_x"], keep="first", inplace=True
        )

        # Retrieve all genes that did not conatin a MANE status
        # This should be removed
        df_genes_Ensemble = df_genes[
            ~df_genes.gene_id.isin(df_transcripts_MANE_Ensemble.gene_id_x)
        ]

        df_genes_Ensemble = df_genes_Ensemble.drop(
            ["feature", "exon_id", "exon_number", "tag"], axis=1
        )
        df_transcripts_MANE_Ensemble = df_transcripts_MANE_Ensemble.drop(
            [
                "gene_id_y",
                "feature",
                "gene_name",
                "GRCh38_chr",
                "tag",
                "chr_start",
                "chr_end",
                "exon_id",
                "exon_number",
            ],
            axis=1,
        )
        df_transcripts_MANE_Ensemble.rename(
            columns={"gene_id_x": "gene_id", "symbol": "gene_name"}, inplace=True
        )

        df_all_genes = pd.concat([df_transcripts_MANE_Ensemble, df_genes_Ensemble])

        df_all_genes = df_all_genes[(df_all_genes["gene_biotype"] != "artifact")]

        return df_all_genes

    def crete_MANE_exons_dataframe(self, df_exons, df_genes) -> pd.DataFrame:
        """
        Function for slecting the exons belonging to the MANE status and Ensemble tagged transcripts
        """

        df_exons = df_exons.drop(
            [
                "chromosome",
                "feature",
                "gene_name",
                "gene_biotype",
                "tag",
                "strand",
                "gene_id",
            ],
            axis=1,
        )
        df_genes = df_genes.drop(["start", "end", "strand", "gene_biotype"], axis=1)

        # Merge genes and exons based on transcriptid
        df_exons_transcript = pd.merge(
            df_exons, df_genes, on="transcript_id", how="right"
        )

        # Remove exons without MANE status and gene_id
        df_exons_transcript = df_exons_transcript.drop(
            df_exons_transcript[df_exons_transcript["MANE_status"].isna()].index
        )
        df_exons_transcript = df_exons_transcript.drop(
            df_exons_transcript[df_exons_transcript["gene_id"].isna()].index
        )

        return df_exons_transcript
