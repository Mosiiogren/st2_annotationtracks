import hdbscan
import itertools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from st2common.runners.base_action import Action

# colorblindr::Okabelto
COLOR = {
    "DEL": "230, 159, 0",
    "DUP": "86, 180, 233",
    "INS": "0, 158, 115",
    "DUP:TANDEM": "240, 228, 66",
    "INV": "0, 114, 178",
    "BND": "213, 94, 0",
    "DUP:INV": "204, 121, 167",
    "Unknown": "153, 153, 153",
}


class ClusteringData(Action):

    def run(
        self,
        variantfile: str,
        genedata: str,
        regulatorydata: str,
        exondata: str,
        outputfileclusters: str,
        outputfilefigure: str,
        outputfolder: str,
    ) -> tuple[bool, str]:

        finaldf = pd.DataFrame()
        finaldf2 = pd.DataFrame()


        df_variant = pd.read_csv(variantfile)
        df_variant = df_variant.rename(
            columns={
                "end_chrom": "chromosomeEND",
                "position": "start",
                "sub_category": "Name",
            }
        )
        df_variant = df_variant[
            df_variant["chromosome"].str.contains("Un|EBV|random|M|KI|GL") == False
        ]

        df_variant["intrachromosomal"] = np.where(
            (df_variant["chromosome"] == df_variant["chromosomeEND"]), "True", "False"
        )

        df_variant = df_variant[df_variant["intrachromosomal"] == "True"]

        df_variant = self.add_color(df_variant)

        df_gene = pd.read_json(genedata, orient="records")
        df_regulatory = pd.read_json(regulatorydata, orient="records")
        df_exon = pd.read_json(exondata, orient="records")

        for SVtype in df_variant["Name"].unique():
            df_SV = self.splitting(df_variant, "Name", SVtype)

            for chromosome in df_SV["chromosome"].unique():
                df = self.splitting(df_SV, "chromosome", chromosome)

                rows, _ = df.shape
                if not rows < 2:

                    df = self.min_max_scaler_chrom(df)
                    X = self.one_hot_encoding(
                        df, df_gene, df_regulatory, df_exon, chromosome, rows
                    )

                    clustrered_df, best_dict = self.hdbscan(df, X, SVtype, chromosome)
                    best_df = pd.DataFrame(best_dict, index=[0])

                    dfmean = self.groupby(clustrered_df, "mean", "mean", SVtype)

                    finaldf = pd.concat([finaldf, best_df])
                    finaldf2 = pd.concat([finaldf2, dfmean])

        self.create_annotationtrack_files(finaldf2, outputfolder)

        df_filtered = finaldf[finaldf["best validity index"] != -10]
        if not df_filtered.empty:
            print(f"Mean score: {df_filtered['best validity index'].mean(axis = 0)}")

            self.stripplot(
                df_filtered,
                x_column="parameter",
                y_column="best validity index",
                organization="chromosome",
                x_label="",
                y_label="Density score",
                title="Density score for each chromosome",
                figname=outputfilefigure,
            )

    def splitting(self, df: pd.DataFrame, column: str, type: str) -> pd.DataFrame:

        return df[df[column] == type]

    def min_max_scaler_chrom(self, df: pd.DataFrame) -> pd.DataFrame:

        df["start_norm"] = (df["start"] - df["start"].min()) / (
            df["start"].max() - df["start"].min()
        )
        df["end_norm"] = (df["end"].astype(int) - df["end"].astype(int).min()) / (
            df["end"].astype(int).max() - df["end"].astype(int).min()
        )

        return df

    def one_hot_encoding(
        self, df_variant, df_gene, df_regulatory, df_exon, chromosome, rows
    ):
        index = 0

        metris_metadata = df_variant[["start_norm", "end_norm"]].to_numpy()

        # Take out only the genes in the chromosome in questions
        df_gene = df_gene[df_gene["chromosome"] == chromosome]
        df_regulatory = df_regulatory[df_regulatory["chromosome"] == chromosome]

        gene_biotype = df_gene["gene_biotype"].unique()
        dict_gene_biotype = {value: index for index, value in enumerate(gene_biotype)}
        metris_biotype = np.zeros((rows, len(gene_biotype)))

        gene_regulator = df_regulatory["Type"].unique()
        dict_gene_regulator = {
            value: index for index, value in enumerate(gene_regulator)
        }
        metris_regulator = np.zeros((rows, len(gene_regulator)))
        metris_exons = np.zeros((rows, 1))

        for row in df_variant.itertuples(index=True):


            df_matches_genes = df_gene[
                (df_gene["chromosome"] == row.chromosome)
                & (
                    ((df_gene["start"] >= row.start) & (df_gene["start"] <= row.end))
                    | ((df_gene["end"] >= row.start) & (df_gene["end"] <= row.end))
                    | ((df_gene["start"] <= row.start) & (df_gene["end"] >= row.end))
                )
            ]

            df_matches_exons = df_exon[
                (df_exon["chromosome"] == row.chromosome)
                & (
                    ((df_exon["start"] >= row.start) & (df_exon["start"] <= row.end))
                    | ((df_exon["end"] >= row.start) & (df_exon["end"] <= row.end))
                    | ((df_exon["start"] <= row.start) & (df_exon["end"] >= row.end))
                )
            ]

            df_matches_regulatory = df_regulatory[
                (df_regulatory["chromosome"] == row.chromosome)
                & (
                    (
                        (df_regulatory["start"] >= row.start)
                        & (df_regulatory["start"] <= row.end)
                    )
                    | (
                        (df_regulatory["end"] >= row.start)
                        & (df_regulatory["end"] <= row.end)
                    )
                    | (
                        (df_regulatory["start"] <= row.start)
                        & (df_regulatory["end"] >= row.end)
                    )
                )
            ]

            if not df_matches_genes.empty:

                for match in df_matches_genes.itertuples(index=True):
                    biotype = match.gene_biotype
                    biotype_pos = dict_gene_biotype[biotype]
                    metris_biotype[index, biotype_pos] = 1

            if not df_matches_exons.empty:
                metris_exons[index] = 1

            if not df_matches_regulatory.empty:

                for match in df_matches_regulatory.itertuples(index=True):
                    regulator = match.Type
                    regulator_pos = dict_gene_regulator[regulator]
                    metris_regulator[index, regulator_pos] = 1

            index += 1

        finalmetadata = np.concatenate(
            (metris_metadata, metris_biotype, metris_regulator, metris_exons), axis=1
        )
 
        # Remove all columns that only consists of zero data (to reduce the dimensionality)
        finalmetadata = finalmetadata[:, ~np.all(finalmetadata == 0, axis=0)]

        return finalmetadata

    def hdbscan(self, df, X, type, chromosome):
        """
        Function to cluster given SV data
        """

        best_dict = self.getdensitybasedscores(X, type, chromosome)

        cluster_model = hdbscan.HDBSCAN(
            min_cluster_size=best_dict["best_min_cluster"],
            min_samples=best_dict["best_min_sample"],
        ).fit(X)
        df["clusters"] = cluster_model.labels_

        return df, best_dict

    def getdensitybasedscores(self, X, parameter, chromosome):
        min_cluster_size = np.arange(2, 3, step=1)
        min_samples = np.arange(1, 2, step=1)
        combinations = list(itertools.product(min_cluster_size, min_samples))

        relative_validity = []
        all_labels_list = []

        # https://towardsdatascience.com/tuning-with-hdbscan-149865ac2970/
        for i, (min_cluster, min_sample) in enumerate(combinations):
            cluster_model = hdbscan.HDBSCAN(
                min_cluster_size=min_cluster,
                min_samples=min_sample,
                gen_min_span_tree=True,
            ).fit(X)
            labels_set = set(cluster_model.labels_)
            num_cluster = len(labels_set)
            if num_cluster < 2:
                relative_validity.append(-10)
                all_labels_list.append("bad")
            else:
                # Relative validity
                # https://hdbscan.readthedocs.io/en/latest/api.html#id126
                relative_validity.append(cluster_model.relative_validity_)
                all_labels_list.append(cluster_model.labels_)

        best_index = np.argmax(relative_validity)
        best_parameters = combinations[best_index]
        best_labels = all_labels_list[best_index]

        if "bad" in best_labels:
            validity_index = -10
        else:
            validity_index = hdbscan.validity.validity_index(X, best_labels)
            if np.isnan(validity_index):
                # The datapoints in cluster is identical -> raises a warning
                # Short term solution
                validity_index = relative_validity[best_index]

        return {
            "best_min_cluster": best_parameters[0],
            "best_min_sample": best_parameters[1],
            "best validity index": validity_index,
            "parameter": parameter,
            "chromosome": chromosome,
        }

    def groupby(
        self, df: pd.DataFrame, par1: str, par2: str, name: str
    ) -> pd.DataFrame:
        """ "
        Function to group given clustered data
        """
        df_noise = df[df.clusters == -1]
        columns = ["chromosome", "start", "end", "Name", "color"]
        df_noise.drop(df_noise.columns.difference(columns), axis=1, inplace=True)
        df_noise["comments"] = "Belonging to no cluster group"

        df = df.groupby(["clusters"]).agg(
            chromosome=("chromosome", "first"),
            start=("start", par1),
            end=("end", par2),
            color=("color", "first"),
            comments=("clusters", "count"),
        )

        # Remove the clustergroup with the noise
        df = df[df.index != -1]

        df["start"] = df["start"].round(0).astype(int)
        df["end"] = df["end"].round(0).astype(int)
        df["Name"] = name

        # Add the noise cluster group
        final_df = pd.concat([df, df_noise])

        return final_df

    def stripplot(
        self,
        df: pd.DataFrame,
        x_column: str,
        y_column: str,
        organization: str,
        x_label: str,
        y_label: str,
        title: str,
        figname: str,
    ):
        sns.set_theme()

        fig, ax = plt.subplots(figsize=(8.2, 6))

        plot = sns.stripplot(
            data=df, x=x_column, y=y_column, hue=organization, palette="viridis", ax=ax
        )
        ax.legend(bbox_to_anchor=(1.01, 0.5), loc="center left", borderaxespad=0)
        plot.set_xlabel(x_label, fontsize=12)
        plot.set_ylabel(y_label, fontsize=12)
        plot.set_title(title, pad=10, fontsize=14)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(figname, dpi=300, bbox_inches="tight")

    def add_color(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Function for adding a color to respectively SV type
        """

        for SVtype in df["Name"].unique():
            if SVtype.upper() in COLOR.keys():
                df.loc[df["Name"] == SVtype, "color"] = COLOR[SVtype.upper()]
            else:
                df.loc[df["Name"] == SVtype, "color"] = COLOR["Unknown"]

        return df

    def create_annotationtrack_files(
        self, df: pd.DataFrame, outputfolder: str
    ) -> list[str]:
        """
        Function for creating tab separate annotation files to be uploaded into Gens
        One file for each SV type is created
        """
        # Remove unwanted chromsomes which cannot be displaced in Gens
        df = df[df["chromosome"].str.contains("Un|EBV|random|M") == False].copy()

        df = df.astype(str)

        df["start"] = df["start"].astype("float").astype("int")
        df["end"] = df["end"].astype("float").astype("int")

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
                        "color",
                        "comments",
                    ]
                ),
                axis=1,
            )
            df_copy.to_csv(
                (outputfolder + f"HDBSCAN_{SVtype}.tsv"),
                sep="\t",
                index=False,
            )
            filenames.append(f"HDBSCAN_{SVtype}.tsv")

        return filenames
