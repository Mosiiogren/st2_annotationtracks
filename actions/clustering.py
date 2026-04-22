import pandas as pd
import numpy as np

from pathlib import Path
from st2common.runners.base_action import Action


class ClusteringData(Action):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.df_clusters = pd.DataFrame(
            columns=[
                "cluster_number",
                "chromosome",
                "chromosomeEND",
                "start",
                "end",
                "score",
                "Name",
                "genes",
            ]
        )

        self.number_of_clusters = 0

        self.SV_range_large = 10000
        self.SV_range_medium = 5000
        self.SV_range_small = 2500

        self.minimum_overlapping_genes_percentage = 0.9

    def run(
        self,
        variantfile: str,
        genedata: str,
        regulatorydata: str,
        exondata: str,
        outputfileclusters: str,
    ) -> tuple[bool, str]:

        # REMOVE WHEN DONE!!!
        if Path(outputfileclusters).exists():
            return (True, "Cluster file already exists!")

        df_variant = pd.DataFrame(variantfile)
        df_gene = pd.read_json(genedata, orient="records")
        df_regulatory = pd.read_json(regulatorydata, orient="records")
        df_exon = pd.read_json(exondata, orient="records")

        # For BND -> Add true/false whether the start and end chromosomes are the same
        df_variant["interchromosomal"] = np.where(
            (df_variant["chromosome"] == df_variant["chromosomeEND"]), "True", "False"
        )

        for SV in df_variant.itertuples(index=False):
            genes = self.get_all_regulatory_elements(
                SV, df_gene, df_regulatory, df_exon
            )
            SV_range = self.check_SV_length(SV.length, SV.Name)
            df_matches = self.find_matching_clusters(SV, SV_range)

            if df_matches.empty:
                self.addcluster(SV, genes)

            else:
                matchingclusters, _ = df_matches.shape
                if matchingclusters > 1:
                    match, cluster_number = self.find_best_macthing_cluster(
                        df_matches, SV, genes, SV_range
                    )

                    if match:
                        df_matches = df_matches[
                            (df_matches["cluster_number"] == cluster_number)
                        ]
                        self.update_cluster(self.get_index(df_matches), SV, genes)
                    else:
                        self.addcluster(SV, genes)

                else:
                    if (not genes) & (not df_matches["genes"].values.tolist()[0]):
                        self.update_cluster(self.get_index(df_matches), SV, genes)

                    elif (not genes) | (not df_matches["genes"].values.tolist()[0]):
                        self.addcluster(SV, genes)

                    else:
                        match = self.check_overlapping_genes(
                            genes, df_matches["genes"].values.tolist()[0]
                        )
                        if match:
                            self.update_cluster(self.get_index(df_matches), SV, genes)
                        else:
                            self.addcluster(SV, genes)

        self.df_clusters.to_json(outputfileclusters, orient="records")

        return (True, f"Clusters are stored in {outputfileclusters}")

    def addcluster(self, SV: tuple, genes: list[str]):
        """
        Function that creates a new cluster and adds a SV to it
        """

        self.number_of_clusters += 1
        self.df_clusters.loc[len(self.df_clusters)] = [
            self.number_of_clusters,
            SV.chromosome,
            SV.chromosomeEND,
            float(SV.start),
            float(SV.end),
            1,
            SV.Name,
            genes,
        ]

    def update_cluster(self, index: int, SV: tuple, genes: list[str]):
        """
        Functiont that updates an existing cluster with the new SVs values
        """

        updates = {
            "start": (self.df_clusters.at[index, "start"] + SV.start) / 2,
            "end": (self.df_clusters.at[index, "end"] + SV.end) / 2,
            "score": self.df_clusters.at[index, "score"] + 1,
            "genes": list(set(self.df_clusters.at[index, "genes"] + genes)),
        }

        for column, value in updates.items():
            self.df_clusters.at[index, column] = value

    def find_matching_clusters(self, SV: tuple, SV_range: int) -> pd.DataFrame:
        """
        Function for returning the clusters a SV is matching
        """

        df_matches = self.df_clusters[
            (self.df_clusters["chromosome"] == SV.chromosome)
            & (self.df_clusters["chromosomeEND"] == SV.chromosomeEND)
            & (self.df_clusters["Name"] == SV.Name)
            & (self.df_clusters["start"] + SV_range >= SV.start)
            & (self.df_clusters["start"] - SV_range <= SV.start)
            & (self.df_clusters["end"] + SV_range >= SV.end)
            & (self.df_clusters["end"] - SV_range <= SV.end)
        ]

        return df_matches

    def find_best_macthing_cluster(
        self, df_matches: pd.DataFrame, SV: tuple, genes: list[str], SV_range: int
    ) -> tuple[bool, int]:
        """
        Function that determines the best macthing cluster for a SV
        """

        min_distance = SV_range * 2
        match = False
        cluster_number = 0

        for cluster in df_matches.itertuples(index=False):
            # Calculate the distance between the start and end positions
            distance = abs(cluster.start - SV.start) + abs(cluster.end - SV.end)

            if (not genes) & (not cluster.genes):
                if distance < min_distance:
                    min_distance = distance
                    cluster_number = cluster.cluster_number
                    match = True

            elif (not genes) | (not cluster.genes):
                # No match since there are no overlapping genes
                continue

            else:
                match = self.check_overlapping_genes(genes, cluster.genes)
                if match:
                    # Check if the distance is lower than other candidates
                    if distance < min_distance:
                        min_distance = distance
                        cluster_number = cluster.cluster_number
                        match = True

        return match, cluster_number

    def get_matches_interchromosomal(
        self,
        SV: tuple,
        df: pd.DataFrame,
        column_match: str,
        column_SV: tuple,
        column_out: str,
    ) -> list[str]:
        """
        Function that checks overlapping regulatory elements for a chromosome
        """

        df_matches = df[
            (df[column_match] == column_SV)
            & (
                ((df["start"] >= SV.start) & (df["start"] <= SV.end))
                | ((df["end"] >= SV.start) & (df["end"] <= SV.end))
                | ((df["start"] <= SV.start) & (df["end"] >= SV.end))
            )
        ]

        return df_matches[column_out].values.tolist()

    def get_matches_non_interchromosomal(
        self, SV: tuple, df: pd.DataFrame, column_out: str
    ) -> list[str]:
        """
        Function that checks overlapping regulatory elements for each chromosome
        """

        df_matches = df[
            (df["chromosome"] == SV.chromosome)
            & (df["start"] <= SV.start)
            & (df["end"] >= SV.start)
        ]
        matches1 = df_matches[column_out].values.tolist()

        df_matches = df[
            (df["chromosome"] == SV.chromosomeEND)
            & (df["start"] <= SV.end)
            & (df["end"] >= SV.end)
        ]
        matches2 = df_matches[column_out].values.tolist()

        return matches1 + matches2

    def get_exons(
        self, SV: tuple, genes: list[str], df_exon: pd.DataFrame
    ) -> list[str]:
        """
        Function that returns the exons the SV is overlapping
        """

        exons = []
        for gene in genes:

            exons_matches = self.get_matches_interchromosomal(
                SV, df_exon, "gene_id", gene, "exon_id"
            )
            exons.extend(exons_matches)

        return exons

    def get_introns(
        self, SV: tuple, genes: list[str], df_gene: pd.DataFrame, df_exon: pd.DataFrame
    ) -> list[str]:
        """
        Function for extracting the introns a SV is overlapping
        The function considers the strands direction in order to retrieve the correct intron number
        """
        introns = []
        for gene in genes:
            all_exons = df_exon[(df_exon["gene_id"] == gene)]
            strand_match = df_gene[(df_gene["gene_id"] == gene)]
            strand = strand_match["strand"].values.tolist()

            for i in range(1, len(all_exons["exon_number"].unique()) + 1):

                if strand[0] == "+":
                    exon_match_start = all_exons[
                        (all_exons["exon_number"] == i) & (all_exons["end"] < SV.start)
                    ]
                    exon_match_end = all_exons[
                        (all_exons["exon_number"] == (i + 1))
                        & (all_exons["start"] > SV.end)
                    ]
                    exon_match_before_start = all_exons[
                        (all_exons["exon_number"] == (i))
                        & (all_exons["start"] > SV.end)
                    ]
                else:
                    exon_match_start = all_exons[
                        (all_exons["exon_number"] == i) & (all_exons["start"] > SV.end)
                    ]
                    exon_match_end = all_exons[
                        (all_exons["exon_number"] == (i + 1))
                        & (all_exons["end"] < SV.start)
                    ]
                    exon_match_before_start = all_exons[
                        (all_exons["exon_number"] == (i))
                        & (all_exons["end"] < SV.start)
                    ]

                if (len(exon_match_start["exon_number"].values.tolist()) != 0) & (
                    len(exon_match_end["exon_number"].values.tolist()) != 0
                ):
                    introns.append(
                        list(
                            gene
                            + "intron"
                            + exon_match_start["exon_number"].astype(str)
                        )[0]
                    )
                    break

                elif (len(exon_match_start["exon_number"].values.tolist()) != 0) & (
                    i == (len(all_exons["exon_number"].unique()))
                ):
                    introns.append(
                        list(
                            gene
                            + "intron"
                            + exon_match_start["exon_number"].astype(str)
                        )[0]
                    )
                    break

                elif len(exon_match_before_start["exon_number"].values.tolist()) != 0:
                    # There are some genes that don't start with an exon since same genes have included there promoter/enhance etc. in the gene region and some only the transcript
                    introns.append(gene + "intron0")
                    break

        return introns

    def get_all_regulatory_elements(
        self,
        SV: tuple,
        df_gene: pd.DataFrame,
        df_regulatory: pd.DataFrame,
        df_exon: pd.DataFrame,
    ) -> list[str]:
        """
        Function that returns all regulatory elements the SV is overlapping
        """

        if SV.interchromosomal == "False":
            genes = self.get_matches_non_interchromosomal(SV, df_gene, "gene_id")
            regulators = self.get_matches_non_interchromosomal(SV, df_regulatory, "ID")
        else:
            genes = self.get_matches_interchromosomal(
                SV, df_gene, "chromosome", SV.chromosome, "gene_id"
            )
            regulators = self.get_matches_interchromosomal(
                SV, df_regulatory, "chromosome", SV.chromosome, "ID"
            )

        if (len(genes) > 0) & (len(genes) < 4):
            exons = self.get_exons(SV, genes, df_exon)

            if len(exons) == 0:
                introns = self.get_introns(SV, genes, df_gene, df_exon)
                return genes + introns + regulators

            else:
                return genes + exons + regulators

        else:
            return genes + regulators

    def check_overlapping_genes(
        self, SVgenes: list[str], clustergenes: list[str]
    ) -> bool:
        """
        Function for checking if a SV and a cluster have enough overlapping genes to be considered a match
        """

        if len(clustergenes) > len(SVgenes):
            minimum_overlapping_genes = (
                len(clustergenes) * self.minimum_overlapping_genes_percentage
            )
        else:
            minimum_overlapping_genes = (
                len(SVgenes) * self.minimum_overlapping_genes_percentage
            )

        overlapping_genes = 0
        for gene in SVgenes:
            if gene in clustergenes:
                overlapping_genes += 1

        if overlapping_genes >= minimum_overlapping_genes:
            return True
        else:
            return False

    def get_index(self, cluster_match: pd.DataFrame | pd.Series) -> int:
        """
        Function that returns the index of the cluster that macthes as SV
        """

        # Since cluster_match can be both a dataframe and a series we can "squezze" it to get the scalar value
        cluster_number = cluster_match["cluster_number"].squeeze()
        index = self.df_clusters.index[
            self.df_clusters["cluster_number"] == cluster_number
        ][0]

        return index

    def check_SV_length(self, length: int, SVtype: str) -> bool:
        """
        Function for setting the range based on the length of the SV
        """

        if length > 10000:
            return self.SV_range_large
        elif ((length > 4000) & (length <= 10000)) | (SVtype == "BND"):
            return self.SV_range_medium
        else:
            return self.SV_range_small
