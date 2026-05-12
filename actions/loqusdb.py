import pandas as pd
import numpy as np

from pathlib import Path
from st2common.runners.base_action import Action


class LoqusDB(Action):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.df_clusters = pd.DataFrame(
            columns=[
                "cluster_number",
                "category",
                "chromosome",
                "chromosomeEND",
                "start",
                "end",
                "score",
                "Name",
            ]
        )

        self.number_of_clusters = 0

        self.SV_range_large = 2000
        self.SV_range_small = 10

        self.test_small = [5, 10, 30, 50, 100, 150, 250, 500]
        self.test_small2 = [100, 200, 300, 500, 1000, 1500, 2500, 5000]
        self.test_large = [1000, 2000, 3000, 5000, 10000, 15000, 25000, 50000]
        self.test_result = []

    def run(
        self,
        variantfile: str,
    ) -> tuple[bool, str]:

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

        # count = df_variant["Name"].value_counts()

        # For BND -> Add true/false whether the start and end chromosomes are the same
        df_variant["interchromosomal"] = np.where(
            (df_variant["chromosome"] == df_variant["chromosomeEND"]), "True", "False"
        )

        df_variant = df_variant[df_variant["chromosome"] == "1"]
        print(len(df_variant))

        # Testing
        for i in range(len(self.test_small)):
            self.df_clusters.drop(self.df_clusters.index, inplace=True)

            for SV in df_variant.itertuples(index=False):

                # SV_range = self.check_SV_length(SV.length, SV.Name)
                #### TESTIN ####
                if SV.length > 10000:
                    SV_range = self.test_large[i]
                else:
                    SV_range = abs(int(self.test_small[i] * SV.length))
                    # SV_range = self.test_small2[i]
                #### TESTING ####

                df_matches = self.find_matching_clusters(SV, SV_range)

                if df_matches.empty:
                    self.addcluster(SV)

                else:
                    matchingclusters, _ = df_matches.shape
                    if matchingclusters > 1:
                        match, cluster_number = self.find_best_macthing_cluster(
                            df_matches, SV, SV_range
                        )

                        if match:
                            df_matches = df_matches[
                                (df_matches["cluster_number"] == cluster_number)
                            ]
                            self.update_cluster(self.get_index(df_matches), SV)
                        else:
                            self.addcluster(SV)

                    else:
                        self.update_cluster(self.get_index(df_matches), SV)

            self.test_result.append(len(self.df_clusters))
        print(self.test_result)

        # self.df_clusters.to_json(outputfileclusters, orient="records")

        return (True, f"Done")

    def addcluster(self, SV: tuple):
        """
        Function that creates a new cluster and adds a SV to it
        """

        self.number_of_clusters += 1
        self.df_clusters.loc[len(self.df_clusters)] = [
            self.number_of_clusters,
            SV.category,
            SV.chromosome,
            SV.chromosomeEND,
            float(SV.start),
            float(SV.end),
            1,
            SV.Name,
        ]

    def update_cluster(self, index: int, SV: tuple):
        """
        Functiont that updates an existing cluster with the new SVs values
        """

        updates = {
            "start": (self.df_clusters.at[index, "start"] + SV.start) / 2,
            "end": (self.df_clusters.at[index, "end"] + SV.end) / 2,
            "score": self.df_clusters.at[index, "score"] + 1,
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
            & (self.df_clusters["category"] == SV.category)
            & (self.df_clusters["Name"] == SV.Name)
            & (self.df_clusters["start"] + SV_range >= SV.start)
            & (self.df_clusters["start"] - SV_range <= SV.start)
            & (self.df_clusters["end"] + SV_range >= SV.end)
            & (self.df_clusters["end"] - SV_range <= SV.end)
        ]

        return df_matches

    def find_best_macthing_cluster(
        self, df_matches: pd.DataFrame, SV: tuple, SV_range: int
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
        else:
            return int(self.SV_range_small * length)
