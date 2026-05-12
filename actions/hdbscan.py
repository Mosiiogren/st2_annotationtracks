import hdbscan
import pandas as pd
import numpy as np

from pathlib import Path
from st2common.runners.base_action import Action


class ClusteringData(Action):

    def run(
        self,
        variantfile: str,
        genedata: str,
        regulatorydata: str,
        exondata: str,
        outputfileclusters: str,
    ) -> tuple[bool, str]:

        df_variant = pd.read_csv(variantfile)
        df_variant = df_variant.rename(
            columns={
                "end_chrom": "chromosomeEND",
                "position": "start",
                "sub_category": "Name",
            }
        )
        df_gene = pd.read_json(genedata, orient="records")
        df_regulatory = pd.read_json(regulatorydata, orient="records")
        df_exon = pd.read_json(exondata, orient="records")

    def one_hot_encoding(self, df_variant, df_gene, df_regulatory, df_exon, chromosome):
        index = 0
        rows, _ = df_variant.shape

        metris_metadata = df_variant[["start", "end"]].to_numpy()

        # Take out only the genes in the chromosome in questions
        df_gene = df_gene[df_gene["chromosome"] == chromosome]
        df_regulatory = df_regulatory[df_regulatory["chromosome"] == chromosome]

        # sv_types = df_variant["Name"].unique()
        # dict_sv_types = {value: index for index, value in enumerate(sv_types)}
        # metris_sv_types = np.zeros((self.rows, len(sv_types)))

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

            # Add the SV type
            # sv_type_pos = dict_sv_types[row.Name]
            # metris_sv_types[self.index, sv_type_pos] = 1

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
        print(
            f"One hotencoding matrix: {finalmetadata}, Size: {np.size(finalmetadata, 1)}"
        )

        # Remove all columns that only consists of zero data (to reduce the dimensionality)
        finalmetadata = finalmetadata[:, ~np.all(finalmetadata == 0, axis=0)]
        print(
            f"One hotencoding matrix2: {finalmetadata}, Size: {np.size(finalmetadata, 1)}"
        )

        return finalmetadata
