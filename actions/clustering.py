import pandas as pd
import numpy as np

from pathlib import Path
from st2common.runners.base_action import Action


class ClusteringData(Action):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.df_clusters = pd.DataFrame(columns=["cluster_number",
                                        "chromosome",
                                        "chromosomeEND",
                                        "start",
                                        "end",
                                        "score",
                                        "Name",
                                        "genes"])
    
        self.number_of_clusters = 0
        
        self.SV_range_large = 10000
        self.SV_range_medium = 5000
        self.SV_range_small = 2500
        self.SV_range = 0

        self.minimum_overlapping_genes_percentage = 0.9

   
    def addcluster (self, row:tuple, genes:list):
        """
            Function that add a SV to a new cluster
        """

        # Check which genes are present and add to the list
        self.number_of_clusters += 1

        self.df_clusters.loc[len(self.df_clusters)] = [self.number_of_clusters, row.chromosome, row.chromosomeEND, float(row.start), float(row.end), 1, row.Name, genes]


        
    def checkSVlength (self, length:int, SVtype:str) -> int:
        """
            Function for setting the range based on the length of the SV
        """
        
        if length > 10000:
            self.SV_range = self.SV_range_large
        elif (length > 4000) & (length <= 10000):
            self.SV_range = self.SV_range_medium
        elif SVtype == "BND":
            self.SV_range = self.SV_range_medium
        else:
            self.SV_range = self.SV_range_small
        

    
    def findMatchingClusters (self, row:tuple) -> pd.DataFrame:
        """
            Function for returning the clusters a SV is matching
        """

        df_matches = self.df_clusters[(self.df_clusters["chromosome"] == row.chromosome) & 
                                        (self.df_clusters["chromosomeEND"] == row.chromosomeEND) & 
                                        (self.df_clusters["Name"] == row.Name) & 
                                        (self.df_clusters["start"] + self.SV_range >= row.start) &
                                        (self.df_clusters["start"] - self.SV_range <= row.start) &
                                        (self.df_clusters["end"] + self.SV_range >= row.end) &
                                        (self.df_clusters["end"] - self.SV_range <= row.end)] 

        return df_matches
    
    def getGENEmatches (self, row:tuple, df:pd.DataFrame, column_match:str, column_row:tuple, column_out:str) -> list:

        df_matches = df[(df[column_match] == column_row) &
                                    (((df["start"] >= row.start) &
                                        (df["start"] <= row.end)) |
                                        ((df["end"] >= row.start) &
                                        (df["end"] <= row.end)) |
                                        ((df["start"] <= row.start) &
                                        (df["end"] >= row.end)))
            ]
        
        return df_matches[column_out].values.tolist()
    

    def getGENEmatchesInterchromosomal (self, row:tuple, df:pd.DataFrame, column_out:str) -> list:
        """
            Function that checks overlapping regulatory elements for each chromosome
        """
        df_matches = df[(df["chromosome"] == row.chromosome) &
                                        (df["start"] <= row.start) &
                                        (df["end"] >= row.start)
                        ]
        matches1 = df_matches[column_out].values.tolist()
        
        df_matches = df[(df["chromosome"] == row.chromosomeEND) &
                                        (df["start"] <= row.end) &
                                        (df["end"] >= row.end)
                        ]
        matches2 = df_matches[column_out].values.tolist()

        return (matches1 + matches2)


    def getallregulatoryfactors (self, row:tuple, df_gene:pd.DataFrame, df_regulatory:pd.DataFrame, df_exon:pd.DataFrame) -> list:
        """
            Function that returns the genes the SV is overlapping
        """
        
        if (row.interchromosomal == "False"):
            genes = self.getGENEmatchesInterchromosomal(row, df_gene, "gene_id")
            regulators = self.getGENEmatchesInterchromosomal(row, df_regulatory, "ID")
            
        else:
            genes = self.getGENEmatches(row, df_gene, "chromosome", row.chromosome, "gene_id")
            regulators = self.getGENEmatches(row, df_regulatory, "chromosome", row.chromosome, "ID")
         

        if (len(genes) > 0) & (len(genes) < 3):
            exons = self.getexons(row, genes, df_exon)

            # Check if the SV overlapped any exons
            if len(exons) == 0:
                # Also need to get the intron data
                introns = self.getintrons (row, genes, df_gene, df_exon)
                return (genes + introns + regulators)
            
            else:
                return (genes + exons + regulators)
            
        else: 
            return (genes + regulators)
    

    def getexons (self, row:tuple, genes:list, df_exon:pd.DataFrame) -> list:
        """
            Function that returns the exons the SV is overlapping
        """

        exons = []
        for gene in genes:
            
            exons_matches = self.getGENEmatches(row, df_exon, "gene_id", gene, "exon_id")
            # All exons as a list that matches the SV
            exons.extend(exons_matches)

        
        return exons

    
    def getintrons (self, row:tuple, genes:list, df_gene:pd.DataFrame, df_exon:pd.DataFrame) -> list:

        introns = []
        for gene in genes:
            all_exons = df_exon[(df_exon["gene_id"] == gene)]
            strand_match = df_gene[(df_gene["gene_id"] == gene)]
            strand = strand_match["strand"].values.tolist()
            
            for i in range(1, len(all_exons["exon_number"].unique()) + 1):
        
                if strand[0] == "+":
                    exon_match_start = all_exons [(all_exons["exon_number"] == i) &
                                            (all_exons["end"] < row.start)
                                            ]
                    exon_match_end = all_exons [(all_exons["exon_number"] == (i+1)) &
                                            (all_exons["start"] > row.end)
                                            ]
                    exon_match_before_start = all_exons [(all_exons["exon_number"] == (i)) &
                                            (all_exons["start"] > row.end)
                                            ]
                else:
                    exon_match_start = all_exons [(all_exons["exon_number"] == i) &
                                            (all_exons["start"] > row.end)
                                            ]
                    exon_match_end = all_exons [(all_exons["exon_number"] == (i+1)) &
                                            (all_exons["end"] < row.start)
                                            ]
                    exon_match_before_start = all_exons [(all_exons["exon_number"] == (i)) &
                                            (all_exons["end"] < row.start)
                                            ]
                
                if (len(exon_match_start["exon_number"].values.tolist()) != 0) & (len(exon_match_end["exon_number"].values.tolist()) != 0):
                    introns.append(list(gene + "intron" + exon_match_start["exon_number"].astype(str))[0])
                    break
                    
                elif (len(exon_match_start["exon_number"].values.tolist()) != 0) & (i == (len(all_exons["exon_number"].unique()))):
                    introns.append(list(gene + "intron" + exon_match_start["exon_number"].astype(str))[0])
                    break

                elif (len(exon_match_before_start["exon_number"].values.tolist()) != 0):
                    # There are some genes that don't start with an exon since same genes have included there promoter/enhance etc. in the gene region and some only the transcript
                    introns.append(gene + "intron0")
                    break
                    
        return introns


    def findbestmacth (self, df_matches:pd.DataFrame, row:tuple, genes:list) -> tuple[bool, int]:
        """
            Function that determines the best macthing cluster for a SV
        """
    
        min_distance = self.SV_range*2
        match = False
        cluster_number = 0
       
        for cluster in df_matches.itertuples(index=False):

            # Calculate the distance between the start and end positions
            distance = abs(cluster.start - row.start) + abs(cluster.end - row.end)

            if (not genes) & (not cluster.genes):
                if distance < min_distance:
                    min_distance = distance
                    cluster_number = cluster.cluster_number
                    match = True

            elif (not genes) | (not cluster.genes):
                # No match since there are no overlapping genes
                continue

            else:
                if len(cluster.genes) > len(genes):
                    minimum_overlapping_genes = len(cluster.genes) * self.minimum_overlapping_genes_percentage
                else:
                    minimum_overlapping_genes = len(genes) * self.minimum_overlapping_genes_percentage
                
                # Calculate how many overlapping genes there are
                overlapping_genes = 0
                for gene in genes:
                    if gene in cluster.genes:
                        overlapping_genes += 1
            
                # Check if the genes overlapping is enough to be considered a match
                if overlapping_genes >= minimum_overlapping_genes:
                    
                    # Check if the distance is lower than other candidates
                    if distance < min_distance:
                        min_distance = distance
                        cluster_number = cluster.cluster_number
                        match = True

        return match, cluster_number
    
    def updatecluster (self, index: int, row: tuple, genes: list):
        """
            Functiont that updates an existing cluster with the new SVs values
        """
        # The columns to be updated
        updates = {
            "start" : (self.df_clusters.at[index, "start"] + row.start)/2,
            "end" : (self.df_clusters.at[index, "end"] + row.end)/2,
            "score" : self.df_clusters.at[index, "score"] + 1,
            "genes" : list(set(self.df_clusters.at[index, "genes"] + genes))       
        }

        for column, value in updates.items():
            self.df_clusters.at[index, column] = value
    

    def getindex (self, cluster_match:pd.DataFrame | pd.Series) -> int:
        """
            Function that returns the index of the cluster that macthes as SV
        """

        # Since cluster_match can be both a dataframe and a series we can "squezze" it to get the scalar value
        cluster_number = cluster_match["cluster_number"].squeeze()
        index = self.df_clusters.index[self.df_clusters["cluster_number"] == cluster_number][0]

        return index


    def run(self, variantfile, genedata, regulatorydata, exondata, outputfileclusters):
        
        # REMOVE WHEN DONE!!!
        if (Path(outputfileclusters).exists()):
          return (True, "Cluster file already exists!")

        # Convert json files/arrays to pandas dataframe
        df_variant = pd.DataFrame(variantfile)
        df_gene = pd.read_json(genedata, orient="records")
        df_regulatory = pd.read_json(regulatorydata, orient="records")
        df_exon = pd.read_json(exondata, orient="records")

        # Check if df_variant is empty
            # Check if all is empty?

        # For BND -> Add if the
        df_variant["interchromosomal"] = np.where((df_variant["chromosome"] == df_variant["chromosomeEND"]), "True", "False")

        for row in df_variant.itertuples(index=False):
            
            # Get all regulatoryfactors the SV is overlapping
            genes = self.getallregulatoryfactors (row, df_gene, df_regulatory, df_exon)

            # Set the range based on the lenght of the SV
            self.checkSVlength(row.length, row.Name)

            # Find matching clusters
            df_matches = self.findMatchingClusters(row)

            if df_matches.empty:
                self.addcluster(row, genes)
            
            else:
                # Check how many clusters the SV is matching
                rows,_ = df_matches.shape
                if rows > 1:
                    # More than one match! Check which is the best fitted and has overlapping genes
                    match, cluster_number = self.findbestmacth(df_matches, row, genes)
                    
                    if match:
                        # Take the cluster the row has the minimum distance to
                        df_matches = df_matches[(df_matches["cluster_number"] == cluster_number)]
                        self.updatecluster(self.getindex(df_matches), row, genes)
                    else:
                        self.addcluster(row, genes)

                else:
                    # Check if the genes in the SV are overlapping the once in the cluster
                    if (not genes) & (not df_matches["genes"].values.tolist()[0]):
                        self.updatecluster(self.getindex(df_matches), row, genes)
                        
                    elif (not genes) | (not df_matches["genes"].values.tolist()[0]):
                        self.addcluster(row, genes)

                    else:
                        
                        if len(df_matches["genes"].values.tolist()[0]) > len(genes):
                            minimum_overlapping_genes = len(df_matches["genes"].values.tolist()[0]) * self.minimum_overlapping_genes_percentage
                        else:
                            minimum_overlapping_genes = len(genes) * self.minimum_overlapping_genes_percentage
                        
                        # Check if the majority of the regulatory factors are overlapping
                        overlapping_genes = 0
                        for gene in genes:
                            # Check if a gene is the same in the SV and the cluster
                            if gene in df_matches["genes"].values.tolist()[0]:
                                overlapping_genes += 1
                        
                        if overlapping_genes >= minimum_overlapping_genes:
                            self.updatecluster(self.getindex(df_matches), row, genes)
                        else:
                            self.addcluster(row, genes)

        # Save to a json file
        self.df_clusters.to_json(outputfileclusters, orient="records")
        
        return (True, f"Clusters are stored in {outputfileclusters}")


    



    