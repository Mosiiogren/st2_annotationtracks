import pandas as pd
import numpy as np

from st2common.runners.base_action import Action


COLOR = {
    "DEL": "255,0,0",
    "DUP": "0,0,255",
    "INS":  "34, 139, 34",
    "DUP:TANDEM": "0,0,0",
    "INV": "79,28,115",
    "BND": "0, 255, 0",
    "DUP:INV": "200, 200, 200",
    "Unknown": "100, 100, 150"
}

class Annotationtracks (Action):

    def run (self, clusteringdata, outputfolder):
        df = pd.read_json(clusteringdata, orient="records")
        
        # Add interchromosomal information 
        df = self.interchromosomal(df)
        
        # Add comments
        df = self.addcomments(df)

        # Add color
        df = self.color(df)

        # Create Annotationfiles
        filenames = self.createannotationfiles(df, outputfolder)
     
        return (True, filenames)

    
    def interchromosomal (self, df: pd.DataFrame) -> pd.DataFrame:
        # Add column interchromosomal to separate between BND SVs that starts and ends on different chromosomes
        df["interchromosomal"] = np.where((df["chromosome"] == df["chromosomeEND"]), "True", "False")
        
        df_interchromosomal = df[df["interchromosomal"] == "True"]
        df_non_interchromosomal = df[df["interchromosomal"] == "False"]

        # BND cannot be placed in GENS, however one can add a cluster that starts at one chromosome
        # And another cluster that ends at another chromosome to be able to visualize them in Gens 
        df_non_interchromosomal_copy = df_non_interchromosomal.copy()
        df_non_interchromosomal_copy.drop = df_non_interchromosomal_copy.drop("end", axis=1)
        df_non_interchromosomal_copy["end"] = df_non_interchromosomal_copy["start"] + 1

        df_non_interchromosomal = df_non_interchromosomal.drop(["start", "chromosome"], axis=1)
        df_non_interchromosomal.loc[:, "start"] = df_non_interchromosomal["end"] - 1
        df_non_interchromosomal["chromosome"] = df_non_interchromosomal["chromosomeEND"]

        # Add the dataframes together again
        df = pd.concat([df_interchromosomal, df_non_interchromosomal, df_non_interchromosomal_copy])

        return df
    

    def addcomments (self, df: pd.DataFrame) -> pd.DataFrame:
        
        df["comments"] = np.where((df["chromosome"] == df["chromosomeEND"]), 
                    "SCORE: " + df["score"].astype(str) + ";" + "SV TYPE: " + df["Name"].astype(str),
                    "SCORE: " + df["score"].astype(str) + ";" + "SV TYPE: " + df["Name"].astype(str) + ";" + "Start Chromosome: " + df["chromosome"].astype(str) + ":" + df["start"].astype(str) + ";" + "End chromosome: " + df["chromosomeEND"] + ":" + df["end"].astype(str)
                    )
        
        return df


    def color (self, df: pd.DataFrame) -> pd.DataFrame:
        
        # If SVtype not in COLOR, then add to color Unknown
        for SVtype in df["Name"].unique():
            if SVtype.upper() in COLOR.keys():
                df.loc[df["Name"] == SVtype, "color"] = COLOR[SVtype.upper()] 
            else:
                df.loc[df["Name"] == SVtype, "color"] = COLOR["Unknown"]

        return df
    

    def createannotationfiles (self, df:pd.DataFrame, outputfolder:str) -> list:
        
        # Remove unwanted chromsomes which cannot be displaced in Gens
        df = df[df["chromosome"].str.contains("Un|EBV|random|M") == False].copy()
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        
        filenames = []

        for SVtype in df["Name"].unique():
            df_copy = df[df["Name"] == SVtype]

            # Remove unwanted columns
            df_copy = df_copy.drop(df_copy.columns.difference(["chromosome", "start", "end", "color", "comments"]), axis=1)

            df_copy.to_csv((outputfolder + f"annotationtrackfiles_{SVtype}.tsv"), sep = "\t", index=False)

            filenames.append(f"annotationtrackfiles_{SVtype}.tsv")

        return filenames