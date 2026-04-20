import re
import io
import requests
import pandas as pd
import numpy as np

from gzip import decompress
from pathlib import Path

from st2common.runners.base_action import Action




class RegulatoryData (Action):

    def run (self, regulatoryfileurl, outputfileregulatory):

        # Check if file already exists
        if (Path(outputfileregulatory).exists()):
          return (True, "Regulator file already exists!")

        df = self.getdata(regulatoryfileurl)
        df = self.filterREGULATORYdata(df)

        df.to_json(outputfileregulatory, orient="records")
          
        return (True, f"Saved regulatory data to {outputfileregulatory}")
    

    def getdata (self, url:str) -> pd.DataFrame:

        response = requests.get(url)
        response = decompress(response.content)

        df = pd.read_csv(
                io.StringIO(response.decode()),
                sep = "\t",
                header=None,
                comment="#",
                dtype = str
              )

        # Check if the column names has changed

        df.columns = columns = ["chromosome", "source", "Type", "start", "end", "col6", "col7", "col8", "col9"]
        df = self.getattributes (df=df,
                                  attributes=["ID", "gene_id"],
                                  last_column= "col9"
                                  )

        return df


    def getattributes (self, df:pd.DataFrame, attributes:list, last_column) -> pd.DataFrame:
        """
          Function for extracting key value data
        """
        
        for attribute in attributes:
          df[attribute] = df[last_column].apply(lambda x: re.findall(rf'{attribute}=([^;]*)', x)[0] if f'{attribute}=' in x else np.nan)

        df.drop(last_column, axis=1, inplace=True)

        return df

    def filterREGULATORYdata (self, df:pd.DataFrame) -> pd.DataFrame:

        df["chromosome"] = df["chromosome"].astype("string")
        df["start"] = df["start"].astype("int")
        df["end"] = df["end"].astype("int")

        df = df[df["Type"] != "open_chromatin_region"]

        return df
