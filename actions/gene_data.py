import pandas as pd
from st2common.runners.base_action import Action

class VariantData (Action):

    def run (self, variantfile, MANEfile):
        print ("Hello from gene_data")

        df_gene = pd.read_csv(
          genefile,
          sep = "\t",
          header=None,
          comment="#",
          dtype = str
        )

        print (df_gene.head)

        df_MANE = pd.read_csv(
          MANEfile,
          sep = "\t",
          na_values="\\N"
          dtype = str
        )

          
        print (df_MANE.head)
        
        return True, "Gene data is working!!!"
