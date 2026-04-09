import pandas as pd
from st2common.runners.base_action import Action

class VariantData (Action):

    def run (self, variantfile):
        print ("Hello from gene_data")

        df = pd.read_csv(
          genefile,
          sep = "\t",
          header=None,
          comment="#",
          dtype = str
        )

          
        print (df.head)
        return True, "Gene data is working!!!"
