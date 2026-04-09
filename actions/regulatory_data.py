import pandas as pd
from st2common.runners.base_action import Action

class VariantData (Action):

    def run (self, regulatoryfile):
        print ("Hello from regulatory_data")

        df = pd.read_csv(
          regulatoryfile,
          sep = "\t",
          header=None,
          comment="#",
          dtype = str
        )

          
        print (df.head)
        return True, "Regulatory data is working!!!"
