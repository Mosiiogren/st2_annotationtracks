import pandas as pd
from st2common.runners.base_action import Action

class VariantData (Action):

    def run (self, variantfile):
        print ("Helloooooo")

        df = pd.read_json(variantfile)
        print (df.head))
        return True, "It is Working!!!"
