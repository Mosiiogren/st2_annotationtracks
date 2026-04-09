import pandas as pd
from st2common.runners.base_action import Action

class VariantData (Action):

    def run (self, variantdata, genedata, regulatorydata):
        print ("Hello from clustering")

        print (variantdata.head)
        print (genedata.head)
        print (regulatorydata.head)
        return True, "Clustering is working!!!"
