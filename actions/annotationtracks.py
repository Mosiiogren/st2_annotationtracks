import pandas as pd
from st2common.runners.base_action import Action

class Annotationtracks (Action):

    def run (self, clusteringdata):
        print ("Hello from Annotationtracks")

        print (clustering.head)
        return True, "Annotationtracks is working!!!"
