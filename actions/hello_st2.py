
from st2common.runners.base_action import Action

class HelloStackStorm (Action):
    def run (self):
        print ("Helloooooo")
        return True, "It is Working!!!"