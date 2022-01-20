from utils import *
from parameters import *

class eventRecord:
    def __init__(self):
        self.inside = False #Initializer 

        self.pos = []
        self.dir = []
        self.length = 0

        self.qDepMap = []
        self.chargeMap = []

        self.weights = []

        self.hitMap = []

        self.pointsPCA = [] 
