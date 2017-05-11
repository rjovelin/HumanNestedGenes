# -*- coding: utf-8 -*-
"""
Created on Thu May 11 13:06:34 2017

@author: RJovelin
"""

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
import json
import random
import copy
import sys
import os
import math
import numpy as np
from scipy import stats
from HsaNestedGenes import *


# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
           'Mouse', 'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
           'Cat', 'Opossum', 'Platypus', 'Shrew']

# make a parallel list of json files with nested genes
NestedFiles = [i + 'NestedGenes.json' for i in Species]
PiggybackFiles = [i + 'PiggyBackGenes.json' for i in Species]
ConvergentFiles = [i + 'ConvergentGenes.json' for i in Species]
DivergentFiles = [i + 'DivergentGenes.json' for i in Species]



def LoadDicts(L):
    Dicts = []
    for i in range(len(L)):
        # load dictionary of overlapping gene pairs
        json_data = open(L[i])
        overlap = json.load(json_data)
        json_data.close()
        Dicts.append(overlap)
    return Dicts


# make a list of dictionaries
AllNestedGenes = LoadDicts(NestedFiles)
AllPbkGenes = LoadDicts(PiggybackFiles)
AllConGenes = LoadDicts(ConvergentFiles)
AllDivGenes = LoadDicts(DivergentFiles)


# get nested pairs 
NestedPairs = [GetHostNestedPairs(AllNestedGenes[i]) for i in range(len(AllNestedGenes))]
PgkPairs = [GetHostNestedPairs(AllPbkGenes[i]) for i in range(len(AllPbkGenes))]
ConPairs = [GetHostNestedPairs(AllConGenes[i]) for i in range(len(AllConGenes))]
DivPairs = [GetHostNestedPairs(AllDivGenes[i]) for i in range(len(AllDivGenes))]

for i in [NestedPairs, PgkPairs, ConPairs, DivPairs]:
    print(list(map(lambda x: len(x), i)))





