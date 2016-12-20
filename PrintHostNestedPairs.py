# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 08:16:18 2016

@author: Richard
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


# load dictionaries of host and nested genes 
with open('HumanNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
# load dictionaries with overlapping genes
with open('HumanOverlappingGenes.json') as human_json_data:
    HumanOverlappingGenes = json.load(human_json_data)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'

# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF]
# make a list of species names
SpeciesNames = ['human', 'chimp', 'gorilla']


# loop over GFF files
for i in range(len(GFFs)):
    print(GFFs[i][:GFFs[i].index('.')], SpeciesNames[i])
    # get expression profile of the species genes
    SpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[i])
    # remove genes wuthout expression
    SpExpression = RemoveGenesLackingExpression(SpExpression)
    # get relative expression
    SpExpression = TransformRelativeExpression(SpExpression)
    # make a list of host-nested gene pairs
    SpHostNestedPairs = GetHostNestedPairs(HostGenes[i])
    # make a list of overlapping gene pairs
    SpOverlappingPairs = GetHostNestedPairs(Overlapping[i])
    # assert each gene pair is found only once in each level of genomic organization
    for pair in SpHostNestedPairs:
        assert SpHostNestedPairs.count(pair) == 1
    for pair in SpOverlappingPairs:
        assert SpOverlappingPairs.count(pair) == 1
    
    print('total number of host-nested pairs', SpeciesNames[i], len(SpHostNestedPairs))
    print('total number of overlapping pairs', SpeciesNames[i], len(SpOverlappingPairs))
    # remove gene pairs with genes lacking expression
    SpHostNestedPairs = FilterGenePairsWithoutExpression(SpHostNestedPairs, SpExpression)
    SpOverlappingPairs = FilterGenePairsWithoutExpression(SpOverlappingPairs, SpExpression)
    print('number of host-nested pairs with expression', SpeciesNames[i], len(SpHostNestedPairs))
    print('number of overlapping pairs with expression', SpeciesNames[i], len(SpOverlappingPairs))
    