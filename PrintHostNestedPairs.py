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
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('ChimpHostNestedGenes.json') as chimp_json_data:
    ChimpHostGenes = json.load(chimp_json_data)
with open('GorillaHostNestedGenes.json') as gorilla_json_data:
    GorillaHostGenes = json.load(gorilla_json_data)
with open('OrangOutanHostNestedGenes.json') as orangoutan_json_data:
    OrangOutanHostGenes = json.load(orangoutan_json_data)
with open('MacaqueHostNestedGenes.json') as macaque_json_data:
    MacaqueHostGenes = json.load(macaque_json_data)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
PabGFF = 'Pongo_abelii.PPYG2.86.gff3'
MmlGFF = 'Macaca_mulatta.Mmul_8.0.1.86.gff3' 

# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF, PabGFF, MmlGFF]
# make a list of species names
SpeciesNames = ['human', 'chimp', 'gorilla', 'orangoutan', 'macaque']
# make a list of host:nested genes dictionaries
HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes, OrangOutanHostGenes, MacaqueHostGenes]



# loop over GFF files, find nested and intronic=nested genes in each species 
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
    print('total number of host-nested pairs', SpeciesNames[i], len(SpHostNestedPairs))
    # remove gene pairs with genes lacking expression
    SpHostNestedPairs = FilterGenePairsWithoutExpression(SpHostNestedPairs, SpExpression)
    print('number of host-nested pairs with expression', SpeciesNames[i], len(SpHostNestedPairs))
    