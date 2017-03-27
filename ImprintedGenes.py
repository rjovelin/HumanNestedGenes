# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:29:11 2017

@author: RJovelin
"""

# use this script to test enrichement of imprinted genes among nested genes



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



# get option from command
Species = sys.argv[1]
assert Species in ['human', 'mouse']

# load dictionaries of overlapping genes
if Species == 'human':
    jsonFile, GFF = 'HumanNestedGenes.json', 'Homo_sapiens.GRCh38.86.gff3'
elif Species == 'mouse':
    jsonFile, GFF = 'MouseNestedGenes.json', 'Mus_musculus.GRCm38.86.gff3'
    
# load dictionary of nested gene pairs
json_data = open(jsonFile)
Nested = json.load(json_data)
json_data.close()

# make pairs of nested genes
NestedPairs = GetHostNestedPairs(Nested)
External, Internal = set(), set()
for pair in NestedPairs:
    External.add(pair[0])
    Internal.add(pair[1])

# get the status of imprinted gene
if Species == 'human':
    Imprinted = ParseImprinted('HumanImprinted.txt')
elif Species == 'mouse':
    Imprinted = ParseImprinted('MouseImprinted.txt')

# map gene names to gene IDs
GeneNames = MapNametoID(GFF)

# make a set of gene IDs that are imprinted
ImprintedIDs = set()
for gene in GeneNames:
    if GeneNames[gene] in Imprinted:
        ImprintedIDs.add(gene)

exttotal = 0
for gene in External:
    if gene in ImprintedIDs:
        exttotal += 1
inttotal = 0
for gene in Internal:
    if gene in ImprintedIDs:
        inttotal += 1

print('ext', len(External), exttotal)
print('int', len(Internal), inttotal)


