# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:29:11 2017

@author: RJovelin
"""

# use this script to count the number of external and internal genes that are imprinted

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
    jsonFile, GFF = 'HumanNestedGenes.json', 'Homo_sapiens.GRCh38.88.gff3'
elif Species == 'mouse':
    jsonFile, GFF = 'MouseNestedGenes.json', 'Mus_musculus.GRCm38.88.gff3'
    
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
 
# make a dict with gene ID: imprinted allele
ImprintedIDs = {}
for gene in GeneNames:
    if GeneNames[gene] in Imprinted:
        ImprintedIDs[gene] = Imprinted[GeneNames[gene]]

# count the number of external genes with imprinted status
ExtImprinted = [gene for gene in External if gene in ImprintedIDs]
# count the number of internal genes with imprinted status
IntImprinted = [gene for gene in Internal if gene in ImprintedIDs]
print('imprinted external genes', len(set(ExtImprinted)), round(len(set(ExtImprinted)) / len(External) * 100, 4))
print('imprinted internal genes', len(set(IntImprinted)), round(len(set(IntImprinted)) / len(Internal) * 100, 4))
 
# get the nested pairs in which both genes are imprinted
ImprintedPairs = []
for pair in NestedPairs:
    if pair[0] in ImprintedIDs and pair[1] in ImprintedIDs:
        ImprintedPairs.append(pair)

if len(ImprintedPairs) != 0:
    for pair in ImprintedPairs:
        print(pair[0], ImprintedIDs[pair[0]], pair[1], ImprintedIDs[pair[1]])
        
        