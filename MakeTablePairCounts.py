# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 21:47:08 2016

@author: Richard
"""

# use this script to generate a table with counts of overlapping gene pairs 

# usage python3 MakeTablePairCounts.py

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


# load dictionary of overlapping gene pairs
json_data = open('HumanOverlappingGenes.json')
OverlappingGenes = json.load(json_data)
json_data.close()

# load dictionary of nested gene pairs
json_data = open('HumanNestedGenes.json')
NestedGenes = json.load(json_data)
json_data.close()

# load dictionary of pibbyback gene pairs
json_data = open('HumanPiggyBackGenes.json')
Piggyback = json.load(json_data)
json_data.close()

# load dictionary of convergent gene pairs
json_data = open('HumanConvergentGenes.json')
Convergent = json.load(json_data)
json_data.close()


# load dictionary of divergent gene pairs
json_data = open('HumanDivergentGenes.json')
Divergent = json.load(json_data)
json_data.close()

# get GFF file
GFF = 'Homo_sapiens.GRCh38.86.gff3'
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(OverlappingGenes)
NestedPairs = GetHostNestedPairs(NestedGenes)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)

# create lists of nested gene pairs with same and opposite directions
same, opposite = [], []
for pair in NestedPairs:
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        opposite.append(pair)
    elif len(set(orientation)) == 1:
        same.append(pair)
  
# create a list of lists of gene pairs
AllPairs = [NestedPairs, same, opposite, PiggybackPairs, ConvergentPairs, DivergentPairs, OverlappingPairs]

# create set of overlapping genes
Nstgenes  = MakeFullPartialOverlapGeneSet(NestedGenes)
Ovlgenes = MakeFullPartialOverlapGeneSet(OverlappingGenes)
Cnvgenes = MakeFullPartialOverlapGeneSet(Convergent)
Divgenes = MakeFullPartialOverlapGeneSet(Divergent)
Piggygenes = MakeFullPartialOverlapGeneSet(Piggyback)
Samegenes = set()
for pair in same:
    SameGenes.add(pair[0])
    SameGenes.add(pair[1])
Oppositegenes = set()
for pair in opposite:
    Oppositegenes.add(pair[0])
    Oppositegenes.add(pair[1])

# create a table with counts of overlapping gene pairs
newfile = open('CountsOverlappingGenePairs.txt', 'w')
header = '\t'.join(['Type', 'Pairs', 'Genes', '% overlappinga', '% totalb'])
line1 = '\t'.join(['Nested', str(len(NestedPairs)), str(len(Nstgenes)), str(round((len(NestedPairs) / len(OverlappingPairs)) * 100, 2)),
                   str(round((len(Nstgenes) / len(GeneCoord)) * 100, 2))])
line2 = '\t'.join(['Nested same', str(len(same)), str(len(SameGenes)), str(round((len(same) / len(OverlappingPairs)) * 100, 2)),
                   str(round((len(SameGenes) / len(GeneCoord)) * 100, 2))])
line3 = '\t'.join(['Piggyback', str(len(PiggybackPairs)), str(len(Piggygenes)), str(round((len(PiggybackPairs) / len(OverlappingPairs)) * 100, 2)),
                   str(round((len(Piggygenes) / len(GeneCoord)) * 100, 2))])
line4 = '\t'.join(['Convergent', str(len(ConvergentPairs)), str(len(Cnvgenes)), str(round((len(ConvergentPairs) / len(OverlappingPairs)) * 100, 2)),
                   str(round((len(Cnvgenes) / len(GeneCoord)) * 100, 2))])
line5 = '\t'.join(['Divergent', str(len(DivergentPairs)), str(len(Divgenes)), str(round((len(DivergentPairs) / len(OverlappingPairs)) * 100, 2)),
                   str(round((len(Divgenes) / len(GeneCoord)) * 100, 2))])
line6 = '\t'.join(['Overlapping', str(len(OverlappingPairs)), str(len(Ovlgenes)), str(round((len(OverlappingPairs) / len(OverlappingPairs)) * 100, 2)),
                   str(round((len(Ovlgenes) / len(GeneCoord)) * 100, 2))])
line7 = '\t'.join(['Total', '_', str(len(GeneCoord)), '_', str(round((len(GeneCoord) / len(GeneCoord)) * 100, 2))])

Ndots = max(list(map(lambda x: len(x), [header, line1, line2, line3, line4, line5, line6, line7]))) 
newfile.write(str(Ndots) * '-' + '\n')
newfile.write(header + '\n')
newfile.write(str(Ndots) * '-' + '\n')
for i in [line1, line2, line3, line4, line5, line6, line7]:
    newfile.write(i + '\n')
newfile.write(str(Ndots) * '-' + '\n')
newfile.close()
      
                  
