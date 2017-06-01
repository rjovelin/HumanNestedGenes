# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:51:25 2017

@author: RJovelin
"""

# use this script to test for enrichement of overlapping genes in regions of high recombination

# usage PlotOverlapRecomb.py

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


# load dictionaries of overlapping genes
JsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',
             'HumanPiggyBackGenes.json', 'HumanConvergentGenes.json',
             'HumanDivergentGenes.json']
# make a list of dictionaries
Overlap = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    Overlap.append(overlapping)

# get GFF file
GFF = 'Homo_sapiens.GRCh38.88.gff3'
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)

# generate gene sets
GeneSets = []
for i in range(len(Overlap)):
    GeneSets.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)


# get the coordinates of recombination hot and cold spots on each chromo
RecombSpots = GetColdHotRecomSport('CS_HRR_autosomes_final.bed')

# assign overlapping genes to hot and cold spots
OverlappingHotColdSpots = []
for i in range(len(GeneSets)):
    OverlappingHotColdSpots.append(AssignGenesToRecombSpots(GeneSets[i], RecombSpots, GeneCoord))
# assign non-overlapping genes to hot and cold spots
NonOverlappingHotColdSpots = AssignGenesToRecombSpots(NonOverlappingGenes, RecombSpots, GeneCoord)

# make lists of hot and cold spots [[overlap hot, non-overlap hot], [overlap cold, non-overlap cold]]
GetL = lambda x: len(x)
AllHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[0]), map(GetL, NonOverlappingHotColdSpots))) 
NestedHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[1]), map(GetL, NonOverlappingHotColdSpots)))
PbkHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[2]), map(GetL, NonOverlappingHotColdSpots)))
ConHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[3]), map(GetL, NonOverlappingHotColdSpots)))
DivHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[4]), map(GetL, NonOverlappingHotColdSpots)))

for L in [AllHotColdSpots, NestedHotColdSpots, PbkHotColdSpots, ConHotColdSpots, DivHotColdSpots]:
    print(stats.fisher_exact(L)[1])

