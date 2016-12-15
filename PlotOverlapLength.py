# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 12:07:55 2016

@author: RJovelin
"""


# use this script to plot the distribution of overlap length between overlapping gene pairs

# usage PlotOverlapLength.py 

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


# load dictionary of overlapping genes
json_data = open('HumanOverlappingGenes.json')
OverlappingGenes = json.load(json_data)
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

# make pairs of overlapping genes
OverlappingPairs = GetHostNestedPairs(OverlappingGenes)

# create a list to store the overlap length
OverlapLength = []
# loop over gene pairs
for pair in OverlappingPairs:
    # get the coordinates of the gene pairs
    coord1 = set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2]))
    coord2 = set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2]))
    OverlapLength.append(len(coord1.intersection(coord2)))

print(min(OverlapLength))
print(max(OverlapLength))
print(np.mean(OverlapLength))
print(np.median(OverlapLength))



# compute the distribution of overlap length
# compute the ratio of overlap to gene length




    


