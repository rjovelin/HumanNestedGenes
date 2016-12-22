# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 17:59:04 2016

@author: RJovelin
"""

# use this script to plot the CDF of the ration of the length of the shorter gene 
# over the length of the longer gene for overlapping genes and non-overlapping gene neighbors

# usage python3 PlotLengthRatio.py

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
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)

# make a set of non-overlapping genes
NonOverlapping = MakeNonOverlappingGeneSet(OverlappingGenes, GeneCoord)

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(OverlappingGenes)
NestedPairs = GetHostNestedPairs(NestedGenes)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)

# create pairs of non-overlapping gene neighbors
NeighborsSame, NeighborsOpposite = [], []
# loop over each chromo
for chromo in OrderedGenes:
    # loop over genes on chromo
    for i in range(len(OrderedGenes[chromo]) -1):
        # get gene neighbors
        gene1, gene2 = OrderedGenes[chromo][i], OrderedGenes[chromo][i+1]
        # check that both genes are non-overlapping
        if gene1 in NonOverlapping and gene2 in NonOverlapping:
            # check gene orientation
            if GeneCoord[gene1][-1] != GeneCoord[gene2][-1]:
                # genes have different orientation
                NeighborsOpposite.append([gene1, gene2])
            elif GeneCoord[gene1][-1] == GeneCoord[gene2][-1]:
                # genes have same orientation
                NeighborsSame.append([gene1, gene2])

# create pairs of nested genes with same and opposite orientation
NestedSame, NestedOpposite = [], []
# loop over nested pairs
for pair in NestedPairs:
    if GeneCoord[pair[0]][-1] != GeneCoord[pair[1]][-1]:
        # opposite direction
        NestedOpposite.append(pair)
    elif GeneCoord[pair[0]][-1] == GeneCoord[pair[1]][-1]:
        # same direction
        NestedSame.append(pair)

# create a list of lists of gene pairs
#AllPairs = [NestedPairs, PiggybackPairs, ConvergentPairs, DivergentPairs]

# create list of gene pairs
AllPairs = [NestedSame, NestedOpposite, PiggybackPairs, ConvergentPairs, DivergentPairs, NeighborsSame, NeighborsOpposite]
# create a parallel list of length ratio
Ratios = []
# Compute ratios of length (in %) of shorter gene over length of longer genes
for i in range(len(AllPairs)):
    # initialize empty list
    Length = []
    for pair in AllPairs[i]:
        # get the coordinates of the gene pairs
        L1 = len(set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2])))
        L2 = len(set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2])))
        if L1 <= L2:
            Length.append((L1 / L2) * 100)
        elif L1 > L2:
            Length.append((L2 / L1) * 100)
    # store lists of length ratio
    Ratios.append(Length)

# sort list of ratios
for i in range(len(Ratios)):
    Ratios[i] = np.sort(Ratios[i])

# compute probabilities
Proba = []
for i in range(len(Ratios)):
    P = np.array(range(len(Ratios[i]))) / len(Ratios[i])
    Proba.append(P)

#NestedSameRatio = np.sort(Ratios[0])
#NestedOppositeRatio = np.sort(Ratios[1])
#PiggyRatio = np.sort(Ratios[2])
#ConvergentRatio = np.sort(Ratios[3])
#DivergentRatio = np.sort(Ratios[4])
#NeighborSameRatio = np.sort(Ratios[5])
#NeighborOppositeRatio = np.sort(Ratios[6])


# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot nested length
Colors = ['#d7191c', '#fdae61', '#abd9e9', '#2c7bb6', 'black', 'grey']


Colors = ['#984ea3', '#33a02c', '#ff7f00', '#2c7bb6', 'black', 'grey']

Colors = ['#756bb1', '#bcbddc', '#33a02c', '#ff7f00', '#2c7bb6', 'black', 'grey']

# plot data
for i in range(len(Ratios)):
    graph = ax.step(Ratios[i], Proba[i], linewidth = 1.2, color = Colors[i], alpha = 0.7)
    if i == 0:
        lns = graph
    else:
        lns += graph


# add label for the Y axis
ax.set_ylabel('Probability', size = 8, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Length ratio (%)', size = 8, ha = 'center', fontname = 'Arial')
# set x axis ticks
plt.xticks([i for i in range(0, 110, 10)], [str(i) for i in range(0, 110, 10)])
    

# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)      

ax.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'on',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 8,
    direction = 'out') # ticks are outside the frame when bottom = 'on

for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# get labels
labs = ['Nested-same', 'Nested-opposite', 'Piggyback', 'Convergent', 'Divergent', 'Neighbors-same', 'Neighbors-oppsite']
# plot legend
ax.legend(lns, labs, loc=4, fontsize = 8, frameon = False)

fig.savefig('truc.pdf', bbox_inches = 'tight')

