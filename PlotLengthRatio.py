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

# Compute ratios of length of shorter gene over length of longer genes
NestedSameRatio, NestedOppositeRatio, PiggyRatio, ConvergentRatio, DivergentRatio, NeighborRatio = [], [], [], [], [], []
 
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
        
            ############## continue here





# create a list of lists of gene pairs
AllPairs = [NestedPairs, PiggybackPairs, ConvergentPairs, DivergentPairs]
# create a parallel list of length ratio
AllLength = []
# loop over lists of gene pairs for each overlapping group
for i in range(len(AllPairs)):
    # initialize empty list
    Length = []
    for pair in AllPairs[i]:
        # get the coordinates of the gene pairs
        L1 = length(set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2])))
        L2 = length(set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2])))
        if L1 <= L2:
            Length.append(L1 / L2)
        elif L1 > L2:
            Length.append(L2 / L1)
    # store lists of length ratio
    AllLength.append(Length)











# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$







# get the list of overlap length for each group
OverlapLength, NestedLength = AllLength[0], AllLength[1]
PiggybackLength, ConvergentLength, DivergentLength = AllLength[2], AllLength[3], AllLength[4]

# convert bp to Kbp
ToKb = lambda x: x / 1000
OverlapLength = list(map(ToKb, OverlapLength))
NestedLength = list(map(ToKb, NestedLength))
PiggybackLength = list(map(ToKb, PiggybackLength))
ConvergentLength = list(map(ToKb, ConvergentLength))
DivergentLength = list(map(ToKb, DivergentLength))

def CombineHighValues(L, cutoff):
    '''
    (list, int) -> list
    Take a list of overlap length and return a modified list with values higher
    than cutoff equal to cutoff
    '''
    for i in range(len(L)):
        if L[i] >= cutoff:
            L[i] = cutoff
    return L

OverlapLength = CombineHighValues(OverlapLength, 200)
NestedLength = CombineHighValues(NestedLength, 200)
PiggybackLength = CombineHighValues(PiggybackLength, 200)
ConvergentLength = CombineHighValues(ConvergentLength, 200)
DivergentLength = CombineHighValues(DivergentLength, 200)

# check that lists are different
assert OverlapLength != NestedLength != PiggybackLength != ConvergentLength != DivergentLength

# sort lists
OverlapLength = np.sort(OverlapLength)
NestedLength = np.sort(NestedLength)
PiggybackLength = np.sort(PiggybackLength)
ConvergentLength = np.sort(ConvergentLength)
DivergentLength = np.sort(DivergentLength)

# compute probabilities
POverlap = np.array(range(len(OverlapLength))) / len(OverlapLength)
PNested = np.array(range(len(NestedLength))) / len(NestedLength) 
PPiggy = np.array(range(len(PiggybackLength))) / len(PiggybackLength)
PConvergent = np.array(range(len(ConvergentLength))) / len(ConvergentLength)
PDivergent = np.array(range(len(DivergentLength))) / len(DivergentLength)

# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot nested length
Colors = ['#d7191c', '#fdae61', '#abd9e9', '#2c7bb6']

# plot nested length
graph1 = ax.step(NestedLength, PNested, linewidth = 1.2, color = '#984ea3', alpha = 0.7)
# plot pibbyback length
graph2 = ax.step(PiggybackLength, PPiggy, linewidth = 1.2, color = '#33a02c', alpha = 0.7)
# plot convergent length
graph3 = ax.step(ConvergentLength, PConvergent, linewidth = 1.2, color = '#ff7f00', alpha = 0.7)
# plot divergent length
graph4 = ax.step(DivergentLength, PDivergent, linewidth = 1.2, color = '#2c7bb6', alpha = 0.7)

# add label for the Y axis
ax.set_ylabel('Probability', size = 8, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Overlap length (Kb)', size = 8, ha = 'center', fontname = 'Arial')
# set x axis ticks
plt.xticks([0, 50, 100, 150, 200], ['0', '50', '100', '150', r'$\geq 200$'])


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

# add lines
lns = graph1+graph2+graph3+graph4
# get labels
labs = ['Nested', 'Piggyback', 'Convergent', 'Divergent']
# plot legend
ax.legend(lns, labs, loc=4, fontsize = 8, frameon = False)

fig.savefig('truc.pdf', bbox_inches = 'tight')

