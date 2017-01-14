# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 23:01:11 2016

@author: Richard
"""


# use this script to plot the proportion of same and opposite strand genes in overlapping and nested pairs

# usage python3 PlotProportionsSameStrand.py


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


# load dictionary of host and nested genes 
with open('HumanNestedGenes.json') as human_json_data:
    NestedGenes = json.load(human_json_data)
# load dictionary of overlapping genes
with open('HumanOverlappingGenes.json') as human_json_data:
    OverlappingGenes = json.load(human_json_data)

# make lists of host-nested gene pairs
HostNestedPairs = GetHostNestedPairs(NestedGenes)
# make a list of overlapping gene pairs
OverlappingPairs = GetHostNestedPairs(OverlappingGenes)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
    
# get gene coordinates in each species
# get gene coordinates on each chromo
GeneChromoCoord = ChromoGenesCoord(HsaGFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(HsaGFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)

# make a set of non-overlapping genes
NonOverlapping = MakeNonOverlappingGeneSet(OverlappingGenes, GeneCoord)

# create pairs of non-overlapping gene neighbors
NeighborPairs = []
# loop over each chromo
for chromo in OrderedGenes:
    # loop over genes on chromo
    for i in range(len(OrderedGenes[chromo]) -1):
        # get gene neighbors
        gene1, gene2 = OrderedGenes[chromo][i], OrderedGenes[chromo][i+1]
        # check that both genes are non-overlapping
        if gene1 in NonOverlapping and gene2 in NonOverlapping:
            NeighborPairs.append([gene1, gene2])
            
# create parallel lists of proportions for same strand and opposite strand 
# for non-overlapping, overlapping and for nested genes 
Same, Opposite = [], []
# populate lists with proportions of non-overlapping genes
Same.append(GetSameOppositeStrandProportions(NeighborPairs, GeneCoord)[0])
Opposite.append(GetSameOppositeStrandProportions(NeighborPairs, GeneCoord)[1])
# populate lists with proportions of all overlapping genes
Same.append(GetSameOppositeStrandProportions(OverlappingPairs, GeneCoord)[0])
Opposite.append(GetSameOppositeStrandProportions(OverlappingPairs, GeneCoord)[1])
# populate lists with proportions of nested genes
Same.append(GetSameOppositeStrandProportions(HostNestedPairs, GeneCoord)[0])
Opposite.append(GetSameOppositeStrandProportions(HostNestedPairs, GeneCoord)[1])


# create figure
fig = plt.figure(1, figsize = (2, 2))

# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
## Create a horizontal bar plot for proportions of opposite strand pairs
ax.bar([0, 0.6, 1.2], Opposite, width = 0.5, label = 'opposite strand', color= 'black', linewidth = 0.7)
# Create a horizontal bar plot for proportions of same strand pairs
ax.bar([0, 0.6, 1.2], Same, width = 0.5, bottom = Opposite, label = 'same strand', color= 'lightgrey', linewidth = 0.7)

# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write label for y and x axis
ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = 8, ha = 'center', **FigFont)
# write label for x axis
plt.xticks([0.25, 0.85, 1.45], ['NoOvl', 'Ovl', 'Nst'], ha = 'center', fontsize = 8, **FigFont)

# limit the y axis value range
plt.ylim([0, 1])   
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)  
 
# do not show ticks
plt.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'on',          
    labelbottom='on', # labels along the bottom edge are on
    colors = 'black',
    labelsize = 8,
    direction = 'out') # ticks are outside the frame when bottom = 'on'  
  
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   
# create a margin around the x axis
plt.margins(0.1)

# add legend
S = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 1, label= 'same')
O = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 1, label= 'opposite')
ax.legend(handles = [S, O], loc = (-0.3, 1.05), fontsize = 8, frameon = False, ncol = 2)


# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('ProportionSameOppositeStrands.pdf', bbox_inches = 'tight')
fig.savefig('ProportionSameOppositeStrands.eps', bbox_inches = 'tight')