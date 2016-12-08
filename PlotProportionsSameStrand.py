# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 23:01:11 2016

@author: Richard
"""


# use this script to plot the proportion of same and opposite strand genes in overlapping and nested pairs

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
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('HumanContainedGenes.json') as human_json_data:
    HumanContainedGenes = json.load(human_json_data)
with open('HumanOverlappingGenes.json') as human_json_data:
    HumanOverlappingGenes = json.load(human_json_data)

# make lists of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
# make a list of host-contained gene pairs
HumanHostContainedPairs = GetHostNestedPairs(HumanContainedGenes)
# make a list of overlapping gene pairs
HumanOverlappingPairs = GetHostNestedPairs(HumanOverlappingGenes)

# remove intronic nested gene pairs from overlapping gene pairs 
HumanOverlappingPairs = RemoveGenePairsFromHigherLevel(HumanOverlappingPairs, HumanHostNestedPairs)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
    
# get gene coordinates in each species
# get gene coordinates on each chromo
HumanGeneChromoCoord = ChromoGenesCoord(HsaGFF)
# map each gene to its mRNA transcripts
HumanMapGeneTranscript = GeneToTranscripts(HsaGFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
HumanGeneChromoCoord = FilterOutGenesWithoutValidTranscript(HumanGeneChromoCoord, HumanMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
HumanGeneCoord = FromChromoCoordToGeneCoord(HumanGeneChromoCoord)

# create parallel lists of proportions for same strand and opposite strand 
# for overlapping and for nested genes 
Same, Opposite = [], []
# populate lists
Same.append(GetSameOppositeStrandProportions(HumanOverlappingPairs, HumanGeneCoord)[0])
Opposite.append(GetSameOppositeStrandProportions(HumanOverlappingPairs, HumanGeneCoord)[1])
Same.append(GetSameOppositeStrandProportions(HumanHostNestedPairs, HumanGeneCoord)[0])
Opposite.append(GetSameOppositeStrandProportions(HumanHostNestedPairs, HumanGeneCoord)[1])

# create figure
fig = plt.figure(1, figsize = (2, 2))

# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# Create a horizontal bar plot for proportions of opposite strand pairs
ax.bar([0, 0.6], Opposite, width = 0.5, label = 'opposite strand', color= '#fc8d59')
# Create a horizontal bar plot for proportions of same strand pairs
ax.bar([0, 0.6], Same, width = 0.5, bottom = Opposite, label = 'same strand', color= '#91bfdb')
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write label for y and x axis
ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = 8, ha = 'center', **FigFont)
# write label for x axis
plt.xticks([0.25, 0.85], ['Overlap', 'Nested'], ha = 'center', fontsize = 10, **FigFont)
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
S = mpatches.Patch(facecolor = '#91bfdb' , edgecolor = 'black', linewidth = 1, label= 'same')
O = mpatches.Patch(facecolor = '#fc8d59' , edgecolor = 'black', linewidth = 1, label= 'opposite')
ax.legend(handles = [S, O], loc = (-0.3, 1.05), fontsize = 8, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('ProportionSameOppositeStrands.pdf', bbox_inches = 'tight')
fig.savefig('ProportionSameOppositeStrands.eps', bbox_inches = 'tight')