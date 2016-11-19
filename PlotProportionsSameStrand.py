# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 23:01:11 2016

@author: Richard
"""


# use this script to plot the proportion of same strand and opposite strand genes in nested pairs

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

# make lists of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
print('host-gene pairs in human', len(HumanHostNestedPairs))
ChimpHostNestedPairs = GetHostNestedPairs(ChimpHostGenes)
print('host-gene pairs in chimp', len(ChimpHostNestedPairs))
GorillaHostNestedPairs = GetHostNestedPairs(GorillaHostGenes)
print('host-gene pairs in gorilla', len(GorillaHostNestedPairs))
OrangOutanHostNestedPairs = GetHostNestedPairs(OrangOutanHostGenes)
print('host-gene pairs in orang-outan', len(OrangOutanHostNestedPairs))
MacaqueHostNestedPairs = GetHostNestedPairs(MacaqueHostGenes)
print('host-gene pairs in macaque', len(MacaqueHostNestedPairs))

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
PabGFF = 'Pongo_abelii.PPYG2.86.gff3'
MmlGFF = 'Macaca_mulatta.Mmul_8.0.1.86.gff3'    
    
    
# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF, PabGFF, MmlGFF]

# make a list of host-nested pai lists
Pairs = [HumanHostNestedPairs, ChimpHostNestedPairs, GorillaHostNestedPairs, OrangOutanHostNestedPairs, MacaqueHostNestedPairs]

# create parallel lists of proportions for same strand and opposite strand 
# for species ordered according to the GFF list
SameStrand, OppositeStrand = [], []

# loop over GFF files, get the gene coordinates for the corresponding species
for i in range(len(GFFs)):
    print(GFFs[i][:GFFs[i].index('.')])
    # get the coordinates of human genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    SpGeneChromoCoord = ChromoGenesCoord(GFFs[i])
    # map each gene to its mRNA transcripts
    SpMapGeneTranscript = GeneToTranscripts(GFFs[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    SpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(SpGeneChromoCoord, SpMapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    SpGeneCoord = FromChromoCoordToGeneCoord(SpGeneChromoCoord)
    # count the number of same strand and opposite strand pairs
    same, opposite = 0, 0
    for pair in Pairs[i]:
        orientation = set(GenePairOrientation(pair, SpGeneCoord))
        if len(orientation) == 1:
            same += 1
        else:
             assert len(orientation) == 2, 'there should be only 2 different signs'
             opposite += 1
    # add proportions to lists
    SameStrand.append(same / (same + opposite))
    OppositeStrand.append(opposite / (same + opposite))


# plot proportions for each species

# make a list of species names
SpeciesNames = ['Hsa', 'Ptr', 'Ggo', 'Pab', 'Mml']

# create figure
fig = plt.figure(1, figsize = (2.5, 2.8))

# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
  
# Create a horizontal bar plot for proportions of opposite strand pairs
ax.bar([0, 0.4, 0.8, 1.2, 1.6], OppositeStrand, width = 0.3, label = 'opposite strand', color= '#fc8d59')
# Create a horizontal bar plot for proportions of same strand pairs
ax.bar([0, 0.4, 0.8, 1.2, 1.6], SameStrand, width = 0.3, bottom = OppositeStrand, label = 'same strand', color= '#91bfdb')


# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write label for y axis
ax.set_ylabel('Proportion of gene pairs', color = 'black', size = 10, ha = 'center', **FigFont)
# write label for x axis
plt.xticks([0.15, 0.55, 0.95, 1.35, 1.75], SpeciesNames, ha = 'center', fontsize = 10, **FigFont)
        
# limit the y axis value range
plt.ylim([0, 1])   
        
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)  
# offset the spines
for spine in ax.spines.values():
    spine.set_position(('outward', 5))
    
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
## hide these grids behind plot objects
#ax.set_axisbelow(True)

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
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on'  
      
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
    
# create a margin around the x axis
plt.margins(0.05)
    
      
# add legend
S = mpatches.Patch(facecolor = '#91bfdb' , edgecolor = 'black', linewidth = 1, label= 'same')
O = mpatches.Patch(facecolor = '#fc8d59' , edgecolor = 'black', linewidth = 1, label= 'opposite')
ax.legend(handles = [S, O], loc = (-0.1, 1), fontsize = 8, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('ProportionSameOppositeStrands.pdf', bbox_inches = 'tight')
fig.savefig('ProportionSameOppositeStrands.eps', bbox_inches = 'tight')