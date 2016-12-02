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
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('ChimpHostNestedGenes.json') as chimp_json_data:
    ChimpHostGenes = json.load(chimp_json_data)
with open('GorillaHostNestedGenes.json') as gorilla_json_data:
    GorillaHostGenes = json.load(gorilla_json_data)
# load dictionaries of host and contained genes
with open('HumanContainedGenes.json') as human_json_data:
    HumanContainedGenes = json.load(human_json_data)
with open('ChimpContainedGenes.json') as chimp_json_data:
    ChimpContainedGenes = json.load(chimp_json_data)
with open('GorillaContainedGenes.json') as gorilla_json_data:
    GorillaContainedGenes = json.load(gorilla_json_data)
# load dictionaries with overlapping genes
with open('HumanOverlappingGenes.json') as human_json_data:
    HumanOverlappingGenes = json.load(human_json_data)
with open('ChimpOverlappingGenes.json') as chimp_json_data:
    ChimpOverlappingGenes = json.load(chimp_json_data)
with open('GorillaOverlappingGenes.json') as gorilla_json_data:
    GorillaOverlappingGenes = json.load(gorilla_json_data)



# make lists of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
print('host-gene pairs in human', len(HumanHostNestedPairs))
ChimpHostNestedPairs = GetHostNestedPairs(ChimpHostGenes)
print('host-gene pairs in chimp', len(ChimpHostNestedPairs))
GorillaHostNestedPairs = GetHostNestedPairs(GorillaHostGenes)
print('host-gene pairs in gorilla', len(GorillaHostNestedPairs))

# make a list of host-contained gene pairs
HumanHostContainedPairs = GetHostNestedPairs(HumanContainedGenes)
print('host-contained pairs in human', len(HumanHostContainedPairs))
ChimpHostContainedPairs = GetHostNestedPairs(ChimpContainedGenes)
print('host-contained pairs in chimp', len(ChimpHostContainedPairs))
GorillaHostContainedPairs = GetHostNestedPairs(GorillaContainedGenes)
print('host-contained pairs in gorilla', len(GorillaHostContainedPairs))

# make a list of overlapping gene pairs
HumanOverlappingPairs = GetHostNestedPairs(HumanOverlappingGenes)
print('human overlapping gene pairs', len(HumanOverlappingPairs))
ChimpOverlappingPairs = GetHostNestedPairs(ChimpOverlappingGenes)
print('chimp overlapping gene pairs', len(ChimpOverlappingPairs))
GorillaOverlappingPairs = GetHostNestedPairs(GorillaOverlappingGenes)
print('gorilla overlapping gene pairs', len(GorillaOverlappingPairs))


# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
    
# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF]
# make a list of host-nested pai lists
NestedPairs = [HumanHostNestedPairs, ChimpHostNestedPairs, GorillaHostNestedPairs]
ContainedPairs = [HumanHostContainedPairs, ChimpHostContainedPairs, GorillaHostContainedPairs]
OverlappingPairs = [HumanOverlappingPairs, ChimpOverlappingPairs, GorillaOverlappingPairs]


# get gene coordinates in each species
# get gene coordinates on each chromo
HumanGeneChromoCoord = ChromoGenesCoord(GFFs[0])
# map each gene to its mRNA transcripts
HumanMapGeneTranscript = GeneToTranscripts(GFFs[0])
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
HumanGeneChromoCoord = FilterOutGenesWithoutValidTranscript(HumanGeneChromoCoord, HumanMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
HumanGeneCoord = FromChromoCoordToGeneCoord(HumanGeneChromoCoord)
# get gene coordinates on each chromo
ChimpGeneChromoCoord = ChromoGenesCoord(GFFs[1])
# map each gene to its mRNA transcripts
ChimpMapGeneTranscript = GeneToTranscripts(GFFs[1])
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
ChimpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(ChimpGeneChromoCoord, ChimpMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
ChimpGeneCoord = FromChromoCoordToGeneCoord(ChimpGeneChromoCoord)
# get gene coordinates on each chromo
GorillaGeneChromoCoord = ChromoGenesCoord(GFFs[2])
# map each gene to its mRNA transcripts
GorillaMapGeneTranscript = GeneToTranscripts(GFFs[2])
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GorillaGeneChromoCoord = FilterOutGenesWithoutValidTranscript(GorillaGeneChromoCoord, GorillaMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GorillaGeneCoord = FromChromoCoordToGeneCoord(GorillaGeneChromoCoord)


# make a list of gene coordinates for each species
GeneCoordinates = [HumanGeneCoord, ChimpGeneCoord, GorillaGeneCoord]

# create parallel lists of proportions for same strand and opposite strand 
# for species ordered according to the GFF list
NestedSame, NestedOpposite = [], []
ContainedSame, ContainedOpposite = [], []
OverlappingSame, OverlappingOpposite = [], []


# populate lists [human, chimp, gorilla]
for i in range(len(NestedPairs)):
    NestedSame.append(GetSameOppositeStrandProportions(NestedPairs[i], GeneCoordinates[i])[0])
    NestedOpposite.append(GetSameOppositeStrandProportions(NestedPairs[i], GeneCoordinates[i])[1])
for i in range(len(ContainedPairs)):
    ContainedSame.append(GetSameOppositeStrandProportions(ContainedPairs[i], GeneCoordinates[i])[0])
    ContainedOpposite.append(GetSameOppositeStrandProportions(ContainedPairs[i], GeneCoordinates[i])[1])
for i in range(len(OverlappingPairs)):
    OverlappingSame.append(GetSameOppositeStrandProportions(OverlappingPairs[i], GeneCoordinates[i])[0])
    OverlappingOpposite.append(GetSameOppositeStrandProportions(OverlappingPairs[i], GeneCoordinates[i])[1])


# plot proportions for each species

# create figure
fig = plt.figure(1, figsize = (3, 2.5))

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, OppositeStrand, SameStrand, XLabel, YAxis):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # Create a horizontal bar plot for proportions of opposite strand pairs
    ax.bar([0, 0.4, 0.8], OppositeStrand, width = 0.3, label = 'opposite strand', color= '#fc8d59')
    # Create a horizontal bar plot for proportions of same strand pairs
    ax.bar([0, 0.4, 0.8], SameStrand, width = 0.3, bottom = OppositeStrand, label = 'same strand', color= '#91bfdb')
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y and x axis
    if YAxis == True:
        ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = 8, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    # write label for x axis
    plt.xticks([0.15, 0.55, 0.95], ['Hsa', 'Ptr', 'Ggo'], ha = 'center', fontsize = 10, **FigFont)
    # limit the y axis value range
    plt.ylim([0, 1])   
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    if YAxis == True:
        ax.spines["left"].set_visible(True)  
    elif YAxis == False:
        ax.spines["left"].set_visible(False)
        
    if YAxis == True:
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
    elif YAxis == False:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='on', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 8,
            labelleft = 'off',
            direction = 'out') # ticks are outside the frame when bottom = 'on'      
    
    if YAxis == True:
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')   
    # create a margin around the x axis
    plt.margins(0.1)
    return ax      


ax1 = CreateAx(3, 1, 1, fig, OverlappingOpposite, OverlappingSame, 'Overlapping', True)
ax2 = CreateAx(3, 1, 2, fig, ContainedOpposite, ContainedSame, 'Contained', False)
ax3 = CreateAx(3, 1, 3, fig, NestedOpposite, NestedSame, 'Nested', False)

# add legend
S = mpatches.Patch(facecolor = '#91bfdb' , edgecolor = 'black', linewidth = 1, label= 'same')
O = mpatches.Patch(facecolor = '#fc8d59' , edgecolor = 'black', linewidth = 1, label= 'opposite')
ax1.legend(handles = [S, O], loc = (0.1, 1), fontsize = 8, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('ProportionSameOppositeStrands.pdf', bbox_inches = 'tight')
fig.savefig('ProportionSameOppositeStrands.eps', bbox_inches = 'tight')