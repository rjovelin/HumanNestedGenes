# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:32:23 2016

@author: RJovelin
"""

# use this function to plot the CDF of the ratios overlap length over gene length
# for the longer and short gene of the overlapping gene pairs 

# usage PlotRatioOverlapLength.py 

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

# create a list of lists of gene pairs
AllPairs = [OverlappingPairs, NestedPairs, PiggybackPairs, ConvergentPairs, DivergentPairs]
# create a parallel list of overlap length ratios for short and long genes
AllLengthShort, AllLengthLong = [], []
# loop over lists of gene pairs for each overlapping group
for i in range(len(AllPairs)):
    # initialize empty list
    LengthShort, LengthLong = [], []
    for pair in AllPairs[i]:
        assert GeneCoord[pair[0]][0] == GeneCoord[pair[1]][0]
        # get the coordinates of the gene pairs
        coord1 = set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2]))
        coord2 = set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2]))
        L = len(coord1.intersection(coord2))
        # check which gene is short and which gene is longer
        if len(coord1) >= len(coord2):
            # gene1 is longer
            LengthLong.append(L / len(coord1))
            LengthShort.append(L / len(coord2))
        elif len(coord1) < len(coord2):
            # gene2 is longer
            LengthLong.append(L / len(coord2))
            LengthShort.append(L / len(coord1))
    # store lists of overlap length
    AllLengthLong.append(LengthLong)
    AllLengthShort.append(LengthShort)
# get the list of overlap length for each group
OverlapLengthLong, NestedLengthLong = AllLengthLong[0], AllLengthLong[1]
PiggybackLengthLong, ConvergentLengthLong, DivergentLengthLong = AllLengthLong[2], AllLengthLong[3], AllLengthLong[4]

OverlapLengthShort, NestedLengthShort = AllLengthShort[0], AllLengthShort[1]
PiggybackLengthShort, ConvergentLengthShort, DivergentLengthShort = AllLengthShort[2], AllLengthShort[3], AllLengthShort[4]

# convert bp to Kbp
ToKb = lambda x: x / 1000
OverlapLengthLong = list(map(ToKb, OverlapLengthLong))
NestedLengthLong = list(map(ToKb, NestedLengthLong))
PiggybackLengthLong = list(map(ToKb, PiggybackLengthLong))
ConvergentLengthLong = list(map(ToKb, ConvergentLengthLong))
DivergentLengthLong = list(map(ToKb, DivergentLengthLong))

OverlapLengthShort = list(map(ToKb, OverlapLengthShort))
NestedLengthShort = list(map(ToKb, NestedLengthShort))
PiggybackLengthShort = list(map(ToKb, PiggybackLengthShort))
ConvergentLengthShort = list(map(ToKb, ConvergentLengthShort))
DivergentLengthShort = list(map(ToKb, DivergentLengthShort))

#def CombineHighValues(L, cutoff):
#    '''
#    (list, int) -> list
#    Take a list of overlap length and return a modified list with values higher
#    than cutoff equal to cutoff
#    '''
#    for i in range(len(L)):
#        if L[i] >= cutoff:
#            L[i] = cutoff + 1
#    return L
#
#OverlapLength = CombineHighValues(OverlapLength, 200)
#NestedLength = CombineHighValues(NestedLength, 200)
#PiggybackLength = CombineHighValues(PiggybackLength, 200)
#ConvergentLength = CombineHighValues(ConvergentLength, 200)
#DivergentLength = CombineHighValues(DivergentLength, 200)


# sort lists
OverlapLengthLong.sort()
NestedLengthLong.sort()
PiggybackLengthLong.sort()
ConvergentLengthLong.sort()
DivergentLengthLong.sort()

OverlapLengthShort.sort()
NestedLengthShort.sort()
PiggybackLengthShort.sort()
ConvergentLengthShort.sort()
DivergentLengthShort.sort()



# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot overlap length
graph1 = ax.step(NestedLengthLong, np.linspace(0, 1, len(NestedLengthLong), endpoint=False), linewidth = 1.2, color = '#d7191c', alpha = 0.7)
# plot pibbyback length
graph2 = ax.step(PiggybackLengthLong, np.linspace(0, 1, len(PiggybackLengthLong), endpoint=False), linewidth = 1.2, color = '#fdae61', alpha = 0.7)
# plot convergent length
graph3 = ax.step(ConvergentLengthLong, np.linspace(0, 1, len(ConvergentLengthLong), endpoint=False), linewidth = 1.2, color = '#abd9e9', alpha = 0.7)
# plot divergent length
graph4 = ax.step(DivergentLengthLong, np.linspace(0, 1, len(DivergentLengthLong), endpoint=False), linewidth = 1.2, color = '#2c7bb6', alpha = 0.7)


graph5 = ax.step(NestedLengthShort, np.linspace(0, 1, len(NestedLengthShort), endpoint=False), linewidth = 1.2, color = '#d7191c', alpha = 0.7)
# plot pibbyback length
graph6 = ax.step(PiggybackLengthShort, np.linspace(0, 1, len(PiggybackLengthShort), endpoint=False), linewidth = 1.2, color = '#fdae61', alpha = 0.7)
# plot convergent length
graph7 = ax.step(ConvergentLengthShort, np.linspace(0, 1, len(ConvergentLengthShort), endpoint=False), linewidth = 1.2, color = '#abd9e9', alpha = 0.7)
# plot divergent length
graph8 = ax.step(DivergentLengthShort, np.linspace(0, 1, len(DivergentLengthShort), endpoint=False), linewidth = 1.2, color = '#2c7bb6', alpha = 0.7)


# add label for the Y axis
ax.set_ylabel('Probability', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Overlap length (Kb)', size = 10, ha = 'center', fontname = 'Arial')

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
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# add lines
lns = graph1+graph2+graph3+graph4
# get labels
labs = ['Nested', 'Piggyback', 'Convergent', 'Divergent']
# plot legend
ax.legend(lns, labs, loc=3, fontsize = 8, frameon = False)

fig.savefig('truc.pdf', bbox_inches = 'tight')


