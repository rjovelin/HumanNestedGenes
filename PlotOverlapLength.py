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
# create a parallel list of overlap length
AllLength = []
# loop over lists of gene pairs for each overlapping group
for i in range(len(AllPairs)):
    # initialize empty list
    Length = []
    for pair in AllPairs[i]:
        assert GeneCoord[pair[0]][0] == GeneCoord[pair[1]][0]
        # get the coordinates of the gene pairs
        coord1 = set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2]))
        coord2 = set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2]))
        L = len(coord1.intersection(coord2))
        Length.append(L)
    # store lists of overlap length
    AllLength.append(Length)
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

# check that lists are different
assert OverlapLength != NestedLength != PiggybackLength != ConvergentLength != DivergentLength

# sort lists
OverlapLength.sort()
NestedLength.sort()
PiggybackLength.sort()
ConvergentLength.sort()
DivergentLength.sort()

# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot nested length
graph1 = ax.step(NestedLength, np.linspace(0, 1, len(NestedLength), endpoint=False), linewidth = 1.2, color = '#d7191c', alpha = 0.7)
# plot pibbyback length
graph2 = ax.step(PiggybackLength, np.linspace(0, 1, len(PiggybackLength), endpoint=False), linewidth = 1.2, color = '#fdae61', alpha = 0.7)
# plot convergent length
graph3 = ax.step(ConvergentLength, np.linspace(0, 1, len(ConvergentLength), endpoint=False), linewidth = 1.2, color = '#abd9e9', alpha = 0.7)
# plot divergent length
graph4 = ax.step(DivergentLength, np.linspace(0, 1, len(DivergentLength), endpoint=False), linewidth = 1.2, color = '#2c7bb6', alpha = 0.7)

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


###################






    
## compute the distribution of overlap length
## compute the ratio of overlap to gene length
#
## create figure
#fig = plt.figure(1, figsize = (2, 2))
## add a plot to figure (N row, N column, plot N)
#ax = fig.add_subplot(1, 1, 1)
#
#NumBins = 100
#ax.hist(OverlapLength, range(0, 100000 + 5000, 5000), color = 'grey', edgecolor = 'black', linewidth = 1)
#
#
#
#
#
#truc = np.histogram(OverlapLength, range(0, max(OverlapLength)+10000, 10000))
#
#print(len(OverlapLength))
#print(sum(truc[0]))
#print(max(OverlapLength))
#
#
#short, medium, large, extra = 0, 0, 0, 0
#for i in OverlapLength:
#    if 0 <= i <= 10000:
#        short += 1
#    elif 10000 < i <= 50000:
#        medium += 1
#    elif 50000 < i <= 100000:
#        large += 1
#    elif i > 100000:
#        extra += 1
#print(short, medium, large, extra)
#
#
#
#
#
## set font for all text in figure
#FigFont = {'fontname':'Arial'}   
## set y axis label
#ax.set_ylabel('Number of gene pairs', size = 7, ha = 'center', **FigFont)
## set x axis label
#ax.set_xlabel('Overlap length', size = 7, ha = 'center', **FigFont)
## do not show lines around figure, keep bottow line  
#ax.spines["top"].set_visible(False)    
#ax.spines["right"].set_visible(False)    
#ax.spines["left"].set_visible(True)      
#ax.spines["bottom"].set_visible(True)
#    
## do not show ticks
#plt.tick_params(
#    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='on',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'on',          
#    labelbottom='on', # labels along the bottom edge are on
#    colors = 'black',
#    labelsize = 7,
#    direction = 'out') # ticks are outside the frame when bottom = 'on'  
#
#
#plt.xticks(rotation = 30)
#  
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')
