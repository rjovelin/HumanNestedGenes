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
    assert GeneCoord[pair[0]][0] == GeneCoord[pair[1]][0]
    # get the coordinates of the gene pairs
    coord1 = set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2]))
    coord2 = set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2]))
    L = len(coord1.intersection(coord2))
    OverlapLength.append(L)
    
# compute the distribution of overlap length
# compute the ratio of overlap to gene length

# create figure
fig = plt.figure(1, figsize = (2, 2))
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)

NumBins = 100
ax.hist(OverlapLength, range(0, 100000 + 5000, 5000), color = 'grey', edgecolor = 'black', linewidth = 1)





truc = np.histogram(OverlapLength, range(0, max(OverlapLength)+10000, 10000))

print(len(OverlapLength))
print(sum(truc[0]))
print(max(OverlapLength))


short, medium, large, extra = 0, 0, 0, 0
for i in OverlapLength:
    if 0 <= i <= 10000:
        short += 1
    elif 10000 < i <= 50000:
        medium += 1
    elif 50000 < i <= 100000:
        large += 1
    elif i > 100000:
        extra += 1
print(short, medium, large, extra)





# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# set y axis label
ax.set_ylabel('Number of gene pairs', size = 7, ha = 'center', **FigFont)
# set x axis label
ax.set_xlabel('Overlap length', size = 7, ha = 'center', **FigFont)
# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)      
ax.spines["bottom"].set_visible(True)
    
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
    labelsize = 7,
    direction = 'out') # ticks are outside the frame when bottom = 'on'  


plt.xticks(rotation = 30)
  
# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
