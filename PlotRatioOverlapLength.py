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
# create a parallel list of overlap length ratios for short and long genes (in %)
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
            LengthLong.append((L / len(coord1)) * 100)
            LengthShort.append((L / len(coord2)) * 100)
        elif len(coord1) < len(coord2):
            # gene2 is longer
            LengthLong.append((L / len(coord2)) * 100)
            LengthShort.append((L / len(coord1)) * 100)
    # store lists of overlap length
    AllLengthLong.append(LengthLong)
    AllLengthShort.append(LengthShort)
# get the list of overlap length for each group
OverlapLengthLong, NestedLengthLong = AllLengthLong[0], AllLengthLong[1]
PiggybackLengthLong, ConvergentLengthLong, DivergentLengthLong = AllLengthLong[2], AllLengthLong[3], AllLengthLong[4]

OverlapLengthShort, NestedLengthShort = AllLengthShort[0], AllLengthShort[1]
PiggybackLengthShort, ConvergentLengthShort, DivergentLengthShort = AllLengthShort[2], AllLengthShort[3], AllLengthShort[4]

# create figure
fig = plt.figure(1, figsize = (3.5, 3.5))

# use list to provide bin boundaries
# for integer values: 
# binBoundaries = range(0, max(data) + binwidth, binwidth)
# for float values:
# binBoundaries = np.arange(0, max(data) + binwidth, binwidth)
# plt.hist(data, bins = binBoundaries)

# create histograms to get the counts in each bin and the bin edges
#for i in range(1, len(AllLengthShort)):
#    data, binedges = np.histogram(AllLengthShort[i], bins = np.arange(0, max(AllLengthShort[i]) + 10, 10))
#    bincenters = 0.5 * (binedges[1:] + binedges[:-1])
#    proportions = [i / sum(data) for i in data]
#    ax.plot(bincenters, proportions, color = ColorScheme[i], linewidth = 1.2, linestyle = '-', marker =  Markers[i], markersize = 3.5, markeredgecolor = ColorScheme[i], markerfacecolor = ColorScheme[i], markeredgewidth = 1.2)
#
#
#for i in range(1, len(AllLengthLong)):
#    data, binedges = np.histogram(AllLengthLong[i], bins = np.arange(0, max(AllLengthLong[i]) + 10, 10))
#    bincenters = 0.5 * (binedges[1:] + binedges[:-1])
#    proportions = [i / sum(data) for i in data]
#    ax.plot(bincenters, proportions, color = ColorScheme[i], linewidth = 1.2, linestyle = '--', marker =  Markers[i], markersize = 3.5, markeredgecolor = ColorScheme[i], markerfacecolor = ColorScheme[i], markeredgewidth = 1.2)


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, Title):
    '''
    (int, int, int, figure_object, list, list, list, list, list, str, str, int)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, 2 lists of data, a list of bar positions, the list of tick
    positions and their labels, a list of colors, a label for the Y axis,
    a maximum value for the Y axis and return an ax instance in the figure
    '''    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data  
    # solid line is long, dotted line is short
    for i in range(len(Data)):
        data, binedges = np.histogram(Data[i], bins = np.arange(0, max(Data[i]) + 10, 10))
        bincenters = 0.5 * (binedges[1:] + binedges[:-1])
        proportions = [(j / sum(data)) * 100 for j in data]
        if i == 0:
            k = ax.plot(bincenters, proportions, color = 'black', linewidth = 1, linestyle = '-', marker =  'o', markersize = 3, markeredgecolor = 'black', markerfacecolor = 'black', markeredgewidth = 1)
            lns = k        
        elif i == 1:
            k = ax.plot(bincenters, proportions, color = 'grey', linewidth = 1, linestyle = '-', marker =  'o', markersize = 3, markeredgecolor = 'grey', markerfacecolor = 'white', markeredgewidth = 1)
            lns += k            
            
    # add label for the Y axis
    ax.set_ylabel('% of genes', size = 7, ha = 'center', fontname = 'Arial')
    # set x axis label
    ax.set_xlabel('Length ratio (%)', size = 7, ha = 'center', fontname = 'Arial')
    # do not show lines around figure, keep bottow line  
    # set title
    FigFont = {'fontname':'Arial'}   
    plt.title(Title, size = 7, color = 'black', ha = 'center', **FigFont )   
    # set lines around graph   
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(True)
    # add white space above x axis    
    for spine in ax.spines.values():
        spine.set_position(('outward', 5))

    # add x axis ticks
    plt.xticks([i for i in range(0, 120, 20)], [str(i) for i in range(0, 120, 20)])
    # add y axis ticks
    plt.yticks([i for i in range(0, 120, 20)], [str(i) for i in range(0, 120, 20)])
    # edit tick parameters
    ax.tick_params(
        axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        right = 'off',
        left = 'on',          
        labelbottom='on', # labels along the bottom edge are off 
        colors = 'black',
        labelsize = 7,
        direction = 'out') # ticks are outside the frame when bottom = 'on
    # use arial font    
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')

    # add legend
    # get labels
    labs = ['longer', 'shorter']
    # plot legend
    ax.legend(lns, labs, loc=1, fontsize = 6, frameon = False, numpoints = 1)
    
    return ax     


ax1 = CreateAx(2, 2, 1, fig, [NestedLengthLong, NestedLengthShort], 'Nested')
ax2 = CreateAx(2, 2, 2, fig, [PiggybackLengthLong, PiggybackLengthShort], 'Piggyback')
ax3 = CreateAx(2, 2, 3, fig, [ConvergentLengthLong, ConvergentLengthShort], 'Convergent')
ax4 = CreateAx(2, 2, 4, fig, [DivergentLengthLong, DivergentLengthShort], 'Divergent')






#for i in range(1, len(AllLengthShort)):
#    data, binedges = np.histogram(AllLengthShort[i], bins = np.arange(0, max(AllLengthShort[i]) + 10, 10))
#    bincenters = 0.5 * (binedges[1:] + binedges[:-1])
#    proportions = [j / sum(data) for j in data]
#    k = ax.plot(bincenters, proportions, color = ColorScheme[i-1], linewidth = 1.2, linestyle = '-')
#    if i == 1:
#        lns = k 
#    else:
#        lns += k
#
#for i in range(1, len(AllLengthLong)):
#    data, binedges = np.histogram(AllLengthLong[i], bins = np.arange(0, max(AllLengthLong[i]) + 10, 10))
#    bincenters = 0.5 * (binedges[1:] + binedges[:-1])
#    proportions = [j / sum(data) for j in data]
#    ax.plot(bincenters, proportions, color = ColorScheme[i-1], linewidth = 1.2, linestyle = '--', dash_capstyle = 'round')
#
#
## add label for the Y axis
#ax.set_ylabel('Number of genes', size = 8, ha = 'center', fontname = 'Arial')
## set x axis label
#ax.set_xlabel('Ratio Overlap / Gene length (%)', size = 8, ha = 'center', fontname = 'Arial')
#
## do not show lines around figure, keep bottow line  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)    
#ax.spines["left"].set_visible(True)      
#
## add x axis ticks
#plt.xticks([i for i in range(0, 110, 10)], [str(i) for i in range(0, 110, 10)])
#
#ax.tick_params(
#    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='on',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'on',          
#    labelbottom='on', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 8,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')
#
#
## get labels
#labs = ['Nested', 'Piggyback', 'Convergent', 'Divergent']
## plot legend
#ax.legend(lns, labs, loc=1, fontsize = 8, frameon = False)


# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')


