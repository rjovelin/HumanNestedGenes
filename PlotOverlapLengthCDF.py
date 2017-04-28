# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:01:43 2017

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



# load dictionaries of overlapping genes
JsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',
             'HumanPiggyBackGenes.json', 'HumanConvergentGenes.json',
             'HumanDivergentGenes.json']
# make a list of dictionaries
Overlap = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    Overlap.append(overlapping)

# get GFF file
GFF = 'Homo_sapiens.GRCh38.88.gff3'

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
NonOverlapping = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

# make pairs of overlapping genes
OverlappingPairs = []
for i in range(len(Overlap)):
    pairs = GetHostNestedPairs(Overlap[i])
    OverlappingPairs.append(pairs)



# 1) plot the distribution of overlap length between overlapping gene pairs

# create a parallel list of overlap length
LengthData = []
# loop over lists of gene pairs for the 4 overlapping groups of interest
for i in range(1, len(OverlappingPairs)):
    # initialize empty list
    Length = []
    for pair in OverlappingPairs[i]:
        assert GeneCoord[pair[0]][0] == GeneCoord[pair[1]][0], 'chromosome should be the same'
        # get the coordinates of the gene pairs
        coord1 = set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2]))
        coord2 = set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2]))
        L = len(coord1.intersection(coord2))
        Length.append(L)
    # store lists of overlap length
    LengthData.append(Length)

# convert bp to Kbp
ToKb = lambda x: x / 1000
for i in range(len(LengthData)):
    LengthData[i] = list(map(ToKb, LengthData[i]))

# bin values > 200 kb into the same bin
for i in range(len(LengthData)):
    LengthData[i] = CombineHighValues(LengthData[i], 200)
    
# sort lists
for i in range(len(LengthData)):
    LengthData[i] = np.sort(LengthData[i])
# compute probabilities
ProbaLength = []
for i in range(len(LengthData)):
    ProbaLength.append(np.array(range(len(LengthData[i]))) / len(LengthData[i]))

# compute P values among groups
for i in range(0, len(LengthData) -1):
    for j in range(i+1, len(LengthData)):
        val, P =  stats.ks_2samp(LengthData[i], LengthData[j])
        print(i, j, val, P)
        

# 2) plot the CDF of the ratio of the length of the shorter gene over the length of the longer gene for overlapping genes and non-overlapping gene neighbors

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
for pair in OverlappingPairs[1]:
    if GeneCoord[pair[0]][-1] != GeneCoord[pair[1]][-1]:
        # opposite direction
        NestedOpposite.append(pair)
    elif GeneCoord[pair[0]][-1] == GeneCoord[pair[1]][-1]:
        # same direction
        NestedSame.append(pair)

# create list of gene pairs
AllPairs = [NestedSame, NestedOpposite, OverlappingPairs[2], OverlappingPairs[3], OverlappingPairs[4], NeighborsSame, NeighborsOpposite]
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

# make lists of data and probabilitties for same-strand and opposite strand overlapping genes
SameStrdData = [Ratios[0], Ratios[2], Ratios[5]]
SameStrdProba = [Proba[0], Proba[2], Proba[5]]
OppositeStrdData = [Ratios[1], Ratios[3], Ratios[4], Ratios[6]]
OppositeStrdProba = [Proba[1], Proba[3], Proba[4], Proba[6]]


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, Proba, Labels, YLabel, XLabel, TickPos, TickLabel, Colors, AddTitle = False, Title = ''):
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
    for i in range(len(Data)):
        graph = ax.step(Data[i], Proba[i], linewidth = 1.2, linestyle = '-', color = Colors[i])
        if i == 0:
            lns = graph
        else:
            lns += graph
    # set title
    FigFont = {'fontname':'Arial'}   
    if AddTitle == True:
        plt.title(Title, size = 7, color = 'black', ha = 'center', **FigFont )     
    
    # add label for the Y axis
    ax.set_ylabel(YLabel, size = 7, ha = 'center', **FigFont)
    # set x axis label
    ax.set_xlabel(XLabel, size = 7, ha = 'center', **FigFont)
    # set x axis ticks
    plt.xticks(TickPos, TickLabel)
    # do not show lines around figure, keep bottow line  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(True)      
    # edit tick parameters
    ax.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
        left = 'on', labelbottom='on', colors = 'black', labelsize = 7, direction = 'out')
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    # get labels
    labs = Labels
    # plot legend
    ax.legend(lns, labs, loc=4, fontsize = 7, frameon = False)
    return ax     



# create figure
fig = plt.figure(1, figsize = (5, 3))


# make list of labels
Label1 = ['Nst', 'Pbk', 'Con', 'Div']

# compare distributions based using Kolmogorov-Smirnov statistic on 2 samples.
Label2 = ['Nst', 'Pbk', 'Not']
for i in range(0, len(SameStrdData) -1):
    val, P =  stats.ks_2samp(SameStrdData[i], SameStrdData[-1])
    if P >= 0.05:
        Label2[i] = Label2[i] + ' NS'
    elif 0.01 <= P < 0.05:
        Label2[i] = Label2[i] + ' *'
    elif 0.001 <= P < 0.01:
        Label2[i] = Label2[i] + ' **'
    elif P < 0.001:
        Label2[i] = Label2[i] + ' ***'

Label3 = ['Nst', 'Con', 'Div', 'Not']
for i in range(0, len(OppositeStrdData) -1):
    val, P =  stats.ks_2samp(OppositeStrdData[i], OppositeStrdData[-1])
    if P >= 0.05:
        Label3[i] = Label3[i] + ' NS'
    elif 0.01 <= P < 0.05:
        Label3[i] = Label3[i] + ' *'
    elif 0.001 <= P < 0.01:
        Label3[i] = Label3[i] + ' **'
    elif P < 0.001:
        Label3[i] = Label3[i] + ' ***'

# create subplots    
ax1 = CreateAx(2, 1, 1, fig, LengthData, ProbaLength, Label1, 'Probability', 'Overlap length (Kb)', [0, 50, 100, 150, 200], ['0', '50', '100', '150', r'$\geq 200$'], ['#f03b20', '#43a2ca', '#fee391', '#74c476'])
ax2 = CreateAx(2, 2, 2, fig, SameStrdData, SameStrdProba, Label2, 'Probability', 'Length ratio (%)', [i for i in range(0, 110, 10)], [str(i) for i in range(0, 110, 10)], ['#f03b20', '#43a2ca', 'lightgrey'], True, 'Same strand')
ax3 = CreateAx(2, 2, 4, fig, OppositeStrdData, OppositeStrdProba, Label3, 'Probability', 'Length ratio (%)', [i for i in range(0, 110, 10)], [str(i) for i in range(0, 110, 10)], ['#f03b20', '#fee391', '#74c476', 'lightgrey'], True, 'Opposite strand')


# add subplot labels
ax1.text(-20, 1.1, 'A', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 8)
ax2.text(-20, 1.1, 'B', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 8)
ax3.text(-20, 1.1, 'C', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 8)
         
# make sure subplots do not overlap
plt.tight_layout()

    
fig.savefig('truc.pdf', bbox_inches = 'tight')
