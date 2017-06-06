# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:39:25 2017

@author: RJovelin
"""

# use this script to plot conservation of overlapping genes


# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
import matplotlib.gridspec as gridspec
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

# make parallel lists of json files with each overlapping classes [human, mouse]
NestedFiles = [i + 'NestedGenes.json' for i in ['Human', 'Mouse']]
PiggybackFiles = [i + 'PiggyBackGenes.json' for i in ['Human', 'Mouse']]
ConvergentFiles = [i + 'ConvergentGenes.json' for i in ['Human', 'Mouse']]
DivergentFiles = [i + 'DivergentGenes.json' for i in ['Human', 'Mouse']]

# make lists of dictionaries for each type of overlapping gene
AllOverlapGenes  = []
# loop over files
AllFiles = [NestedFiles, PiggybackFiles, ConvergentFiles, DivergentFiles]
for i in range(len(AllFiles)):
    OvDicts = []
    for j in range(len(AllFiles[i])):
        # load dictionary from json file
        json_data = open(AllFiles[i][j])
        OvDicts.append(json.load(json_data))
        json_data.close()
    AllOverlapGenes.append(OvDicts)
for i in range(len(AllOverlapGenes)):
    assert len(AllOverlapGenes[i]) == 2
    

# 1) determine the proportion of gene pairs that are conserved in mouse

# for each overlapping gene class, count the proportion of gene pairs with both genes
# in the same configuration, 1 gene only or none

HumanMouseConservation = []

# get nested pairs 
NstPairs = [GetHostNestedPairs(AllOverlapGenes[0][i]) for i in range(len(AllOverlapGenes[0]))]
PbkPairs = [GetHostNestedPairs(AllOverlapGenes[1][i]) for i in range(len(AllOverlapGenes[1]))]
ConPairs = [GetHostNestedPairs(AllOverlapGenes[2][i]) for i in range(len(AllOverlapGenes[2]))]
DivPairs = [GetHostNestedPairs(AllOverlapGenes[3][i]) for i in range(len(AllOverlapGenes[3]))]

# get orthologs between human and mouse
Orthologs = MatchOrthologs('HumanMouseOrthologs.txt')
# reverse dictionary 
MouseOrthologs = {}
for gene in Orthologs:
    for ortho in Orthologs[gene]:
        if ortho not in MouseOrthologs:
            MouseOrthologs[ortho] = [gene]
        else:
            MouseOrthologs[ortho].append(gene)

# make pairs of human orthologs for each mouse gene pair
PairsOrthos = []
for L in [NstPairs, PbkPairs, ConPairs, DivPairs]:
    # create a list to store the gene pairs
    GenePairs = []
    # loop over mouse gene pairs
    for pair in L[1]:
        # count only pairs in which both genes have orthos
        if pair[0] in MouseOrthologs and pair[1] in MouseOrthologs:
            for ortho1 in MouseOrthologs[pair[0]]:
                for ortho2 in MouseOrthologs[pair[1]]:
                    # remove order
                    GenePairs.append(set([ortho1, ortho2]))
    PairsOrthos.append(GenePairs)


# make a list with pair counts for each overlapping gene type
PairCounts = []
AllPairs = [NstPairs, PbkPairs, ConPairs, DivPairs]

# loop over human overlapping gene classes
for i in range(len(AllPairs)):
    # create a list with pair counts [both gene conserved, 1 conserved, non conserved]
    counts = [0, 0, 0]
    # make a set of orthologous overlapping gene
    orthosovlp = set()
    for pair in PairsOrthos[i]:
        pair = list(pair)
        for item in pair:
            orthosovlp.add(item)
    # loop over human gene pairs
    for pair in AllPairs[i][0]:
        if set(pair) in PairsOrthos[i]:
            counts[0] += 1
        else:
            if pair[0] in orthosovlp or pair[1] in orthosovlp:
                counts[1] += 1
            elif pair[0] not in orthosovlp and pair[1] not in orthosovlp:
                counts[2] += 1
    # divide by number of gene pairs to get proportions    
    assert sum(counts) == len(AllPairs[i][0])
    for j in range(len(counts)):
        counts[j] = counts[j] / len(AllPairs[i][0])
    PairCounts.append(counts)
    assert sum(counts) == 1


fig = plt.figure(1, figsize = (1.3, 1.3))

# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
a, b, c = [i[0] for i in PairCounts], [i[1] for i in PairCounts], [i[2] for i in PairCounts]
# make a list for added values for a and b
d = [a[i] + b[i] for i in range(len(a))]

## Create a bar plot for proportions of conserved gene pairs
ax.bar([0, 0.4, 0.8, 1.2], a, width = 0.3, label = '2 conserved', color= '#9e9ac8', linewidth = 0.7)
ax.bar([0, 0.4, 0.8, 1.2], b, width = 0.3, bottom = a, label = '1 conserved', color= '#fd8d3c', linewidth = 0.7)
ax.bar([0, 0.4, 0.8, 1.2], c, width = 0.3, bottom = d, label = '0 conserved', color= '#78c679', linewidth = 0.7)


LabelSize = 7

# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write label for y and x axis
ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = LabelSize, ha = 'center', **FigFont)
# write label for x axis
plt.xticks([0.15, 0.55, 0.95, 1.35], ['Nst', 'Pbk', 'Con', 'Div'], ha = 'center', fontsize = LabelSize, **FigFont)

# limit the y axis value range
plt.ylim([0, 1])   
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)  

# do not show ticks
plt.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
                left = 'on', labelbottom='on', colors = 'black', labelsize = LabelSize, direction = 'out')  
  
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   
# create a margin around the x axis
plt.margins(0.1)

# add legend
Two = mpatches.Patch(facecolor = '#9e9ac8' , edgecolor = 'black', linewidth = 0.7, label= '2')
One = mpatches.Patch(facecolor = '#fd8d3c' , edgecolor = 'black', linewidth = 0.7, label= '1')
Zero = mpatches.Patch(facecolor = '#78c679' , edgecolor = 'black', linewidth = 0.7, label= '0')
ax.legend(handles = [Two, One, Zero], loc = (-0.3, 1.05), fontsize = LabelSize, frameon = False, ncol = 3)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
