# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 22:25:45 2017

@author: Richard
"""

# use this script to compare expression divergence between same strand and opposite strand nested gene pairs

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


# load dictionary of nested gene pairs
json_data = open('HumanNestedGenes.json')
Nested = json.load(json_data)
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

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseGTEXExpression('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)


# create lists of gene pairs
NestedPairs = GetHostNestedPairs(Nested)
# create lists of nested gene pairs with same and opposite directions
same, opposite = [], []
for pair in NestedPairs:
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        opposite.append(pair)
    elif len(set(orientation)) == 1:
        same.append(pair)

# filter gene pairs lacking expression
AllPairs = [same, opposite]
for i in range(len(AllPairs)):
    AllPairs[i] = FilterGenePairsWithoutExpression(AllPairs[i], ExpressionProfile)

# compute expression divergence between pairs of genes
Divergence = []
for i in range(len(AllPairs)):
    Div = ComputeExpressionDivergenceGenePairs(AllPairs[i], ExpressionProfile)
    Divergence.append(Div)

# make a list of gene category names parallel to the list of gene pairs
GeneCats = ['same', 'opposite']

# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(Divergence)):
    MeanExpDiv.append(np.mean(Divergence[i]))
    SEMExpDiv.append(np.std(Divergence[i]) / math.sqrt(len(Divergence[i])))

## create figure
#fig = plt.figure(1, figsize = (3, 2))
#
## add a plot to figure (N row, N column, plot N)
#ax = fig.add_subplot(1, 1, 1)
## set colors
#colorscheme = ['grey', 'lightgrey']
## plot nucleotide divergence
#ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
#       edgecolor = 'black', linewidth = 0.5,
#       error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
## set font for all text in figure
#FigFont = {'fontname':'Arial'}   
## write y axis label
#ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
## add ticks and lebels
#plt.xticks([0.15, 0.45], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
## add a range for the Y and X axes
#plt.ylim([0, 1])
## add a margin
#plt.margins(0.1)
## do not show lines around figure  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(True)  
## edit tick parameters    
#plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                right = 'off', left = 'on', labelbottom='on',
#                colors = 'black', labelsize = 7, direction = 'out')  
## Set the tick labels font name
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')   
      

# perform statistical tests between gene categories
P = stats.ranksums(Divergence[0], Divergence[1])[1]
## convert p-values to star significance level
#if P >= 0.05:
#    P = ''
#elif P < 0.05 and P >= 0.01:
#    P = '*'
#elif P < 0.01 and P >= 0.001:
#    P = '**'
#elif P < 0.001:
#    P = '***'

print(MeanExpDiv)
print(P)








def permutation_resampling(case, control, num_samples, statistic):
    """Returns p-value that statistic for case is different
    from statistc for control."""

    observed_diff = abs(statistic(case) - statistic(control))
    num_case = len(case)

    combined = np.concatenate([case, control])
    diffs = []
    for i in range(num_samples):
        xs = np.random.permutation(combined)
        diff = np.mean(xs[:num_case]) - np.mean(xs[num_case:])
        diffs.append(diff)

    pval = (np.sum(diffs > observed_diff) +
            np.sum(diffs < -observed_diff))/float(num_samples)
    return pval, observed_diff, diffs

P, obs, diff = permutation_resampling(Divergence[0], Divergence[1], 10000, np.mean)


print(P)




















## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
#ypos = [0.55] * 8
#xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
#for i in range(len(Diff)):
#    ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
    
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')

