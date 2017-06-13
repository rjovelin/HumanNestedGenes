# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 21:26:13 2017

@author: Richard
"""

# use this script to plot expression divergence between host and nested genes  of same and opposite orientation 


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


# load dictionaries of host and nested genes 
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanNestedGenes.json') as human_json_data:
    Nested = json.load(human_json_data)

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

# make a a list of host, nested gene pairs
NestedPairs = GetHostNestedPairs(Nested)

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# generate gene pairs of external and internal genes with same and opposite orientation
Same, Opposite = [], []
for pair in NestedPairs:
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        Opposite.append(pair)
    elif len(set(orientation)) == 1:
        Same.append(pair)

# make a list of lists of gene pairs
OrientationPairs = [Same, Opposite]
# remove gene pairs if any gene in the pair lacls expression
for i in range(len(OrientationPairs)):
    OrientationPairs[i] = FilterGenePairsWithoutExpression(OrientationPairs[i], ExpressionProfile, 'strict')

# compute expression divergence between pairs of genes
ExpDivergOrientation = []
for i in range(len(OrientationPairs)):
    Div = ComputeExpressionDivergenceGenePairs(OrientationPairs[i], ExpressionProfile)
    ExpDivergOrientation.append(Div)
print('computed divergence')

# make a list of gene category names parallel to the list of gene pairs
GeneCatOrientation = ['Same', 'Opp']

# create lists with means and SEM for each gene category
MeanOrientation, SEMOrientation = [], []
for i in range(len(ExpDivergOrientation)):
    MeanOrientation.append(np.mean(ExpDivergOrientation[i]))
    SEMOrientation.append(np.std(ExpDivergOrientation[i]) / math.sqrt(len(ExpDivergOrientation[i])))



# create figure
fig = plt.figure(1, figsize = (2, 2))

# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# set colors
colorscheme = ['#225ea8', '#e31a1c']
# plot nucleotide divergence
ax.bar([0.05, 0.35], MeanOrientation, 0.2, yerr = SEMOrientation, color = colorscheme,
        edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.15, 0.45], GeneCatOrientation, size = 7, color = 'black', ha = 'center', **FigFont)
# add title
ax.set_xlabel('Orientation', color = 'black', size = 7, ha = 'center', **FigFont)    
# add a range for the Y and X axes
plt.ylim([0, 0.8])
plt.xlim([0, 0.6])
# edit y axis ticks
plt.yticks(np.arange(0, 1, 0.2)) 
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)
# edit tick parameters    
plt.tick_params(axis='both', which='both', bottom='on', top='off',
                right = 'off', left = 'on', labelbottom='on',
                colors = 'black', labelsize = 7, direction = 'out')  
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')    

# perform statistical tests between gene categories
PValsOrientation = [PermutationResampling(ExpDivergOrientation[0], ExpDivergOrientation[1], 1000, statistic = np.mean)]
# convert P values to stars
PValsOrientation = ConvertPToStars(PValsOrientation)[0]

# annotate figure to add significance
if PValsOrientation != '':
    ax = AddSignificanceToBars(ax, PValsOrientation, 0.15, 0.45, 0.68, 0.3, 0.72)

  
# save figure
fig.savefig('ExpDivNestedOrientation.pdf', bbox_inches = 'tight')
fig.savefig('ExpDivNestedOrientation.eps', bbox_inches = 'tight')
