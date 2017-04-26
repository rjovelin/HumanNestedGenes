# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 18:38:17 2016

@author: Richard
"""

# use this script to plot expression divergence between gene pairs as a function
# of the distance between the 2 genes


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

# create a list of list of gene pairs
AllPairs = []
for i in range(len(Overlap)):
    AllPairs.append(GetHostNestedPairs(Overlap[i]))


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

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# generate lists of gene pairs separated by distance 
Proximal, Moderate, Intermediate, Distant = GenerateSetsGenePairsDistance(GeneCoord, OrderedGenes, ExpressionProfile)

# make a list of lists of gene pairs with expression
GenePairs = [copy.deepcopy(AllPairs[i]) for i in range(1, len(AllPairs))]
for i in range(len(GenePairs)):
    # remove pairs if any gene is lacking expression
    GenePairs[i] = FilterGenePairsWithoutExpression(GenePairs[i], ExpressionProfile, 'strict')
    
# add gene pairs to Genepairs list
for L in [Proximal, Moderate, Intermediate, Distant]:
    AllPairs.append(L)

# compute expression divergence between pairs of genes
Divergence = []
for i in range(len(GenePairs)):
    Div = ComputeExpressionDivergenceGenePairs(GenePairs[i], ExpressionProfile)
    Divergence.append(Div)

# make a list of gene category names parallel to the list of gene pairs
GeneCats = ['Nst', 'Pbk', 'Conv', 'Div', 'Prox', 'Mod', 'Int', 'Dist']

# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(Divergence)):
    MeanExpDiv.append(np.mean(Divergence[i]))
    SEMExpDiv.append(np.std(Divergence[i]) / math.sqrt(len(Divergence[i])))

# create figure
fig = plt.figure(1, figsize = (3, 2))

# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# set colors
colorscheme = ['#f03b20', '#43a2ca', '#fee391', '#74c476', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']

# plot nucleotide divergence
ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
       edgecolor = 'black', linewidth = 0.5,
       error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
# add a range for the Y and X axes
plt.ylim([0, 0.61])
plt.xlim([0, 2.45])
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
# create list to store the p-values
PValues = []
# loop over inner list, compare gene categories
for i in range(0, len(Divergence) -1):
    for j in range(i+1, len(Divergence)):
        P = stats.ranksums(Divergence[i], Divergence[j])[1]
        PValues.append(P)
# print p values
for p in PValues:
    print(p)

# convert p-values to star significance level
Significance = []
for pvalue in PValues:
    if pvalue >= 0.05:
        Significance.append('')
    elif pvalue < 0.05 and pvalue >= 0.01:
        Significance.append('*')
    elif pvalue < 0.01 and pvalue >= 0.001:
        Significance.append('**')
    elif pvalue < 0.001:
        Significance.append('***')


# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
ypos = [0.55] * 8
xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
for i in range(len(Diff)):
    ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
    
# save figure
fig.savefig('ExpressionDivergenceDistance.pdf', bbox_inches = 'tight')
fig.savefig('ExpressionDivergenceDistance.eps', bbox_inches = 'tight')
