# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:40:24 2017

@author: RJovelin
"""

# use this script to plot expression breadth for overlapping and non-overlapping genes

# usage python3 PlotExpressionBreadth.py [options]
# -[breadth/specificity]: compare the number of tissues (breadth) or tissue specificity (specificity)

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


# get the type of variable to plot
ExpBreadth = sys.argv[1]
assert ExpBreadth in ['specificity', 'breadth']

# load dictionary of overlapping gene pairs
json_data = open('HumanOverlappingGenes.json')
Overlapping = json.load(json_data)
json_data.close()
# load dictionary of nested gene pairs
json_data = open('HumanNestedGenes.json')
Nested = json.load(json_data)
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
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)
NestedPairs = GetHostNestedPairs(Nested)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)

# generate gene sets
NestedGenes  = MakeFullPartialOverlapGeneSet(Nested)
OverlappingGenes = MakeFullPartialOverlapGeneSet(Overlapping)
ConvergentGenes = MakeFullPartialOverlapGeneSet(Convergent)
DivergentGenes = MakeFullPartialOverlapGeneSet(Divergent)
PiggyBackGenes = MakeFullPartialOverlapGeneSet(Piggyback)

# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)

# create lists of nested gene pairs with same and opposite directions
same, opposite = [], []
for pair in NestedPairs:
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        opposite.append(pair)
    elif len(set(orientation)) == 1:
        same.append(pair)
# create sets of internal and external nested genes depending on orientation 
InternalSameGenes, InternalOppositeGenes, ExternalSameGenes, ExternalOppositeGenes = set(), set(), set(), set()
for pair in same:
    ExternalSameGenes.add(pair[0])
    InternalSameGenes.add(pair[1])
for pair in opposite:
    ExternalOppositeGenes.add(pair[0])
    InternalOppositeGenes.add(pair[1])


# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseGTEXExpression('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# check if expression breadth is measured by the number of tissues or tissue specificity
if ExpBreadth == 'breadth':
    # compute the expression breadth (number of tissues in which a gene is expressed)
    Breadth = ExpressionBreadth(ExpressionProfile)
elif ExpBreadth == 'specificity':
    Breadth = ExpressionSpecificity(ExpressionProfile)    


# make a list of all gene sets
AllGeneSets = [NonOverlappingGenes, InternalSameGenes, InternalOppositeGenes, ExternalSameGenes,
               ExternalOppositeGenes, PiggyBackGenes, ConvergentGenes,
               DivergentGenes]
# make a parallel list of lists of gene breadth
GeneBreadth = []
for i in range(len(AllGeneSets)):
    # loop over the gene in the given gene set and record its expression breadth
    expbreadth = [Breadth[gene] for gene in AllGeneSets[i] if gene in Breadth]
    GeneBreadth.append(expbreadth)
               

# create a function to get the mean and SEM of items in a list
def GetMeanSEM(L):
    '''
    (list) -> (list, list)
    Take a list of inner lists of numbers and return a list with mean values
    and a parallel list with SEM values for each item of the outter list
    '''
    # create lists of mean and SEM
    MeanVal, SEMVal = [], []
    # loop over the outter ist
    for i in range(len(L)):
        MeanVal.append(np.mean(L[i]))
        SEMVal.append(np.std(L[i]) / math.sqrt(len(L[i])))
    return MeanVal, SEMVal

# make a list of means and SEM for each gene category
MeanBreadth, SEMBreadth = GetMeanSEM(GeneBreadth)

# compare each overlapping gene category to the non-overlapping genes
PVal = []
# loop over gene breadth for overlapping genes only
for i in range(1, len(GeneBreadth)):
    P = stats.ranksums(GeneBreadth[i], GeneBreadth[0])[1]
    PVal.append(P)
# convert p-values to star significance level
Significance = []
for pvalue in PVal:
    if pvalue >= 0.05:
        Significance.append('')
    elif pvalue < 0.05 and pvalue >= 0.01:
        Significance.append('*')
    elif pvalue < 0.01 and pvalue >= 0.001:
        Significance.append('**')
    elif pvalue < 0.001:
        Significance.append('***')
  

# create figure
fig = plt.figure(1, figsize = (2.5, 1.5))
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)

# plot variable 
BarPos = [0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05]
Colors = ['black','lightgrey','lightgrey', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
ax.bar(BarPos, MeanBreadth, 0.1, yerr = SEMBreadth, color = Colors, edgecolor = 'black', linewidth = 1,
       error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write label for y
if ExpBreadth == 'breadth':
    ax.set_ylabel('Expression breadth', color = 'black',  size = 7, ha = 'center', **FigFont)
elif ExpBreadth == 'specificity':
    ax.set_ylabel('Expression specificity', color = 'black',  size = 7, ha = 'center', **FigFont)


# add a range for the Y axis
if ExpBreadth == 'breadth':
    plt.ylim([0, 35])
elif ExpBreadth == 'specificity':
    plt.ylim([0, 1])

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)  
# edit tick paramters
plt.tick_params(axis='both', which='both', bottom='on', top='off', right='off',
                left='on', labelbottom='on', colors='black', labelsize=7, direction='out')  
# add ticks on the x axis
TickPos = [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95, 1.10]
Labels = ['NoOv', 'CisInt', 'TransInt', 'CisExt', 'TransExt', 'Pbk', 'Conv', 'Div']
plt.xticks(TickPos, Labels, rotation = 30, ha = 'right')
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   

StarPos = [0.2, 0.35, 0.5, 0.65, 0.8, 0.95, 1.10]
if ExpBreadth == 'breadth':
    YPos = [32] * 7
elif ExpBreadth == 'specificity':
    YPos = [1] * 7

# add stars for significance
for i in range(len(Significance)):
    ax.text(StarPos[i], YPos[i], Significance[i], horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 6)

# check if specificity or breadth is computed
if ExpBreadth == 'breadth':
    outputfile = 'ExpressionBreadth'
elif ExpBreadth == 'specificity':
    outputfile = 'ExpressionSpecificity'
fig.savefig(outputfile + '.pdf', bbox_inches = 'tight')
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
