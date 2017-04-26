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

# create lists of gene pairs
GenePairs = []
for i in range(len(Overlap)):
    GenePairs.append(GetHostNestedPairs(Overlap[i]))

# generate gene sets
GeneSets = []
for i in range(len(Overlap)):
    GeneSets.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

# create sets of internal and external nested genes 
Internal, External = set(), set()
for pair in GenePairs[1]:
    External.add(pair[0])
    Internal.add(pair[1])

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
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
# [Not, Nst, Int, Ext, Pgk, Con, Div]
AllGeneSets = [NonOverlappingGenes, GeneSets[1], Internal, External]
for i in range(2, len(GeneSets)):
    AllGeneSets.append(GeneSets[i])
  
# make a parallel list of lists of gene breadth
GeneBreadth = []
for i in range(len(AllGeneSets)):
    # loop over the gene in the given gene set and record its expression breadth
    expbreadth = [Breadth[gene] for gene in AllGeneSets[i] if gene in Breadth]
    GeneBreadth.append(expbreadth)
               
# create lists with means and SEM for each gene category
MeanBreadth, SEMBreadth = [], []
for i in range(len(GeneBreadth)):
    MeanBreadth.append(np.mean(GeneBreadth[i]))
    SEMBreadth.append(np.std(GeneBreadth[i]) / math.sqrt(len(GeneBreadth[i])))

# compare each overlapping gene category to the non-overlapping genes using Kolmogorov-Smirnof test
# create list to store the p-values
PVal = []
for i in range(1, len(GeneBreadth)):
    # compare each gene category to non-overlapping genes
    P = PermutationResampling(GeneBreadth[0], GeneBreadth[i], 1000, statistic = np.mean)    
    PVal.append(P)
# convert p-values to star significance level
PVal = ConvertPToStars(PVal)


# create figure
fig = plt.figure(1, figsize = (2, 1.5))
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)

# plot variable 
BarPos = [0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9]
Colors = ['lightgrey', '#f03b20', '#fd8d3c', '#feb24c', '#43a2ca', '#fee391', '#74c476']
ax.bar(BarPos, MeanBreadth, 0.1, yerr = SEMBreadth, color = Colors, edgecolor = 'black', linewidth = 0.7,
       error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
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
TickPos = [0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95]
Labels = ['Not', 'Nst', 'Int', 'Ext', 'Pbk', 'Con', 'Div']
plt.xticks(TickPos, Labels, rotation = 0, ha = 'center')
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   

StarPos = [0.2, 0.35, 0.5, 0.65, 0.8, 0.95]
if ExpBreadth == 'breadth':
    YPos = [29, 28, 31, 30, 31, 31]
elif ExpBreadth == 'specificity':
    YPos = [0.82, 0.9, 0.8, 0.79, 0.78, 0.77]

# add stars for significance
for i in range(len(PVal)):
    ax.text(StarPos[i], YPos[i], PVal[i], horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 6)

# add some space around x axis
plt.margins(0.05)

outputfile = 'Expression' + ExpBreadth.title()
fig.savefig(outputfile + '.pdf', bbox_inches = 'tight')
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
