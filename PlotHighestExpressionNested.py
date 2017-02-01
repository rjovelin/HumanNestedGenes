# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 14:34:07 2017

@author: RJovelin
"""

# use this script to plot the proportion of external and internal nested genes
# with highest expression in each tissue

# usage python3 PlotHighestExpressionNested.py 


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
Overlapping = json.load(json_data)
json_data.close()
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

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)
NestedPairs = GetHostNestedPairs(Nested)

# generate gene sets
NestedGenes  = MakeFullPartialOverlapGeneSet(Nested)
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

# make a list of all gene sets
AllGeneSets = [NonOverlappingGenes, InternalSameGenes, InternalOppositeGenes, ExternalSameGenes, ExternalOppositeGenes]

# remove genes that do not have expression profile
for i in range(len(AllGeneSets)):
    # make a list of genes to remove
    to_remove = [gene for gene in AllGeneSets[i] if gene not in ExpressionProfile]
    for gene in to_remove:
        AllGeneSets[i].remove(gene)

# make a list of tissues
infile = open('GTEX_Median_Normalized_FPKM.txt')
header = infile.readline().rstrip().split('\t')
Tissues = header[1:]
# replace spaces in tissue names
for i in range(len(Tissues)):
    if ' ' in Tissues[i]:
        Tissues[i] = Tissues[i].replace(' ', '_')

# make a list of genes with expression
Expressed = list(ExpressionProfile.keys())
# check that all genes have the same number of tissues
for gene in Expressed:
    assert len(ExpressionProfile[gene]) == len(Tissues)


# make a list of gene category names parallel to the list of gene sets
GeneCats = ['NoOv', 'CisInt', 'TransInt', 'CisExt', 'TransExt']


# create a dict to count the number of genes in each category with highest expression in each tissue
HighestExpression = {}
# inititialize dict with list of 0s
for i in GeneCats:
    HighestExpression[i] = [0] * len(Tissues)

# loop over the gene sets
for i in range(len(AllGeneSets)):
    # loop over each gene in given set
    for gene in AllGeneSets[i]:
        # check if gene has expression profile
        assert gene in ExpressionProfile
        # get the index of the maximum expression value
        pos = ExpressionProfile[gene].index(max(ExpressionProfile[gene]))
        # check that there is a single maximum expression value
        assert ExpressionProfile[gene].count(max(ExpressionProfile[gene])) == 1, '> 1 max value'
        # update counter at pos index
        HighestExpression[GeneCats[i]][pos] += 1

# divide by the total number of genes in each category to get proportions
for i in range(len(GeneCats)):
    for j in range(len(HighestExpression[GeneCats[i]])):
        HighestExpression[GeneCats[i]][j] = round(HighestExpression[GeneCats[i]][j] / len(AllGeneSets[i]), 4)



# create a dictionary with tissue as key and a list of gene proportions for each gene category as value
Proportions = {}
for i in range(len(Tissues)):
    Proportions[Tissues[i]] = []
    for j in range(len(GeneCats)):
        Proportions[Tissues[i]].append(HighestExpression[GeneCats[j]][i])


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YMax, YLabel):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['#fb9a99', '#9ecae1','#3182bd', '#a1d99b','#31a354'] 
    
    # plot nucleotide divergence
    ax.bar([0, 0.2, 0.4, 0.6, 0.8], Data, 0.2, color = colorscheme,
           edgecolor = 'black', linewidth = 0.5,
           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    if YLabel == True:
        ax.set_ylabel('Proportion', color = 'black',  size = 7, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    
    # add a range for the Y axis
    plt.ylim([0, YMax])
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)
    if YLabel == True:
        ax.spines["left"].set_visible(True)  
    elif YLabel == False:
        ax.spines['left'].set_visible(False)
    
    # edit tick parameters    
    if YLabel == True:
        plt.tick_params(axis='both', which='both', bottom='off', top='off',
                        right = 'off', left = 'on', labelbottom='off', 
                        colors = 'black', labelsize = 7, direction = 'out')  
    elif YLabel == False:
        plt.tick_params(axis='both', which='both', bottom='off', top='off',
                        right = 'off', left = 'off', labelleft='off', labelbottom='off', 
                        colors = 'black', labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    return ax      


# create figure
fig = plt.figure(1, figsize = (5, 3))

# plot data
j = 1
for i in range(len(Tissues)):
    tissue = Tissues[i]
    if j == 1 or j == 11 or j == 21:
        YLabel = True
    else:
        YLabel = False
    ax = CreateAx(10, 3, j, fig, Proportions[tissue], tissue.lower().replace('_', '\n'), 0.31, YLabel)
    j += 1

    # create legend
    if j == 10:
        No = mpatches.Patch(facecolor = '#fb9a99', edgecolor = 'black', linewidth = 0.5, label= 'NoOv')
        cisint = mpatches.Patch(facecolor = '#9ecae1', edgecolor = 'black', linewidth = 0.5, label= 'CisInt')
        transint = mpatches.Patch(facecolor = '#3182bd', edgecolor = 'black', linewidth = 0.5, label= 'TransInt')
        cisopp = mpatches.Patch(facecolor = '#a1d99b', edgecolor = 'black', linewidth = 0.5, label= 'CisExt')
        transopp = mpatches.Patch(facecolor = '#31a354', edgecolor = 'black', linewidth = 0.5, label= 'TransExt')
        ax.legend(handles = [No, cisint, transint, cisopp, transopp], bbox_to_anchor=(-10, 0.6), loc = 3, fontsize = 6, frameon = False, ncol = 5)
             
# adjust padding between subplots
# pad controls the padding around the figure border
# hpad and wpad control the padding between subplots
plt.tight_layout(pad=0.2, w_pad=0, h_pad=0.5)

fig.savefig('HighestExpressionNested.pdf', bbox_inches = 'tight')
fig.savefig('HighestExpressionNested.eps', bbox_inches = 'tight')

