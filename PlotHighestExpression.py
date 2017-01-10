# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:14:49 2016

@author: RJovelin
"""

# use this script to plot the proportion of overlapping and non-overlapping genes
# with highest expression in each tissue

# usage python3 PlotHighestExpression.py 


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


# make a list of all gene sets
AllGeneSets = [NonOverlappingGenes, InternalSameGenes, InternalOppositeGenes, ExternalSameGenes,
               ExternalOppositeGenes, PiggyBackGenes, ConvergentGenes, DivergentGenes]

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
GeneCats = ['NoOv', 'CisInt', 'TransInt', 'CisExt', 'TransExt', 'Pbk', 'Conv', 'Div']


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
        if gene in ExpressionProfile:
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



for i in range(len(Tissues)):
    print(Tissues[i], HighestExpression['NoOv'][i], HighestExpression['CisInt'][i], HighestExpression['TransInt'][i], HighestExpression['CisExt'][i], HighestExpression['TransExt'][i], HighestExpression['Pbk'][i], HighestExpression['Conv'][i], HighestExpression['Div'][i], sep = '\t')

# create a dictionary with tissue as key and a list of gene proportions for each gene category as value
Proportions = {}
for i in range(len(Tissues)):
    Proportions[Tissues[i]] = []
    for j in range(len(GeneCats)):
        Proportions[Tissues[i]].append(HighestExpression[GeneCats[j]][i])


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Proportions, Title, YMax, YLabel):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['lightgrey'] * 8
    # plot nucleotide divergence
    ax.bar([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], Proportions, 0.1, color = colorscheme,
           edgecolor = 'black', linewidth = 1,
           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    if YLabel == True:
        ax.set_ylabel('Proportions', color = 'black',  size = 7, ha = 'center', **FigFont)
    
    #ax.set_xlabel(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
        
    plt.title(Title, size = 7, color = 'black', ha = 'center', **FigFont )     
    
    
#    # write x ticks
#    plt.xticks([0.05, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45], ['br', 'cb', 'ht', 'kd', 'lv', 'ts'], ha = 'center', fontsize = 8, **FigFont)
    
    # add a range for the Y axis
    plt.ylim([0, YMax])
       
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(True)  
        
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
         
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
     
    return ax      


# create figure
fig = plt.figure(1, figsize = (8, 5))


j = 1
for i in range(len(Tissues)):
    if i == 0 or i == 16:
        YLabel = True
    else:
        YLabel = False
    ax = CreateAx(15, 2, j, fig, Proportions[Tissues[i]], Tissues[i].lower().replace('_', ' '), 0.3, YLabel)
    j += 1


#black_line = mlines.Line2D([], [], color='black', marker='', linewidth = 1.2, label = Label1)
#grey_line = mlines.Line2D([], [], color='red', marker='', linewidth = 1.2, label = 'Overlapping')
#plt.legend(handles=[black_line, grey_line], bbox_to_anchor=(-12, 0.8), loc = 3, ncol = 2, fontsize = 10, frameon = False, borderaxespad = 0.)
#
#
## add legend relative to ax1 using ax1 coordinates
#H = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'Hosts')
#N = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 1, label= 'Nested')
#U = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'Control')
#ax1.legend(handles = [H, N, U], loc = (0, 1), fontsize = 8, frameon = False, ncol = 3)

# make sure subplots do not overlap
plt.tight_layout()

fig.savefig('truc.pdf', bbox_inches = 'tight')
#fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
