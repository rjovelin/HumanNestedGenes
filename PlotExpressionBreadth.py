# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 16:40:24 2017

@author: RJovelin
"""

# use this script to plot expression breadth for overlapping and non-overlapping genes

# usage python3 PlotExpressionBreadth.py 


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
# compute the expression breadth (number of tissues in which a gene is expressed
Breadth = ExpressionBreadth(ExpressionProfile)


# make a list of all gene sets
AllGeneSets = [InternalSameGenes, InternalOppositeGenes, ExternalSameGenes,
               ExternalOppositeGenes, PiggyBackGenes, ConvergentGenes,
               DivergentGenes, NonOverlappingGenes]
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

print(MeanBreadth)



## perform statistical tests between gene categories in all species
## create dict to store results {species: [P_host-nested, P_host-unnested, P_nested-unnested]}
#PValues = {}
## loop over inner lists in data list
#for i in range(len(AllData)):
#    # initialize dict with empty list
#    PValues[SpeciesNames[i]] = []
#    # loop over inner list, compare gene categories
#    for j in range(0, len(AllData[i]) -1):
#        for k in range(j+1, len(AllData[i])):
#            P = stats.ranksums(AllData[i][j], AllData[i][k])[1]
#            PValues[SpeciesNames[i]].append(P)
## print p values
#for sp in PValues:
#    print(sp, PValues[sp])
#
#print(HumanMeans)
#print(ChimpMeans)
#print(GorillaMeans)
#print(OrangutanMeans)
#print(MacaqueMeans)
#
#
## create a dict with significance level as stars
#Significance = {}
#for species in SpeciesNames:
#    # initialize dict with empty list
#    Significance[species] = [] 
#    # get the significance level
#    for pval in PValues[species]:
#        if pval >= 0.05:
#            Significance[species].append('')
#        elif pval < 0.05 and pval >= 0.01:
#            Significance[species].append('*')
#        elif pval < 0.01 and pval >= 0.001:
#            Significance[species].append('**')
#        elif pval < 0.001:
#            Significance[species].append('***')
#  
#
#
## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, figure, Means, SEM, XLabel, YLabel, YMax, YAxis):
#    '''
#    (int, int, int, list, figure_object, str, int, list, list)
#    Take the number of a column, and rows in the figure object and the position of
#    the ax in figure, a list of data, a title, a maximum value for the Y axis,
#    a list with species names and list of X axis tick positions and return an
#    ax instance in the figure
#    '''    
#    
#    # create subplot in figure
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#    # set colors
#    colorscheme = ['#a6cee3','#1f78b4','#b2df8a']
#    # plot nucleotide divergence
#    ax.bar([0, 0.2, 0.4], Means, 0.2, yerr = SEM, color = colorscheme,
#           edgecolor = 'black', linewidth = 1,
#           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
#
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write label for y and x axis
#    if YAxis == True:
#        ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
#    ax.set_xlabel(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
#    
#    # add a range for the Y axis
#    plt.ylim([0, YMax])
#       
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)    
#    if YAxis == True:
#        ax.spines["left"].set_visible(True)  
#    elif YAxis == False:
#        ax.spines["left"].set_visible(False)
#    
#    if YAxis == True:
#        # do not show ticks
#        plt.tick_params(
#            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#            which='both',      # both major and minor ticks are affected
#            bottom='off',      # ticks along the bottom edge are off
#            top='off',         # ticks along the top edge are off
#            right = 'off',
#            left = 'on',          
#            labelbottom='off', # labels along the bottom edge are on
#            colors = 'black',
#            labelsize = 8,
#            direction = 'out') # ticks are outside the frame when bottom = 'on'  
#    elif YAxis == False:
#        # do not show ticks
#        plt.tick_params(
#            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#            which='both',      # both major and minor ticks are affected
#            bottom='off',      # ticks along the bottom edge are off
#            top='off',         # ticks along the top edge are off
#            right = 'off',
#            left = 'off',          
#            labelbottom='off', # labels along the bottom edge are on
#            colors = 'black',
#            labelsize = 8,
#            labelleft = 'off',
#            direction = 'out') # ticks are outside the frame when bottom = 'on'      
#     
#    if YAxis == True:
#        # Set the tick labels font name
#        for label in ax.get_yticklabels():
#            label.set_fontname('Arial')   
#     
#    # create a margin around the x axis
#    plt.margins(0.1)
#    return ax      
#
#
## use this function to annotate the graph with significance levels
#def AddSignificance(ax, SignificanceLevel, XLine1, XLine2, YLine, XText, YText):
#    '''
#    (ax, str, num, num, num, num, num) -> ax
#    Take a matplotlib ax object, the significance level (as stars), the positions
#    of the bracket and star and return the ax with annotated significance level
#    '''
#    ax.annotate("", xy=(XLine1, YLine), xycoords='data', xytext=(XLine2, YLine), textcoords='data',
#                 arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
#    # add stars for significance
#    ax.text(XText, YText, SignificanceLevel, horizontalalignment='center', verticalalignment='center',
#            color = 'grey', fontname = 'Arial', size = 6)
#    return ax
#
#
## create figure
#fig = plt.figure(1, figsize = (5, 2.5))
#
## plot data
#ax1 = CreateAx(5, 1, 1, fig, HumanMeans, HumanSEM, 'Human', 'Expression specificity', 1, True)
#ax2 = CreateAx(5, 1, 2, fig, ChimpMeans, ChimpSEM, 'Chimp', 'Expression specificity', 1, False)
#ax3 = CreateAx(5, 1, 3, fig, GorillaMeans, GorillaSEM, 'Gorilla', 'Expression specificity', 1, False)
#ax4 = CreateAx(5, 1, 4, fig, OrangutanMeans, OrangutanSEM, 'Orangutan', 'Expression specificity', 1, False)
#ax5 = CreateAx(5, 1, 5, fig, MacaqueMeans, MacaqueSEM, 'Macaque', 'Expression specificity', 1, False)
#
## make lists with bracket and star positions
#XPos = [[0.1, 0.28, 0.8, 0.2, 0.85], [0.1, 0.5, 0.9, 0.3, 0.95], [0.32, 0.5, 0.8, 0.4, 0.85]]
#
## annotate figure to add significance
#for i in range(len(Significance['human'])):
#    if Significance['human'][i] != '':
#        ax1 = AddSignificance(ax1, Significance['human'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
#for i in range(len(Significance['chimp'])):
#    if Significance['chimp'][i]  != '':
#        ax2 = AddSignificance(ax2, Significance['chimp'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
#for i in range(len(Significance['gorilla'])):
#    if Significance['gorilla'][i] != '':
#        ax3 = AddSignificance(ax3, Significance['gorilla'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
#for i in range(len(Significance['orangoutan'])):
#    if Significance['orangoutan'][i] != '':
#        ax4 = AddSignificance(ax4, Significance['orangoutan'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
#for i in range(len(Significance['macaque'])):
#    if Significance['macaque'][i] != '':
#        ax5 = AddSignificance(ax5, Significance['macaque'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
# 
#
## add legend relative to ax1 using ax1 coordinates
#H = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'Hosts')
#N = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 1, label= 'Nested')
#U = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'Control')
#ax1.legend(handles = [H, N, U], loc = (0.5, 1), fontsize = 8, frameon = False, ncol = 3)
#
## make sure subplots do not overlap
#plt.tight_layout()
#
#outputfile = ''
#fig.savefig('truc.pdf', bbox_inches = 'tight')
#
