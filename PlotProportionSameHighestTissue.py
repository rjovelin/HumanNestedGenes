# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:15:18 2017

@author: RJovelin
"""

# use this script to plot proportions of genes for which each member have the 
# same highest expressed tissues or not

# usage puthon3 PlotProportionSameHighestTissue.py

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


# make a list of json files
files = ['OverlappingGenes.json', 'NestedGenes.json']
JsonFiles = list(map(lambda x: x[0] + x[1], zip(['Human'] * len(files) , files)))

# load the dictionaries of overlapping and nested genes
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
# Map Transcript names to gene names {transcript: gene}
MapTranscriptGene = TranscriptToGene(GFF)
# get the coordinates of all exons    
ExonCoord = GeneExonCoord(GFF)
ExonCoord = CleanGeneFeatureCoord(ExonCoord, MapTranscriptGene)
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
IntronCoord = GeneIntronCoord(ExonCoord)
IntronCoord = CleanGeneFeatureCoord(IntronCoord, MapTranscriptGene)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
TranscriptCoordinates = TranscriptsCoord(GFF)
# map genes to their longest transcript {gene: longest_transcript}
GeneLongestTranscript = LongestTranscript(TranscriptCoordinates, MapGeneTranscript)
# match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
Matches = MatchHostTranscriptWithNestedTranscript(Overlap[1], MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# make a set of un-nested genes
UnNestedGenes = MakeNonOverlappingGeneSet(Overlap[1], GeneCoord)
# make a set of non-overlapping genes
NonOverlapping = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)


# make a list of non-overlapping neighboring genes
Neighbors = []
# loop over each chromo
for chromo in OrderedGenes:
    # loop over genes on chromo
    for i in range(len(OrderedGenes[chromo]) -1):
        # get gene neighbors
        gene1, gene2 = OrderedGenes[chromo][i], OrderedGenes[chromo][i+1]
        # check that both genes are non-overlapping
        if gene1 in NonOverlapping and gene2 in NonOverlapping:
            Neighbors.append([gene1, gene2])
  
# make a list of nested gene pairs
NestedPairs = GetHostNestedPairs(Overlap[1])

# make a list of host and testd transcript pairs
HostNestedTSPairs = GetHostNestedPairs(Matches)
# make list of nested gene pairs for intronless and with intron internal genes
IntronlessPairs, WithIntronPairs = [], []
for i in range(len(HostNestedTSPairs)):
    # get the external and internal transcripts of the pair
    external, internal = HostNestedTSPairs[i][0], HostNestedTSPairs[i][1]
    # check if internal transcript has intron
    if internal in IntronCoord:
        WithIntronPairs.append([MapTranscriptGene[external], MapTranscriptGene[internal]])
    else:
        IntronlessPairs.append([MapTranscriptGene[external], MapTranscriptGene[internal]])

# make lists of nested gene pairs for genes in the same or opposite direction
SamePairs, OppositePairs = [], []
for pair in NestedPairs:
    # get orientation of genes in the pair
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        OppositePairs.append(pair)
    elif len(set(orientation)) == 1:
        SamePairs.append(pair)
        

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# remove pairs with non-expressed genes
Neighbors = FilterGenePairsWithoutExpression(Neighbors, ExpressionProfile, 'strict')
NestedPairs = FilterGenePairsWithoutExpression(NestedPairs, ExpressionProfile, 'strict')
IntronlessPairs = FilterGenePairsWithoutExpression(IntronlessPairs, ExpressionProfile, 'strict')
WithIntronPairs = FilterGenePairsWithoutExpression(WithIntronPairs, ExpressionProfile, 'strict')
SamePairs = FilterGenePairsWithoutExpression(SamePairs, ExpressionProfile, 'strict')
OppositePairs = FilterGenePairsWithoutExpression(OppositePairs, ExpressionProfile, 'strict')


# 1) plot the distribution of tissue expression overlap between external and internal genes

# make a list to store the overlap
TissueOverlap = []
for pair in NestedPairs:
    # set overlap
    Overlap = 0
    # loop over the expression profile of external gene
    for i in range(len(ExpressionProfile[pair[0]])):
        # compare expression in given tissue between external and internal gene
        if ExpressionProfile[pair[0]][i] != 0 and ExpressionProfile[pair[1]][i] != 0:
            Overlap += 1
    TissueOverlap.append(Overlap)
TissueOverlap.sort()
print(TissueOverlap)

# 2) plot the proportions of gene pairs for which external and internal genes
# have highest expression in same tissue 

# make a list of list with proportions for each gene category
# [non-overlapping, all nested, same strand, opposite strand, internal with intron, intronless internal]
HighestProportions = []
AllGenes = [Neighbors, NestedPairs, SamePairs, OppositePairs, WithIntronPairs, IntronlessPairs]
for i in range(len(AllGenes)):
    HighestProportions.append([0,0])
assert len(AllGenes) == len(HighestProportions)
# loop over gene categories
for i in range(len(AllGenes)):
    # loop over gene pairs in the current gene category
    for pair in AllGenes[i]:
        # get index of tissue with maximum expression for each gene
        pos1 = ExpressionProfile[pair[0]].index(max(ExpressionProfile[pair[0]]))
        pos2 = ExpressionProfile[pair[1]].index(max(ExpressionProfile[pair[1]]))
        if pos1 == pos2:
            HighestProportions[i][0] += 1
        else:
            HighestProportions[i][1] += 1
    assert sum(HighestProportions[i]) == len(AllGenes[i])

# divide by the total number of genes in each category to get proportions
for i in range(len(HighestProportions)):
    total = sum(HighestProportions[i])
    for j in range(len(HighestProportions[i])):
        HighestProportions[i][j] = HighestProportions[i][j] / total     
    assert sum(HighestProportions[i]) == 1
# 3) plot the proportion of gene pairs for which external and internal genes
# have highest expression respectively in brain and in testis

# make a list with proportions of gene pairs for which the external gene has
# maximum expression in brain and internal gene has minimum expression in testis
RepulsiveProportions = []
for i in range(6):
    RepulsiveProportions.append([0,0])
# get indices of brain and testis tissues
infile = open('GTEX_Median_Normalized_FPKM.txt')
header = infile.readline().rstrip().split('\t')
BrainPos, TestisPos = header.index('Brain'), header.index('Testis')
infile.close()
# loop over gene categories
for i in range(len(AllGenes)):
    # loop over gene pairs in the current gene category
    for pair in AllGenes[i]:
        # get index of tissue with maximum expression for each gene
        pos1 = ExpressionProfile[pair[0]].index(max(ExpressionProfile[pair[0]]))
        pos2 = ExpressionProfile[pair[1]].index(max(ExpressionProfile[pair[1]]))
        # check if expression is maximum in brain for external gene
        if pos1 == BrainPos and pos2 == TestisPos:
                RepulsiveProportions[i][0] += 1
        else:
            RepulsiveProportions[i][1] += 1
        
# divide by the total number of genes in each category to get proportions
for i in range(len(RepulsiveProportions)):
    total = sum(RepulsiveProportions[i])    
    for j in range(len(RepulsiveProportions[i])):
        RepulsiveProportions[i][j] = RepulsiveProportions[i][j] / total     

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, YLabel, XLabel, YMax, GraphType):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # check graph type
    if GraphType == 'histo':
        # create a histogram
        ax.hist(Data, bins = range(0, len(Data), 2), linewidth = 0.7, histtype = 'bar', fill = False, facecolor = '#253494', edgecolor = 'black', alpha = 1)
    elif GraphType == 'highest':
        # Create a bar plot for proportions of gene pairs having highest expression in same tissue
        ax.bar([0, 0.2, 0.4, 0.6, 0.8, 1], Data[0], width = 0.2, label = 'same tissues', color= '#253494', linewidth = 0.7)
        # Create a bar plot for proportions of gene pairs having highest expression in different tissues
        ax.bar([0, 0.2, 0.4, 0.6, 0.8, 1], Data[1], width = 0.2, bottom = Data[0], label = 'different tissues', color= '#e31a1c', linewidth = 0.7)
    

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # add a range for the Y axis
    plt.ylim([0, YMax])
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
       
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right = 'off', left = 'on', labelbottom='off', 
                    colors = 'black', labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    return ax      


# create figure
fig = plt.figure(1, figsize = (5, 3))
# plot data
ax1 = CreateAx(3, 1, 1, fig, TissueOverlap, 'Pair counts', 'Tissue overlap', 300, 'histo')
ax2 = CreateAx(3, 1, 2, fig, [[i[0] for i in HighestProportions], [i[1] for i in HighestProportions]], 'Proportions', 'Tissue overlap', 1, 'highest')
ax3 = CreateAx(3, 1, 3, fig, [[i[0] for i in RepulsiveProportions], [i[1] for i in RepulsiveProportions]], 'Proportions', 'Tissue overlap', 1, 'highest')

# adjust padding between subplots
# pad controls the padding around the figure border
# hpad and wpad control the padding between subplots
plt.tight_layout(pad=0.2, w_pad=0, h_pad=0.5)

fig.savefig('truc.pdf', bbox_inches = 'tight')
