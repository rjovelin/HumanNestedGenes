# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 21:26:13 2017

@author: Richard
"""

# use this script to plot expression divergence between host and nested genes 
# separately for intron-containing and intronless nested genes, and for nested genes of same and opposite orientation 


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
 
# find nested and intronic-nested genes 
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
Matches = MatchHostTranscriptWithNestedTranscript(Nested, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)

# list all host, nested transcript pairs [[host, nested]]
HostNestedTSPairs = GetHostNestedPairs(Matches)
# make a a list of host, nested gene pairs
Nestedpairs = GetHostNestedPairs(Nested)

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# generate lists of gene pairs separated by distance 
Proximal, Moderate, Intermediate, Distant = GenerateSetsGenePairsDistance(GeneCoord, OrderedGenes, ExpressionProfile)
print('made list of pairs')

# generate gene pairs of external and internal genes for intronless and intron-containing internal genes
PairsWithIntrons, PairsNoIntrons = [], []
# loop over host-nested transcript pairs
for i in range(len(HostNestedTSPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedTSPairs[i][0] in MapTranscriptGene and HostNestedTSPairs[i][0] in TranscriptCoordinates    
    assert HostNestedTSPairs[i][1] in MapTranscriptGene and HostNestedTSPairs[i][1] in TranscriptCoordinates    
    if HostNestedTSPairs[i][1] in IntronCoord:
        # internal gene has introns, add gene pairs to list 
            PairsWithIntrons.append([MapTranscriptGene[HostNestedTSPairs[i][0]], MapTranscriptGene[HostNestedTSPairs[i][1]]])
    else:
        # internal gene is intronless, add gene pair to list
        PairsNoIntrons.append([MapTranscriptGene[HostNestedTSPairs[i][0]], MapTranscriptGene[HostNestedTSPairs[i][1]]])

# make a list of lists of gene pairs
IntronPairs = [PairsWithIntrons, PairsNoIntrons]
# remove gene pairs if any gene in the pair lacks expression
for i in range(len(IntronPairs)):
    IntronPairs[i] = FilterGenePairsWithoutExpression(IntronPairs[i], ExpressionProfile, 'strict')
# add gene pairs defined by distance to list
for L in [Proximal, Moderate, Intermediate, Distant]:
    IntronPairs.append(L)

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
# add gene pairs defined by distance to list
for L in [Proximal, Moderate, Intermediate, Distant]:
    OrientationPairs.append(L)


# compute expression divergence between pairs of genes
ExpDivergIntron = []
for i in range(len(IntronPairs)):
    Div = ComputeExpressionDivergenceGenePairs(IntronPairs[i], ExpressionProfile)
    ExpDivergIntron.append(Div)
ExpDivergOrientation = []
for i in range(len(OrientationPairs)):
    Div = ComputeExpressionDivergenceGenePairs(OrientationPairs[i], ExpressionProfile)
    ExpDivergOrientation.append(Div)
print('computed divergence')

# make a list of gene category names parallel to the list of gene pairs
GeneCatIntrons = ['Introns', 'Intronless', 'Prox', 'Mod', 'Int', 'Dist']
GeneCatOrientation = ['Same', 'Opposite', 'Prox', 'Mod', 'Int', 'Dist']

# create lists with means and SEM for each gene category
MeanIntron, SEMIntron = [], []
for i in range(len(ExpDivergIntron)):
    MeanIntron.append(np.mean(ExpDivergIntron[i]))
    SEMIntron.append(np.std(ExpDivergIntron[i]) / math.sqrt(len(ExpDivergIntron[i])))
MeanOrientation, SEMOrientation = [], []
for i in range(len(ExpDivergOrientation)):
    MeanOrientation.append(np.mean(ExpDivergOrientation))
    SEMOrientation.append(np.std(ExpDivergOrientation[i]) / math.sqrt(len(ExpDivergOrientation[i])))


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, GeneCats, YRange, YMax):
    '''
    return an ax object part of figure
    '''

    # add a plot to figure (N row, N column, plot N)
    ax = fig.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['#225ea8', '#e31a1c', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
    # plot nucleotide divergence
    ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55], Data[0], 0.2, yerr = Data[1], color = colorscheme,
           edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
    # add ticks and lebels
    plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
    # add a range for the Y and X axes
    plt.ylim([0, YMax])
    plt.xlim([0, 1.8])
    # edit y axis ticks
    plt.yticks(YRange) 
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
    return ax  


# create figure
fig = plt.figure(1, figsize = (3, 2))

ax1 = CreateAx(1, 2, 1, fig, [MeanIntron, SEMIntron], GeneCatIntrons, np.arange(0, 1.1, 0.1), 1)
ax2 = CreateAx(1, 2, 2, fig, [MeanOrientation, SEMOrientation], GeneCatOrientation, np.arange(0, 1.1, 0.1), 1)

## perform statistical tests between gene categories
#PValsIntron = []
## loop over inner list, compare gene categories
#for i in range(0, len(ExpDivergIntron) -1):
#    for j in range(i+1, len(ExpDivergIntron)):
#        P = PermutationResampling(ExpDivergIntron[i], ExpDivergIntron[j], 1000, statistic = np.mean)
#        print('intron', i, j, P)
#        PValsIntron.append(P)
#PValsOrientation = []
## loop over inner list, compare gene categories
#for i in range(0, len(ExpDivergOrientation) -1):
#    for j in range(i+1, len(ExpDivergOrientation)):
#        P = PermutationResampling(ExpDivergOrientation[i], ExpDivergOrientation[j], 1000, statistic = np.mean)
#        print('orientation', i, j, P)
#        PValsOrientation.append(P)
#
## convert p-values to star significance level
#Significance = ConvertPToStars(PValues)
#
## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
#ypos = [0.55] * 6
#xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65]
#for i in range(len(Diff)):
#     ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
    
# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
#fig.savefig('ExpressionDivergenceDistance.pdf', bbox_inches = 'tight')
#fig.savefig('ExpressionDivergenceDistance.eps', bbox_inches = 'tight')


