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

##### continue here






# make a list of gene category names parallel to the list of gene pairs
PairsCats = ['WithIntrons', 'Intronless']
# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(ExpressDivergence)):
    MeanExpDiv.append(np.mean(ExpressDivergence[i]))
    SEMExpDiv.append(np.std(ExpressDivergence[i]) / math.sqrt(len(ExpressDivergence[i])))

# perform statistical tests between gene categories
PExpressDiv= PermutationResampling(ExpressDivergence[0], ExpressDivergence[1], 10000, np.mean)
# convert p-values to star significance level
if PExpressDiv >= 0.05:
    PExpressDiv = ''
elif PExpressDiv < 0.05 and PExpressDiv >= 0.01:
    PExpressDiv = '*'
elif PExpressDiv < 0.01 and PExpressDiv >= 0.001:
    PExpressDiv = '**'
elif PExpressDiv < 0.001:
    PExpressDiv = '***'



# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XticksLabel, XticksPos, BarPos, YLabel, YRange, YMax, Colors):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data
    ax.bar(BarPos, Data[0], 0.2, yerr = Data[1], color = Colors,
           edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # add ticks and lebels
    plt.xticks(XticksPos, XticksLabel, rotation = 0, size = 7, color = 'black', ha = 'center', **FigFont)
    # edit y axis ticks
    plt.yticks(YRange)   
    # add a range for the Y and X axes
    plt.ylim([0, YMax])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right = 'off', left = 'on', labelbottom='on',
                    colors = 'black', labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # add margins
    plt.margins(0.1)
    return ax



# make a figure with expression and sequence divergence

# create figure
fig = plt.figure(1, figsize = (4.5, 2))
# plot data
ax1 = CreateAx(3, 1, 1, fig, [MeanExpDiv, SEMExpDiv], ['external'], [0.25], [0, 0.3], 'Expression divergence', np.arange(0, 0.9, 0.1), 0.8, ['black', 'lightgrey'])
ax2 = CreateAx(3, 1, 2, fig, [MeanOmega, SEMOmega], ['external', 'internal'], [0.23, 0.97], [0, 0.3, 0.7, 1], 'Nucleotide divergence $\it{dN/dS}$', np.arange(0, 1.2, 0.2), 1, ['black', 'lightgrey', 'black', 'lightgrey'])
ax3 = CreateAx(3, 1, 3, fig, [MeanOrthoExpDiv, SEMOrthoExpDiv], ['external', 'internal'], [0.23, 0.97], [0, 0.3, 0.7, 1], 'Expression divergence', np.arange(0, 0.50, 0.1), 0.4, ['black', 'lightgrey', 'black', 'lightgrey'])

# annotate graphs with significance level
if PExpressDiv != '':
    ax1 = AddSignificance(ax1, PExpressDiv, 0.1, 0.4, 0.61, 0.25, 0.7)
if POmega[0] != '':
    ax2 = AddSignificance(ax2, POmega[0], 0.1, 0.4, 0.5, 0.25, 0.55)
if POmega[1] != '':
    ax2 = AddSignificance(ax2, POmega[1], 0.8, 1.1, 0.9, 0.95, 0.95)
if POrthosExprx[0] != '':
    ax3 = AddSignificance(ax3, POrthosExprx[0], 0.1, 0.4, 0.23, 0.25, 0.25)
if POrthosExprx[1] != '':
    ax3 = AddSignificance(ax3, POrthosExprx[1], 0.8, 1.1, 3.1, 0.95, 3.4)


# add legend
NoH = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'Intronless internal genes')
WiH = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'Intron-containing internal genes')
ax1.legend(handles = [WiH, NoH], loc = (0, 1.1), fontsize = 7, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure to file
fig.savefig('truc.pdf', bbox_inches = 'tight')



###########################
$$$$$$$$$$$$$$$$$$$$$$$$$$$




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
        P = PermutationResampling(Divergence[i], Divergence[j], 1000, statistic = np.mean)
        print(i, j, P)
        PValues.append(P)
# print p values
for p in PValues:
    print(p)

# convert p-values to star significance level
Significance = ConvertPToStars(PValues)

# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
ypos = [0.55] * 8
xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
for i in range(len(Diff)):
    ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
    
# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
#fig.savefig('ExpressionDivergenceDistance.pdf', bbox_inches = 'tight')
#fig.savefig('ExpressionDivergenceDistance.eps', bbox_inches = 'tight')

