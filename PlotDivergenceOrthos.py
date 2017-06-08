# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 20:59:27 2017

@author: Richard
"""

# use this script to plot sequence and expression divergence between orthologs 

# usage python3 PlotDivergenceOrthos.py [options]
# -[chimp/mouse]: use human-chimp or human-mouse comparisons

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

Species = sys.argv[1]
assert Species in ['chimp', 'mouse']


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

# generate gene sets
OverlappingGeneSets = []
for i in range(len(Overlap)):
    OverlappingGeneSets.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

# create sets of internal and external nested gene pairs
NestedPairs = GetHostNestedPairs(Overlap[1])
InternalGenes, ExternalGenes = set(), set()
for pair in NestedPairs:
    ExternalGenes.add(pair[0])
    InternalGenes.add(pair[1])

# get orthologs
if Species == 'chimp':
    Orthos = MatchOrthologs('HumanChimpOrthologs.txt')
elif Species == 'mouse':
    Orthos = MatchOrthologs('HumanMouseOrthologs.txt')

# create a list of gene categories
GeneCats = ['Not', 'Nst', 'Int', 'Ext', 'Pbk', 'Con', 'Div']

# create lists of orthologous pairs for each gene category 
AllPairs = []
AllGenes = [NonOverlappingGenes, OverlappingGeneSets[1], InternalGenes, ExternalGenes,
            OverlappingGeneSets[2], OverlappingGeneSets[3], OverlappingGeneSets[4]] 
# loop over gene sets
for i in range(len(AllGenes)):
    # create a list of gene pairs
    orthologs = []
    # loop over genes in given gene set
    for gene in AllGenes[i]:
        # check that gene has ortholog
        if gene in Orthos:
            for homolog in Orthos[gene]:
                orthologs.append([gene, homolog])
    AllPairs.append(orthologs)


# 1) plot sequence divergence between orthologs for genes in each category

# create a dict with divergence values {human_gene: {ortholog: divergence}}
SeqDiv = {}
if Species == 'chimp':
    infile = open('HumanChimpSeqDiverg.txt')
elif Species == 'mouse':
    infile = open('HumanMouseSeqDiverg.txt')
infile.readline()
for line in infile:
    if line.startswith('ENSG'):
        line = line.rstrip().split('\t')
        if line[-1] != 'NA':
            if line[0] not in SeqDiv:
                SeqDiv[line[0]] = {}
            SeqDiv[line[0]][line[1]] = float(line[-1])
infile.close()            
    
# make list of divergence for each gene category
Divergence = []
for i in range(len(AllPairs)):
    nucldiv = []
    for pair in AllPairs[i]:
        if pair[0] in SeqDiv:
            if pair[1] in SeqDiv[pair[0]]:
                nucldiv.append(SeqDiv[pair[0]][pair[1]])
    Divergence.append(nucldiv)
    
 # create lists with means and SEM for divergence for each gene category
MeanDiverg, SEMDiverg = [], []
for i in range(len(Divergence)):
    MeanDiverg.append(np.mean(Divergence[i]))
    SEMDiverg.append(np.std(Divergence[i]) / math.sqrt(len(Divergence[i])))

# perform statistical tests between gene categories using permutation tests
# create list to store the p-values
PValDiverg = []
for i in range(1, len(Divergence)):
    # compare each gene category to non-overlapping genes
    P = PermutationResampling(Divergence[0], Divergence[i], 1000, statistic = np.mean)
    PValDiverg.append(P)
# replace P values with significance
PValDiverg = ConvertPToStars(PValDiverg)


# 2) compare expression divergence between orthologs

# get expression profile of human genes and genes of sister-species
if Species == 'chimp':
    HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
    SisterSpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
elif Species == 'mouse':
    HumanExpression = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
    SisterSpExpression = ParseExpressionFile('Mouse_Median_Normalized_FPKM.txt')
    # match expression profiles between mouse and human 
    HumanExpression = MatchHumanToMouseExpressionProfiles(HumanExpression)
# remove genes without expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
SisterSpExpression = RemoveGenesLackingExpression(SisterSpExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
SisterSpExpression = TransformRelativeExpression(SisterSpExpression)

# remove gene pairs if genes lack expression 
for i in range(len(AllPairs)):
    to_remove = [pair for pair in AllPairs[i] if pair[0] not in HumanExpression or pair[1] not in SisterSpExpression]
    for pair in to_remove:
        AllPairs[i].remove(pair)

# compute expression divergence between orthologs for each gene category
ExpressionDivergence = []
for i in range(len(AllPairs)):
    Div = ComputeExpressionDivergenceOrthologs(AllPairs[i], HumanExpression, SisterSpExpression)
    ExpressionDivergence.append(Div)
    
# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(ExpressionDivergence)):
    MeanExpDiv.append(np.mean(ExpressionDivergence[i]))
    SEMExpDiv.append(np.std(ExpressionDivergence[i]) / math.sqrt(len(ExpressionDivergence[i])))

# perform statistical tests between gene categories using permutation test
PValExpDiv = []
for i in range(1, len(ExpressionDivergence)):
    # compare each gene category to non-overlapping genes
    P = PermutationResampling(ExpressionDivergence[0], ExpressionDivergence[i], 1000, statistic = np.mean)
    PValExpDiv.append(P)
# convert p-values to star significance level
PValExpDiv = ConvertPToStars(PValExpDiv)

    
# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YLabel, YMax):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['lightgrey', '#f03b20', '#fd8d3c', '#feb24c', '#43a2ca', '#fee391', '#74c476']
    # plot divergence
    ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[0], 0.2, yerr = Data[1], color = colorscheme,
           edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    # set label size for all labels
    LabelSize = 6.5
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel(YLabel, color = 'black',  size = LabelSize, ha = 'center', **FigFont)
    # add ticks and lebels
    plt.xticks([0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9], XLabel, rotation = 0, size = LabelSize, color = 'black', ha = 'center', **FigFont)
    # add a range for the Y and X axes
    plt.ylim([0, YMax])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='on', top='off',
                    right = 'off', left = 'on', labelbottom='on',
                    colors = 'black', labelsize = LabelSize, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # add margins
    plt.margins(0.1)
    return ax


# plot sequence divergence, proportion of homologs, expression divergence in a single figure
fig = plt.figure(1, figsize = (4, 1.8))
# plot data
if Species == 'chimp':
    YMaxSeq, YMaxExp = 0.505, 0.305
elif Species == 'mouse':
    YMaxSeq, YMaxExp = 0.305, 0.505
    
ax1 = CreateAx(2, 1, 1, fig, [MeanDiverg, SEMDiverg], GeneCats, 'Nucleotide divergence (dN/dS)', YMaxSeq)
ax2 = CreateAx(2, 1, 2, fig, [MeanExpDiv, SEMExpDiv], GeneCats, 'Expression divergence', YMaxExp)

# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
if Species == 'chimp':
    yposSeq = [0.47, 0.50, 0.47, 0.47, 0.45, 0.45]
    yposExp = [0.25, 0.28, 0.25, 0.23, 0.22, 0.21]
elif Species == 'mouse':
    yposSeq = [0.23, 0.25, 0.23, 0.23, 0.22, 0.22]
    yposExp = [0.36, 0.38, 0.36, 0.32, 0.30, 0.28]    
for i in range(len(PValDiverg)):
    ax1.text(xpos[i], yposSeq[i], PValDiverg[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
for i in range(len(PValExpDiv)):
    ax2.text(xpos[i], yposExp[i], PValExpDiv[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)

# add subplot labels
ax1.text(-1, YMaxSeq + 0.025, 'A', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 7)
ax1.text(2.52, YMaxSeq + 0.025, 'B', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 7)

# make sure subplots do not overlap
plt.tight_layout()

# build figure name with option
FigName = 'DivergenceOrthos' + Species.title()
# save figure
fig.savefig(FigName + '.pdf', bbox_inches = 'tight')
fig.savefig(FigName + '.eps', bbox_inches = 'tight')

