# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 21:26:13 2017

@author: Richard
"""



# use this script to plot expression divergence between host and nested genes 
# separately for intron-containing and intronless nested genes 


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
GFF = 'Homo_sapiens.GRCh38.86.gff3'
 
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
HostNestedPairs = GetHostNestedPairs(Matches)


# make pairs of host-nested genes with intronless and intron-containing internal genes
WithIntrons, NoIntrons = [], []
# loop over host-nested transcript pairs
for i in range(len(HostNestedPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
    if HostNestedPairs[i][1] in IntronCoord:
        # internal gene has introns, add gene pairs to list 
            WithIntrons.append([MapTranscriptGene[HostNestedPairs[i][0]], MapTranscriptGene[HostNestedPairs[i][1]]])
    else:
        # internal gene is intronless, add gene pair to list
        NoIntrons.append([MapTranscriptGene[HostNestedPairs[i][0]], MapTranscriptGene[HostNestedPairs[i][1]]])


# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseGTEXExpression('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# filter gene pairs lacking expression
AllPairs = [WithIntrons, NoIntrons]
for i in range(len(AllPairs)):
    AllPairs[i] = FilterGenePairsWithoutExpression(AllPairs[i], ExpressionProfile)

# compute expression divergence between pairs of genes
Divergence = []
for i in range(len(AllPairs)):
    Div = ComputeExpressionDivergenceGenePairs(AllPairs[i], ExpressionProfile)
    Divergence.append(Div)

# make a list of gene category names parallel to the list of gene pairs
GeneCats = ['WithIntrons', 'Intronless']

# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(Divergence)):
    MeanExpDiv.append(np.mean(Divergence[i]))
    SEMExpDiv.append(np.std(Divergence[i]) / math.sqrt(len(Divergence[i])))



print(MeanExpDiv)


# perform statistical tests between gene categories
# create list to store the p-values
PVal = stats.ranksums(Divergence[0], Divergence[1])[1]
print(PVal)

P = PermutationResampling(Divergence[0], Divergence[1], 10000, np.mean)
print(P)


P = stats.ks_2samp(Divergence[0], Divergence[1])[1]
print(P)



## convert p-values to star significance level
#Significance = []
#for pvalue in PValues:
#    if pvalue >= 0.05:
#        Significance.append('')
#    elif pvalue < 0.05 and pvalue >= 0.01:
#        Significance.append('*')
#    elif pvalue < 0.01 and pvalue >= 0.001:
#        Significance.append('**')
#    elif pvalue < 0.001:
#        Significance.append('***')
#
#
#
#
#
#
#
## create figure
#fig = plt.figure(1, figsize = (3, 2))
#
## add a plot to figure (N row, N column, plot N)
#ax = fig.add_subplot(1, 1, 1)
## set colors
#colorscheme = ['grey','grey','grey','grey', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
## plot nucleotide divergence
#ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
#       edgecolor = 'black', linewidth = 0.5,
#       error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
## set font for all text in figure
#FigFont = {'fontname':'Arial'}   
## write y axis label
#ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
## add ticks and lebels
#plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
## add a range for the Y and X axes
#plt.ylim([0, 0.61])
#plt.xlim([0, 2.45])
## do not show lines around figure  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(True)  
## edit tick parameters    
#plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                right = 'off', left = 'on', labelbottom='on',
#                colors = 'black', labelsize = 7, direction = 'out')  
## Set the tick labels font name
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')   
#      
#
#
#
## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
#ypos = [0.55] * 8
#xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
#for i in range(len(Diff)):
#    ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#    
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')



