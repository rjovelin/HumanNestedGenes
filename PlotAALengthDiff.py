# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 16:36:26 2016

@author: RJovelin
"""

# use this script to compare amino acid length between overlapping and non-overlapping genes

# usage python3 CompareAAlength.py

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
# get the transcript coordinates
TSCoord = TranscriptsCoord(GFF)
# map the longest transcript to each gene
Longest = LongestTranscript(TSCoord, MapGeneTranscript)
# get the CDS of all transcripts
CDS = ConvertCDSToFasta('Homo_sapiens.GRCh38.cds.all.fa')
# get the protein sequence of the longest transcript for each gene
ProtSeq = {}
for gene in Longest:
    if Longest[gene] in CDS:
        ProtSeq[gene] = TranslateCDS(CDS[Longest[gene]])
        
# make sets of genes
OverlapGenes = MakeFullPartialOverlapGeneSet(Overlapping)
NestedGenes = MakeFullPartialOverlapGeneSet(Nested)
PiggybackGenes = MakeFullPartialOverlapGeneSet(Piggyback)
ConvergentGenes = MakeFullPartialOverlapGeneSet(Convergent)
DivergentGenes = MakeFullPartialOverlapGeneSet(Divergent)
NonOverlapGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)

# create a list of gene sets
AllGenes = [NonOverlapGenes, NestedGenes, PiggybackGenes, ConvergentGenes, DivergentGenes]
# create a parallel list of protein length length
AllLength = []
# loop over lists of genes for each overlapping group
for i in range(len(AllGenes)):
    Length = [len(ProtSeq[gene]) for gene in AllGenes[i]]
    AllLength.append(Length)

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

# create lists with means and with SEM
LengthMeans, LengthSEM = GetMeanSEM(AllLength)

# perform statistical tests between non-overlapping genes and each overlapping category
# create list to store the P values
PValues = []
for i in range(1, len(AllLength)):
    P = stats.ranksums(AllLength[0], AllLength[i])[1]
    PValues.append(P)
    
# create a list with significance level as stars
Significance = []
for i in PValues:
    if i >= 0.05:
        Significance.append('')
    elif i < 0.05 and i >= 0.01:
        Significance.append('*')
    elif i < 0.01 and i >= 0.001:
        Significance.append('**')
    elif i < 0.001:
        Significance.append('***')
  
# create figure
fig = plt.figure(1, figsize = (3, 2))
# create subplot in figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)

# plot variable 
BarPos = [0, 0.15, 0.3, 0.45, 0.6]
Colors = ['black','lightgrey','lightgrey', 'lightgrey', 'lightgrey']
ax.bar(BarPos, LengthMeans, 0.1, yerr = LengthSEM, color = Colors, edgecolor = 'black', linewidth = 1,
       error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write label for y
ax.set_ylabel('Protein length', color = 'black',  size = 8, ha = 'center', **FigFont)
# add a range for the Y axis
plt.ylim([0, 700])
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)  
# edit tick paramters
plt.tick_params(axis='both', which='both', bottom='on', top='off', right='off',
                left='on', labelbottom='on', colors='black', labelsize=8, direction='out')  

# add ticks on the x axis
TickPos = [0.05, 0.2, 0.35, 0.5, 0.65]
Labels = ['NonOvl', 'Nst', 'Pbk', 'Conv', 'Div']
plt.xticks(TickPos, Labels)
    
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   

StarPos = [0.2, 0.35, 0.5, 0.65]
YPos = [450, 550, 600, 550]
for i in range(len(Significance)):
    # add stars for significance
    ax.text(StarPos[i], YPos[i], Significance[i], horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 6)

# save figure
fig.savefig('ProteinLengthDiff.pdf', bbox_inches = 'tight')
fig.savefig('ProteinLengthDiff.eps', bbox_inches = 'tight')