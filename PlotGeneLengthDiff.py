# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 16:36:26 2016

@author: RJovelin
"""

# use this script to compare amino acid and gene length between overlapping and non-overlapping genes

# usage python3 PlotGeneLengthDiff.py

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
        
# generate gene sets
OverlappingGeneSets = []
for i in range(len(Overlap)):
    OverlappingGeneSets.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)


# create a list of gene sets
AllGenes = [NonOverlappingGenes]
for i in range(1, len(OverlappingGeneSets)):
    AllGenes.append(OverlappingGeneSets[i])

# create a parallel list of protein length length
ProtLength, GeneLength = [], []
# loop over lists of genes for each overlapping group
for i in range(len(AllGenes)):
    aaLength = [len(ProtSeq[gene]) for gene in AllGenes[i]]
    ProtLength.append(aaLength)
    dnaLength = []
    for gene in AllGenes[i]:
        dnaLength.append(GeneCoord[gene][2] - GeneCoord[gene][1])
    GeneLength.append(dnaLength)
    
# convert bp to Kbp
for i in range(len(GeneLength)):
    GeneLength[i] = list(map(lambda x: x/1000, GeneLength[i]))

# create lists with means and with SEM
ProtMeans, ProtSEM = [], []
for i in range(len(ProtLength)):
    ProtMeans.append(np.mean(ProtLength[i]))
    ProtSEM.append(np.std(ProtLength[i]) / math.sqrt(len(ProtLength[i])))
DnaMeans, DnaSEM = [], []
for i in range(len(GeneLength)):
    DnaMeans.append(np.mean(GeneLength[i]))
    DnaSEM.append(np.std(GeneLength[i]) / math.sqrt(len(GeneLength[i])))

# perform statistical tests between non-overlapping genes and each overlapping category
#using Kolmogorov-Smirnof test
# create list to store the P values
ProtPValues = []
for i in range(1, len(ProtLength)):
    P = PermutationResampling(ProtLength[0], ProtLength[i], 1000, statistic = np.mean)    
    ProtPValues.append(P)
DnaPValues = []
for i in range(1, len(GeneLength)):
    P = PermutationResampling(GeneLength[0], GeneLength[i], 1000, statistic = np.mean)    
    DnaPValues.append(P)

# convert p values to star significance
ProtPValues = ConvertPToStars(ProtPValues)
DnaPValues = ConvertPToStars(DnaPValues)

 
# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Means, SEM, YLabel, YMax):
    '''
    (int, int, int, figure_object, list, list, int)
    Take the number of column, and rows of the figure and the ax position, 
    2 lists of data, a maximum value for the Y axis and return an ax instance in the figure
    '''    
        
    # create subplot in figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot variable 
    BarPos = [0, 0.15, 0.3, 0.45, 0.6]
    Colors = ['lightgrey', '#f03b20', '#43a2ca', '#fee391', '#74c476']    
    ax.bar(BarPos, Means, 0.1, yerr = SEM, color = Colors, edgecolor = 'black', linewidth = 0.7,
           error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y
    ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    # add a range for the Y axis
    plt.ylim([0, YMax])
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
    Labels = ['Not', 'Nst', 'Pbk', 'Con', 'Div']
    plt.xticks(TickPos, Labels)
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   

    # add margins
    plt.margins(0.05)

    return ax      


# create figure
fig = plt.figure(1, figsize = (4, 2))

# plot protein and dna length    
ax1 = CreateAx(2, 1, 1, fig, ProtMeans, ProtSEM, 'Protein length', 700)
ax2 = CreateAx(2, 1, 2, fig, DnaMeans, DnaSEM, 'Gene length (Kb)', 120)


StarPos = [0.2, 0.35, 0.5, 0.65]
ProtYPos = [450, 550, 600, 550]
DnaYPos = [105, 70, 70, 70]

# add stars for significance
for i in range(len(ProtPValues)):
    ax1.text(StarPos[i], ProtYPos[i], ProtPValues[i], horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 6)
for i in range(len(DnaPValues)):
    ax2.text(StarPos[i], DnaYPos[i], DnaPValues[i], horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 6)

# add subplot labels
ax1.text(-0.30, 705, 'A', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 9)
ax1.text(0.85, 705, 'B', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 9)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
for extension in ['.pdf', '.eps', '.png']:
    fig.savefig('LengthDifferences' + extension, bbox_inches = 'tight')
