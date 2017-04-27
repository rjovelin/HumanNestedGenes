# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:41:55 2016

@author: RJovelin
"""


# use this script to plot the number of introns, intron length between host,
# nested and un-nested genes and intron length of host genes for gene-containing introns
# and introns without genes

# usage python3 PlotIntronComparison.py 

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
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

# load dictionaries of host and nested genes 
with open('HumanNestedGenes.json') as human_json_data:
    HostGenes = json.load(human_json_data)

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
Matches = MatchHostTranscriptWithNestedTranscript(HostGenes, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# make a set of un-nested genes
UnNestedGenes = MakeNonOverlappingGeneSet(HostGenes, GeneCoord)

   
# get the number of introns for the longest transcript of the host, nested and un-nested genes 
HostNum, NestedNum, OthersNum = CollectIntronNumbers(UnNestedGenes, Matches, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# get intron length for host transcripts
WithGene, WithoutGene = CollectHostGeneIntronLength(Matches, TranscriptCoordinates, IntronCoord)
HostLength = WithGene + WithoutGene
# get intron length for nested transcripts
NestedLength = CollectNestedGeneIntronLength(Matches, TranscriptCoordinates, IntronCoord)
# get intron length for un-nested genes
OthersLength = CollectUnNestedGeneIntronLength(UnNestedGenes,  GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# convert bp to Kbp
HostLength = list(map(lambda x: x/1000, HostLength))
NestedLength = list(map(lambda x: x/1000, NestedLength))
OthersLength = list(map(lambda x: x/1000, OthersLength))
WithGene =  list(map(lambda x: x/1000, WithGene))
WithoutGene = list(map(lambda x: x/1000, WithoutGene))


# create a list of lists with intron numbers/length for hosts, nested and un-nested genes
IntronNumbers = [HostNum, NestedNum, OthersNum]
# create a list of lists with intron length for hosts, nested and un-nested genes
IntronLength = [HostLength, NestedLength, OthersLength]
# create a list of lists with intron length of host genes with and without nested genes
HostIntrons = [WithGene, WithoutGene]


# create lists with means and SEM
NumMeans, NumSEM = [], []
for i in range(len(IntronNumbers)):
    NumMeans.append(np.mean(IntronNumbers[i]))
    NumSEM.append(np.std(IntronNumbers[i]) / math.sqrt(len(IntronNumbers[i])))
LengthMeans, LengthSEM = [], []
for i in range(len(IntronLength)):
    LengthMeans.append(np.mean(IntronLength[i]))
    LengthSEM.append(np.std(IntronLength[i]) / math.sqrt(len(IntronLength[i])))
HostIntronMeans, HostIntronSEM = [], []
for i in range(len(HostIntrons)):
    HostIntronMeans.append(np.mean(HostIntrons[i]))
    HostIntronSEM.append(np.std(HostIntrons[i]) / math.sqrt(len(HostIntrons[i])))

# perform statistical tests between gene categories
# create dict to store results
# {number or length: [P_host-nested, P_host-unnested, P_nested-unnested], host: [P_withgene_nogene]
PValues = {}
# initialize dict with empty list
PValues['number'] = []
# loop list of intron numbers
for i in range(0, len(IntronNumbers) -1):
    for j in range(i+1, len(IntronNumbers)):
        P = PermutationResampling(IntronNumbers[i], IntronNumbers[j], 1000, statistic = np.mean)
        PValues['number'].append(P)
PValues['length'] = []
# loop list of intron length
for i in range(0, len(IntronLength) -1):
    for j in range(i+1, len(IntronLength)):
        P = PermutationResampling(IntronLength[i], IntronLength[j], 1000, statistic = np.mean)
        PValues['length'].append(P)
PValues['host'] = []
# loop list of host intron length
for i in range(0, len(HostIntrons) -1):
    for j in range(i+1, len(HostIntrons)):
        P = PermutationResampling(HostIntrons[i], HostIntrons[j], 1000, statistic = np.mean)        
        PValues['host'].append(P)

# convert p values to signififance level
for i in PValues:
    PValues[i] = ConvertPToStars(PValues[i])


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Means, SEM, BarPos, TickPos, Ticklabel, ColorScheme, YLabel, YMax):
    '''
    (int, int, int, figure_object, list, list, list, list, list, str, str, int)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, 2 lists of data, a list of bar positions, the list of tick
    positions and their labels, a list of colors, a label for the Y axis,
    a maximum value for the Y axis and return an ax instance in the figure
    '''    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot variable
    ax.bar(BarPos, Means, 0.2, yerr = SEM, color = ColorScheme,
           edgecolor = 'black', linewidth = 0.5,
           error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
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
    plt.tick_params(axis='both', which='both', bottom='on', top='off', 
                    right = 'off', left = 'on', labelbottom='on', colors = 'black',
                    labelsize = 8, direction = 'out')  
    # add ticks on the x axis
    plt.xticks(TickPos, Ticklabel)    
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # create a margin around the x axis
    plt.margins(0.1)
    return ax      


# create figure
fig = plt.figure(1, figsize = (4.5, 2.2))

# plot data for intron numner
ax1 = CreateAx(3, 1, 1, fig, NumMeans, NumSEM, [0, 0.2, 0.4], [0.1, 0.3, 0.5], ['Ext', 'Int', 'Not'], ['grey','black','white'], 'N introns / gene', 20)
ax2 = CreateAx(3, 1, 2, fig, LengthMeans, LengthSEM, [0, 0.2, 0.4], [0.1, 0.3, 0.5], ['Ext', 'Int', 'Not'], ['grey','black','white'], 'Intron length (Kbp)', 20)
ax3 = CreateAx(3, 1, 3, fig, HostIntronMeans, HostIntronSEM, [0, 0.2], [0.1, 0.3], ['With', 'None'], ['grey','black'], 'Intron length (Kbp)', 70)

# make lists with bracket and star positions
XPosNum = [[0.1, 0.28, 16, 0.2, 16.7], [0.1, 0.5, 17, 0.3, 18.2], [0.32, 0.5, 16, 0.4, 16.7]]
XPosLength = [[0.1, 0.28, 16, 0.2, 16.7], [0.1, 0.5, 17, 0.3, 18.2], [0.32, 0.5, 16, 0.4, 16.7]]
XPosHost = [[0.1, 0.28, 63, 0.2, 67]]
    
# annotate figure to add significance
for i in range(len(PValues['number'])):
    if PValues['number'][i] != '':
        ax1 = AddSignificanceToBars(ax1, PValues['number'][i], XPosNum[i][0], XPosNum[i][1], XPosNum[i][2], XPosNum[i][3], XPosNum[i][4])
for i in range(len(PValues['length'])):
    if PValues['length'][i]  != '':
        ax2 = AddSignificanceToBars(ax2, PValues['length'][i], XPosLength[i][0], XPosLength[i][1], XPosLength[i][2], XPosLength[i][3], XPosLength[i][4])
for i in range(len(PValues['host'])):
    if PValues['host'][i] != '':
        ax3 = AddSignificanceToBars(ax3, PValues['host'][i], XPosHost[i][0], XPosHost[i][1], XPosHost[i][2], XPosHost[i][3], XPosHost[i][4])


# add subplot labels
ax1.text(-0.35, 21.5, 'A', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 9)
ax1.text(0.8, 21.5, 'B', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 9)
ax1.text(2.1, 21.5, 'C', horizontalalignment='center', verticalalignment='center',
         color = 'black', fontname = 'Arial', size = 9)

# make sure subplots do not overlap
plt.tight_layout()

# one can control padding between subplots with w_pad and h_pad 
#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

## padding between subplots can also be controlled with gridspec
#gs = gridspec.GridSpec(1, 3) # N rows and columns
#gs.update(wspace=0.3, hspace=0) # set the spacing between axes. 

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

#fig.savefig('IntronDifferences.pdf', bbox_inches = 'tight')
#fig.savefig('IntronDifferences.eps', bbox_inches = 'tight')
