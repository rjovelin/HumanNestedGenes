# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 08:47:41 2017

@author: Richard
"""



# use this script to plot the distribution of intron location for external genes
# containing internal genes

# usage python3 PlotIntronLocation.py 

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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
    Nested = json.load(human_json_data)

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
Matches = MatchHostTranscriptWithNestedTranscript(Nested, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# make a set of un-nested genes
UnNestedGenes = MakeNonOverlappingGeneSet(Nested, GeneCoord)
# list all host, nested transcript pairs [[host, nested]]
HostNestedPairs = GetHostNestedPairs(Matches)


# 1) plot the distribution of intron number per external gene 

# make a set of external transcripts
ExternalTS = list(Matches.keys())
# make a set of external transcripts with intronless internal genes
# make a set of external transcripts with intron-containing genes
ExternalIntronless, ExternalWithIntron = set(), set()
for i in range(len(HostNestedPairs)):
    # get the external and internal transcripts of the pair
    external, internal = HostNestedPairs[i][0], HostNestedPairs[i][1]
    # check if internal transcript has intron
    if internal in IntronCoord:
        ExternalWithIntron.add(external)
    else:
        ExternalIntronless.add(external)

# count the number of intron per external gene for all genes, external genes
# with intronless internal genes and external genes with intron-containing internal genes
TotalCount, IntronlessCount, WithIntronCount = [], [], []
for gene in ExternalTS:
    TotalCount.append(len(IntronCoord[gene]))
for gene in ExternalIntronless:
    IntronlessCount.append(len(IntronCoord[gene]))
for gene in ExternalWithIntron:
    WithIntronCount.append(len(IntronCoord[gene]))


# 2) plot the distribution of intron position for gene-containing introns of external genes
# 3) plot the distribution of intron position divided by the total number of introns

IntronlessPos, WithIntronPos = [], []
IntronlessPosNorm, WithIntronPosNorm = [], []

for i in range(len(HostNestedPairs)):
    # get the external and internal transcripts of the pair
    external, internal = HostNestedPairs[i][0], HostNestedPairs[i][1]
    assert TranscriptCoordinates[external][0] == TranscriptCoordinates[internal][0]
    # find the position of the internal gene in its host gene
    # get the start position of the internal gene
    InternalStart = TranscriptCoordinates[internal][1]
    # set boolean to check if intron is found
    FoundIntron = False    
    # check orientation of the external gene
    assert TranscriptCoordinates[external][-1] in ['-', '+']    
    if TranscriptCoordinates[external][-1] == '-':
        # 1st intron of the list is last intron of transcript
        # need to loop trhough introns in reversed order        
        for j in range(len(IntronCoord[external]) -1, -1, -1):
            # check if intron contains the internal gene
            if InternalStart in range(IntronCoord[external][j][0], IntronCoord[external][j][1]):
                # record position (1-based)
                # check if internal gene is ontronless or not
                if internal in IntronCoord:
                    WithIntronPos.append(len(IntronCoord[external]) -j)
                    WithIntronPosNorm.append((len(IntronCoord[external]) -j) / len(IntronCoord[external]))                    
                else:
                    IntronlessPos.append(len(IntronCoord[external]) -j)
                    IntronlessPosNorm.append((len(IntronCoord[external]) -j) / len(IntronCoord[external]))
            # update boolean
            FoundIntron = True
            # verify that internal gene is not in tron once intron has been found
            if FoundIntron == False:
                assert InternalStart not in range(IntronCoord[external][j][0], IntronCoord[external][j][1])
    if TranscriptCoordinates[external][-1] == '+':
        # loop over introns of the external gene    
        for j in range(len(IntronCoord[external])):
            # check if intron contains the internal gene
            if InternalStart in range(IntronCoord[external][j][0], IntronCoord[external][j][1]):
                # record position (1-based)
                # check if internal gene is ontronless or not
                if internal in IntronCoord:
                    WithIntronPos.append(j+1)
                    WithIntronPosNorm.append((j+1) / len(IntronCoord[external]))
                else:
                    IntronlessPos.append(j+1)
                    IntronlessPosNorm.append((j+1) / len(IntronCoord[external]))                    
                # update boolean
                FoundIntron = True
            # verify that internal gene is not in tron once intron has been found
            if FoundIntron == False:
                assert InternalStart not in range(IntronCoord[external][j][0], IntronCoord[external][j][1])
  

# sort lists
WithIntronPosNorm = np.sort(WithIntronPosNorm)
IntronlessPosNorm = np.sort(IntronlessPosNorm)
# compute probabilities
PWithPosNorm = np.array(range(len(WithIntronPosNorm))) / len(WithIntronPosNorm)
PNonePosNorm = np.array(range(len(IntronlessPosNorm))) / len(IntronlessPosNorm) 

# compare distributions of intron numbers, intron position
PVals = []
Data = [IntronlessCount, WithIntronCount, IntronlessPos, WithIntronPos, IntronlessPosNorm, WithIntronPosNorm]
for i in range(0, len(Data), 2):
    val, P =  stats.ks_2samp(Data[i], Data[i+1])
    PVals.append(P)

print(PVals)

# get significance levels
for i in range(len(PVals)):
    if PVals[i] >= 0.05:
        PVals[i] = ''
    elif PVals[i] < 0.05 and PVals[i] >= 0.01:
        PVals[i] = '$\it{P}$ < 0.05'
    elif PVals[i] < 0.01 and PVals[i] >= 0.001:
        PVals[i] = '$\it{P}$ < 0.01'
    elif PVals[i] < 0.001:
        PVals[i] = '$\it{P}$ < 0.001'

print(PVals)  

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, YLabel, XLabel, YRange, YMax, XRange, Colors, GraphType):
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
    # plot data    
    if GraphType == 'histo':
        ax.hist(Data[0], bins = np.arange(0, max([max(Data[0]), max(Data[1])]) + 1, 1), linewidth = 0.7, histtype='step', fill = True, facecolor = Colors[0], edgecolor = Colors[0], alpha = 0.5, stacked = False)    
        ax.hist(Data[1], bins = np.arange(0, max([max(Data[0]), max(Data[1])]) + 1, 1), linewidth = 0.7, histtype='step', fill = True, facecolor = Colors[1], edgecolor = Colors[1], alpha = 0.5, stacked = False)
    elif GraphType == 'cdf':
        ax.step(Data[0], Data[1], linewidth = 1, linestyle = '-', color = Colors[0], alpha = 0.5)
        ax.step(Data[2], Data[3], linewidth = 1, linestyle = '-', color = Colors[1], alpha = 0.5)
        
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label axis
    ax.set_ylabel(YLabel, color = 'black',  size = 6.5, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 6.5, ha = 'center', **FigFont)
    # set a limit to y axis
    plt.ylim([0, YMax])
    # add a range to the axis
    plt.xticks(XRange, size = 6.5, color = 'black', ha = 'center', **FigFont)
    plt.yticks(YRange, size = 6.5, color = 'black', ha = 'right', **FigFont)        
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(True)  
    # edit tick paramters
    plt.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
                    left = 'on', labelbottom='on', colors = 'black', labelsize = 7, direction = 'out')  
    # add ticks on the x axis
    #plt.xticks(TickPos, Ticklabel)    
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # create a margin around the x axis
    plt.margins(0.05)
    return ax      

# create figure
fig = plt.figure(1, figsize = (2.5, 4.5))

# plot distribution of intron number
YMax1 = max(max([WithIntronCount.count(i) for i in WithIntronCount]),  max([IntronlessCount.count(i) for i in IntronlessCount])) + 10
ax1 = CreateAx(1, 3, 1, fig, [WithIntronCount, IntronlessCount], 'N external genes', 'Number of introns', np.arange(0, YMax1 + 10, 10), YMax1, np.arange(0, 80, 10), ['orange', 'blue'], 'histo')
# plot distribution of intron position
YMax2 = max(max([WithIntronPos.count(i) for i in WithIntronPos]), max([IntronlessPos.count(i) for i in IntronlessPos])) + 10
ax2 = CreateAx(1, 3, 2, fig, [WithIntronPos, IntronlessPos], 'N external genes', 'Intron position', np.arange(0, YMax2 + 20, 20), YMax2, np.arange(0, 60, 10), ['orange', 'blue'], 'histo')
# plot distribution of intron position / intron number
ax3 = CreateAx(1, 3, 3, fig, [WithIntronPosNorm, PWithPosNorm, IntronlessPosNorm, PNonePosNorm] , 'Probability', 'Intron position / intron number', np.arange(0, 1.1, 0.1), 1, np.arange(0, 1.1, 0.1), ['orange', 'blue'], 'cdf')

# make sure subplots do not overlap
plt.tight_layout()

## annotate graphs with legends
# create patches 
Title = mpatches.Patch(facecolor = 'none', edgecolor = 'none', linewidth = 0, label= 'Internal genes:', alpha = 0)
Yes = mpatches.Patch(facecolor = 'orange', edgecolor = 'none', linewidth = 0.5, label= 'with introns', alpha = 0.5)
No = mpatches.Patch(facecolor = 'blue', edgecolor = 'black', linewidth = 0.5, label= 'intronless', alpha = 0.5)
ax1.legend(handles = [Title, Yes, No], bbox_to_anchor=(0.2, 0.4), loc = 3, fontsize = 5, frameon = False, ncol = 1)
ax2.legend(handles = [Title, Yes, No], bbox_to_anchor=(0.2, 0.4), loc = 3, fontsize = 5, frameon = False, ncol = 1)
# create lines
orange_line = mlines.Line2D([], [], color='orange', marker='', markersize=15, label='with introns', alpha = 0.5)
blue_line = mlines.Line2D([], [], color='blue', marker='', markersize=15, label='intronless', alpha = 0.5)
ax3.legend(handles = [Title, orange_line, blue_line], bbox_to_anchor=(0.05, 0.5), loc = 3, fontsize = 5, frameon = False, ncol = 1)

# annotate graphs with p values
for i in range(len(PVals)):
    if PVals[i] != '':
        if i == 0:
            ax1.text(25, 20, PVals[i], color = 'black', size = 5, fontname = 'Arial')
        elif i == 1:
            ax2.text(15, 60, PVals[i], color = 'black', size = 5, fontname = 'Arial')
        elif i == 2:
            ax3.text(0.4, 0.3, PVals[i], color = 'black', size = 5, fontname = 'Arial')

# save figure
for extension in ['.pdf', '.eps', '.png']:
    fig.savefig('IntronPositionInternalGenes' + extension, bbox_inches = 'tight')
