# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 21:47:14 2016

@author: Richard
"""

# use this script to plot intron length of host genes for gene-containing introns
# and introns without genes

# usage python3 PlotIntronLengthHostGenes.py

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
import matplotlib.gridspec as gridspec
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
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('ChimpHostNestedGenes.json') as chimp_json_data:
    ChimpHostGenes = json.load(chimp_json_data)
with open('GorillaHostNestedGenes.json') as gorilla_json_data:
    GorillaHostGenes = json.load(gorilla_json_data)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'

# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF]
# make a list of species names
SpeciesNames = ['human', 'chimp', 'gorilla']
# make a list of host:nested genes dictionaries
HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes]

# create lists to store intron length of gene-containing and without genes for
# host genes in each speciest [[human], [chimp], [gorilla]]
ContainingIntron, NoGeneIntron = [], []

# loop over GFF files, find nested and intronic=nested genes in each species 
for i in range(len(GFFs)):
    print(GFFs[i][:GFFs[i].index('.')], SpeciesNames[i])
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    SpGeneChromoCoord = ChromoGenesCoord(GFFs[i])
    # map each gene to its mRNA transcripts
    SpMapGeneTranscript = GeneToTranscripts(GFFs[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    SpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(SpGeneChromoCoord, SpMapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    SpGeneCoord = FromChromoCoordToGeneCoord(SpGeneChromoCoord)
    # Map Transcript names to gene names {transcript: gene}
    SpMapTranscriptGene = TranscriptToGene(GFFs[i])
    # get the coordinates of all exons    
    SpExonCoord = GeneExonCoord(GFFs[i])
    SpExonCoord = CleanGeneFeatureCoord(SpExonCoord, SpMapTranscriptGene)
    # get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
    SpIntronCoord = GeneIntronCoord(SpExonCoord)
    SpIntronCoord = CleanGeneFeatureCoord(SpIntronCoord, SpMapTranscriptGene)
    # get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
    SpTranscriptCoordinates = TranscriptsCoord(GFFs[i])
    # map genes to their longest transcript {gene: longest_transcript}
    SpGeneLongestTranscript = LongestTranscript(SpTranscriptCoordinates, SpMapGeneTranscript)
    # match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
    SpMatches = MatchHostTranscriptWithNestedTranscript(HostGenes[i], SpMapGeneTranscript, SpGeneLongestTranscript, SpTranscriptCoordinates, SpIntronCoord)
    # get intron length for host transcripts
    WithGene, WithoutGene = CollectHostGeneIntronLength(SpMatches, SpTranscriptCoordinates, SpIntronCoord)
    # convert bp to Kbp
    WithGene = list(map(lambda x: x/1000, WithGene))
    WithoutGene = list(map(lambda x: x/1000, WithoutGene))
    # populate lists
    ContainingIntron.append(WithGene)
    NoGeneIntron.append(WithoutGene)

# create a list of lists with intron length for gene-containing introns and introns without genes
HumanIntron = [ContainingIntron[0], NoGeneIntron[0]]
ChimpIntron = [ContainingIntron[1], NoGeneIntron[1]]
GorillaIntron = [ContainingIntron[2], NoGeneIntron[2]]

# create a list of list with all data
AllData = [HumanIntron, ChimpIntron, GorillaIntron]

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
HumanMeans, HumanSEM = GetMeanSEM(HumanIntron)
ChimpMeans, ChimpSEM = GetMeanSEM(ChimpIntron)
GorillaMeans, GorillaSEM = GetMeanSEM(GorillaIntron)

# perform statistical tests between gene categories in all species
# create dict to store results {species: [P_host-nested, P_host-unnested, P_nested-unnested]}
PValues = {}
# loop over inner lists in data list
for i in range(len(AllData)):
    # initialize dict with empty list
    PValues[SpeciesNames[i]] = []
    # loop over inner list, compare gene categories
    for j in range(0, len(AllData[i]) -1):
        for k in range(j+1, len(AllData[i])):
            P = stats.ranksums(AllData[i][j], AllData[i][k])[1]
            PValues[SpeciesNames[i]].append(P)
# print p values
for sp in PValues:
    print(sp, PValues[sp])

# create a dict with significance level as stars
Significance = {}
for species in SpeciesNames:
    # initialize dict with empty list
    Significance[species] = [] 
    # get the significance level
    for pval in PValues[species]:
        if pval >= 0.05:
            Significance[species].append('')
        elif pval < 0.05 and pval >= 0.01:
            Significance[species].append('*')
        elif pval < 0.01 and pval >= 0.001:
            Significance[species].append('**')
        elif pval < 0.001:
            Significance[species].append('***')


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Means, SEM, XLabel, YLabel, YMax, YAxis):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['#a6cee3','#1f78b4']
    # plot nucleotide divergence
    ax.bar([0, 0.2], Means, 0.2, yerr = SEM, color = colorscheme,
           edgecolor = 'black', linewidth = 1,
           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y and x axis
    if YAxis == True:
        ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    
    # add a range for the Y axis
    plt.ylim([0, YMax])
       
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    if YAxis == True:
        ax.spines["left"].set_visible(True)  
    elif YAxis == False:
        ax.spines["left"].set_visible(False)
    
    if YAxis == True:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'on',          
            labelbottom='off', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 8,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
    elif YAxis == False:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='off', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 8,
            labelleft = 'off',
            direction = 'out') # ticks are outside the frame when bottom = 'on'      
     
    if YAxis == True:
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')   
     
    # create a margin around the x axis
    plt.margins(0.1)
    return ax      


# use this function to annotate the graph with significance levels
def AddSignificance(ax, SignificanceLevel, XLine1, XLine2, YLine, XText, YText):
    '''
    (ax, str, num, num, num, num, num) -> ax
    Take a matplotlib ax object, the significance level (as stars), the positions
    of the bracket and star and return the ax with annotated significance level
    '''
    ax.annotate("", xy=(XLine1, YLine), xycoords='data', xytext=(XLine2, YLine), textcoords='data',
                 arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
    # add stars for significance
    ax.text(XText, YText, SignificanceLevel, horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 6)
    return ax


# create figure
fig = plt.figure(1, figsize = (2.5, 1.5))

# plot data
ax1 = CreateAx(3, 1, 1, fig, HumanMeans, HumanSEM, 'Human', 'Intron length (Kbp)', 70, True)
ax2 = CreateAx(3, 1, 2, fig, ChimpMeans, ChimpSEM, 'Chimp', 'Intron length (Kbp)', 70, False)
ax3 = CreateAx(3, 1, 3, fig, GorillaMeans, GorillaSEM, 'Gorilla', 'Intron length (Kbp)', 70, False)

# make lists with bracket and star positions
XPos = [[0.1, 0.28, 68, 0.2, 73]]
    
# annotate figure to add significance
for i in range(len(Significance['human'])):
    if Significance['human'][i] != '':
        ax1 = AddSignificance(ax1, Significance['human'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
for i in range(len(Significance['chimp'])):
    if Significance['chimp'][i]  != '':
        ax2 = AddSignificance(ax2, Significance['chimp'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
for i in range(len(Significance['gorilla'])):
    if Significance['gorilla'][i] != '':
        ax3 = AddSignificance(ax3, Significance['gorilla'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])

# add legend relative to ax1 using ax1 coordinates
W = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'With genes')
Wo = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 1, label= 'Without genes')
ax1.legend(handles = [W, Wo], loc = (-0.4, 1.1), fontsize = 8, frameon = False, ncol = 2)

# make sure subplots do not overlap
#plt.tight_layout()
# one can control padding between subplots with w_pad and h_pad 
#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

# padding between subplots can also be controlled with gridspec
gs = gridspec.GridSpec(1, 5) # N rows and columns
gs.update(wspace=0.3, hspace=0) # set the spacing between axes. 

# save figure
fig.savefig('IntronLengthHostgenes.pdf', bbox_inches = 'tight')
fig.savefig('IntronLengthHostgenes.eps', bbox_inches = 'tight')

