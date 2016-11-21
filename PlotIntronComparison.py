# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:41:55 2016

@author: RJovelin
"""


# use this script to plot the number of introns, or the intron length between host,
# nested and un-nested genes

# usage python3 PlotIntronComparison.py [option]:
# -['IntronNumber'/'Intronlength']: compare intron number or intron length

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


# get the variable to record from the comand line
PlottingVariable = sys.argv[1]


# load dictionaries of host and nested genes 
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('ChimpHostNestedGenes.json') as chimp_json_data:
    ChimpHostGenes = json.load(chimp_json_data)
with open('GorillaHostNestedGenes.json') as gorilla_json_data:
    GorillaHostGenes = json.load(gorilla_json_data)
with open('OrangOutanHostNestedGenes.json') as orangoutan_json_data:
    OrangOutanHostGenes = json.load(orangoutan_json_data)
with open('MacaqueHostNestedGenes.json') as macaque_json_data:
    MacaqueHostGenes = json.load(macaque_json_data)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
PabGFF = 'Pongo_abelii.PPYG2.86.gff3'
MmlGFF = 'Macaca_mulatta.Mmul_8.0.1.86.gff3' 

# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF, PabGFF, MmlGFF]
# make a list of species names
SpeciesNames = ['human', 'chimp', 'gorilla', 'orangoutan', 'macaque']
# make a list of host:nested genes dictionaries
HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes, OrangOutanHostGenes, MacaqueHostGenes]


# create lists to store the number of introns for the longest transcript of the host,
# nested and un-nested genes in each species [[human], [chimp], [gorilla], [orang-outan], [macaque]]
HostIntron, NestedIntron, OthersIntron = [], [], []


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
    # get expression profile of the species genes
    SpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[i])
    # remove genes wuthout expression
    SpExpression = RemoveGenesLackingExpression(SpExpression)
    # get relative expression
    SpExpression = TransformRelativeExpression(SpExpression)
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
    # make a set of un-nested genes
    SpUnNestedGenes = MakeHostNestedGeneSet(HostGenes[i])
    
    # check if intron number or intron length is to be recorded
    if PlottingVariable == 'IntronNumber':
        # get the number of introns for the longest transcript of the host, nested and un-nested genes 
        HostNum, NestedNum, OthersNum = CollectIntronNumbers(SpUnNestedGenes, SpMatches, SpGeneLongestTranscript, SpTranscriptCoordinates, SpIntronCoord)
    elif PlottingVariable == 'IntronLength':
        # get intron length for host transcripts
        GeneContainingIntron, NoGeneIntron = CollectHostGeneIntronLength(SpMatches, SpTranscriptCoordinates, SpIntronCoord)
        HostNum = GeneContainingIntron + NoGeneIntron
        # get intron length for nested transcripts
        NestedNum = CollectNestedGeneIntronLength(SpMatches, SpTranscriptCoordinates, SpIntronCoord)
        # get intron length for un-nested genes
        OthersNum = CollectUnNestedGeneIntronLength(SpUnNestedGenes,  SpGeneLongestTranscript, SpTranscriptCoordinates, SpIntronCoord)
    # populate lists
    HostIntron.append(HostNum)
    NestedIntron.append(NestedNum)
    OthersIntron.append(OthersNum)
    
# create a list of lists with intron numbers/length for hosts, nested and un-nested genes for each species
HumanIntron = [HostIntron[0], NestedIntron[0], OthersIntron[0]]
ChimpIntron = [HostIntron[1], NestedIntron[1], OthersIntron[1]]
GorillaIntron = [HostIntron[2], NestedIntron[2], OthersIntron[2]]
OrangutanIntron = [HostIntron[3], NestedIntron[3], OthersIntron[3]]
MacaqueIntron = [HostIntron[4], NestedIntron[4], OthersIntron[4]]

# create a list of list with all data
AllData = [HumanIntron, ChimpIntron, GorillaIntron, OrangutanIntron, MacaqueIntron]


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
OrangOutanMeans, OrangOutanSEM = GetMeanSEM(OrangutanIntron)
MacaqueMeans, MacaqueSEM = GetMeanSEM(MacaqueIntron)

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
def CreateAx(Columns, Rows, Position, figure, Means, SEM, XLabel, YLabel, YMax):
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
    colorscheme = ['#a6cee3','#1f78b4','#b2df8a']
    # plot nucleotide divergence
    ax.bar([0, 0.2, 0.4], Means, 0.2, yerr = SEM, color = colorscheme,
           edgecolor = 'black', linewidth = 1,
           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y and x axis
    ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, style = 'italic', color = 'black',  size = 8, ha = 'center', **FigFont)
    
    # add a range for the Y axis
    plt.ylim([0, YMax])
       
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(True)  
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 5))
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)

    # do not show ticks
    plt.tick_params(
        axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        right = 'off',
        left = 'on',          
        labelbottom='on', # labels along the bottom edge are on
        colors = 'black',
        labelsize = 8,
        direction = 'out') # ticks are outside the frame when bottom = 'on'  
     
    # write label for x axis
    ax.set_xticks([0.1, 0.3, 0.5])
    ax.set_xticklabels(['Hosts', 'Nested', 'Others'], rotation = 30, ha = 'right', fontsize = 8, **FigFont)   
     
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
fig = plt.figure(1, figsize = (5, 3))

# plot data
if PlottingVariable == 'IntronNumber':
    ax1 = CreateAx(5, 1, 1, fig, HumanMeans, HumanSEM, 'Human', 'Number of introns per gene', 50)
    ax2 = CreateAx(5, 1, 2, fig, ChimpMeans, ChimpSEM, 'Chimp', 'Number of introns per gene', 50)
    ax3 = CreateAx(5, 1, 3, fig, GorillaMeans, GorillaSEM, 'Gorilla', 'Number of introns per gene', 50)
    ax4 = CreateAx(5, 1, 4, fig, OrangutanMeans, OrangutanSEM, 'Orangutan', 'Number of introns per gene', 50)
    ax5 = CreateAx(5, 1, 5, fig, MacaqueMeans, MacaqueSEM, 'Macaque', 'Number of introns per gene', 50)





#elif PlottingVariable == 'IntronLength':
#    ax1 = CreateAx(3, 1, 1, fig, CelMeans, CelSEM, 'C. elegans', 'Intron length (bp)', 800)
#    ax2 = CreateAx(3, 1, 2, fig, CbrMeans, CbrSEM, 'C. briggsae', 'Intron length (bp)', 1000)
#    ax3 = CreateAx(3, 1, 3, fig, CrmMeans, CrmSEM, 'C. remanei', 'Intron length (bp)', 600)

## make lists with bracket and star positions
#if PlottingVariable == 'IntronNumber':
#    CelPos = [[0.1, 0.28, 11.5, 0.2, 12], [0.1, 0.5, 13, 0.3, 14], [0.32, 0.5, 11.5, 0.4, 12]]
#    CbrPos = [[0.1, 0.28, 11.5, 0.2, 12], [0.1, 0.5, 13, 0.3, 14], [0.32, 0.5, 11.5, 0.4, 12]]
#    CrmPos = [[0.1, 0.28, 12.8, 0.2, 13.3], [0.1, 0.5, 14, 0.3, 15], [0.32, 0.5, 12.8, 0.4, 13.3]]
#elif PlottingVariable == 'IntronLength':
#    CelPos = [[0.1, 0.28, 720, 0.2, 740], [0.1, 0.5, 760, 0.3, 800], [0.32, 0.5, 720, 0.4, 740]]
#    CbrPos = [[0.1, 0.28, 940, 0.2, 970], [0.1, 0.5, 980, 0.3, 1030], [0.32, 0.5, 940, 0.4, 970]]
#    CrmPos = [[0.1, 0.28, 460, 0.2, 480], [0.1, 0.5, 490, 0.3, 530], [0.32, 0.5, 460, 0.4, 480]]
#
## annotate figure to add significance
#for i in range(len(Significance['Cel'])):
#    if Significance['Cel'][i] != '':
#        ax1 = AddSignificance(ax1, Significance['Cel'][i], CelPos[i][0], CelPos[i][1], CelPos[i][2], CelPos[i][3], CelPos[i][4])
#for i in range(len(Significance['Cbr'])):
#    if Significance['Cbr'][i]  != '':
#        ax2 = AddSignificance(ax2, Significance['Cbr'][i], CbrPos[i][0], CbrPos[i][1], CbrPos[i][2], CbrPos[i][3], CbrPos[i][4])
#for i in range(len(Significance['Crm'])):
#    if Significance['Crm'][i] != '':
#        ax3 = AddSignificance(ax3, Significance['Crm'][i], CrmPos[i][0], CrmPos[i][1], CrmPos[i][2], CrmPos[i][3], CrmPos[i][4])
    
 
# make sure subplots do not overlap
plt.tight_layout()

## build outputfile with arguments
#if PlottingVariable == 'IntronNumber':
#    outputfile = 'IntronNumberDiffAllGenes'
#elif PlottingVariable == 'IntronLength':
#    outputfile = 'IntronLengthDiffAllGenes'
# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
#fig.savefig(outputfile + '.png', bbox_inches = 'tight')

