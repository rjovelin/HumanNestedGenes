# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:14:49 2016

@author: RJovelin
"""

# use this script to plot the proportion of host, nested and control genes with
# with highest expression in each tissue

# usage python3 PlotHighestExpression.py 


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

# make parallel lists to store expression specificity of host and nested genes [[human], [chimp], [gorilla], [orangutan], [macaque]]
HostSpecificity, NestedSpecificity, ControlSpecificity = [], [], []


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
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    SpOrderedGenes = OrderGenesAlongChromo(SpGeneChromoCoord)
    # get expression profile of the species genes
    SpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[i])
    # remove genes wuthout expression
    SpExpression = RemoveGenesLackingExpression(SpExpression)
    # get relative expression
    SpExpression = TransformRelativeExpression(SpExpression)
    # compute expression specificity for all genes
    SpTau = ExpressionSpecificity(SpExpression)    
    # make a set of host and nested genes (include all host and nested even if not expressed)    
    SpNestedConformation = MakeHostNestedGeneSet(HostGenes[i])    
    # generate a dict of expressed genes on each chromo to randomly draw from
    SpGenesToDrawFrom = GenerateAllUnNestedGenes(SpNestedConformation, SpOrderedGenes, SpExpression)
    # make a list of host-nested gene pairs
    SpHostNestedPairs = GetHostNestedPairs(HostGenes[i])
    # remove gene pairs with genes lacking expression
    SpHostNestedPairs = FilterGenePairsWithoutExpression(SpHostNestedPairs, SpExpression)
    # create a list of control genes matching chromosomes of the host and nested pairs
    SpControl = []
    # loop over pairs of host and nested genes
    for pair in SpHostNestedPairs:
        # get the chromosome of the host and nested genes
        chromo = SpGeneCoord[pair[0]][0]
        assert chromo == SpGeneCoord[pair[1]][0], 'chromosome of host and nested genes do not match'        
        # draw 2 genes at random on chromo to match tthe host and nested genes
        for j in range(2):
            k = random.randint(0, len(SpGenesToDrawFrom[chromo]) -1)
            SpControl.append(SpGenesToDrawFrom[chromo][k])
    # create lists to store specificity for host and nested genes for that species
    tauhost, taunested, taucontrol = [], [], []
    # loop over gene pairs:
    for pair in SpHostNestedPairs:
        # populate lists
        tauhost.append(SpTau[pair[0]])
        taunested.append(SpTau[pair[1]])
    # loop over control genes
    for control in SpControl:
        taucontrol.append(SpTau[control])
    HostSpecificity.append(tauhost)
    NestedSpecificity.append(taunested)
    ControlSpecificity.append(taucontrol)

# make lists of gene specificity for each species     
HumanTau = [HostSpecificity[0], NestedSpecificity[0], ControlSpecificity[0]]
ChimpTau = [HostSpecificity[1], NestedSpecificity[1], ControlSpecificity[1]]
GorillaTau = [HostSpecificity[2], NestedSpecificity[2], ControlSpecificity[2]]
OrangutanTau = [HostSpecificity[3], NestedSpecificity[3], ControlSpecificity[3]]
MacaqueTau = [HostSpecificity[4], NestedSpecificity[4], ControlSpecificity[4]]

print(list(map(lambda x: len(x), HumanTau)))
print(list(map(lambda x: len(x), ChimpTau)))
print(list(map(lambda x: len(x), GorillaTau)))
print(list(map(lambda x: len(x), OrangutanTau)))
print(list(map(lambda x: len(x), MacaqueTau)))






# make a list with all the list data
AllData = [HumanTau, ChimpTau, GorillaTau, OrangutanTau, MacaqueTau]

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
HumanMeans, HumanSEM = GetMeanSEM(HumanTau)
ChimpMeans, ChimpSEM = GetMeanSEM(ChimpTau)
GorillaMeans, GorillaSEM = GetMeanSEM(GorillaTau)
OrangutanMeans, OrangutanSEM = GetMeanSEM(OrangutanTau)
MacaqueMeans, MacaqueSEM = GetMeanSEM(MacaqueTau)

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

print(HumanMeans)
print(ChimpMeans)
print(GorillaMeans)
print(OrangutanMeans)
print(MacaqueMeans)


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
    colorscheme = ['#a6cee3','#1f78b4','#b2df8a']
    # plot nucleotide divergence
    ax.bar([0, 0.2, 0.4], Means, 0.2, yerr = SEM, color = colorscheme,
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
fig = plt.figure(1, figsize = (5, 2.5))

# plot data
ax1 = CreateAx(5, 1, 1, fig, HumanMeans, HumanSEM, 'Human', 'Expression specificity', 1, True)
ax2 = CreateAx(5, 1, 2, fig, ChimpMeans, ChimpSEM, 'Chimp', 'Expression specificity', 1, False)
ax3 = CreateAx(5, 1, 3, fig, GorillaMeans, GorillaSEM, 'Gorilla', 'Expression specificity', 1, False)
ax4 = CreateAx(5, 1, 4, fig, OrangutanMeans, OrangutanSEM, 'Orangutan', 'Expression specificity', 1, False)
ax5 = CreateAx(5, 1, 5, fig, MacaqueMeans, MacaqueSEM, 'Macaque', 'Expression specificity', 1, False)

# make lists with bracket and star positions
XPos = [[0.1, 0.28, 0.8, 0.2, 0.85], [0.1, 0.5, 0.9, 0.3, 0.95], [0.32, 0.5, 0.8, 0.4, 0.85]]

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
for i in range(len(Significance['orangoutan'])):
    if Significance['orangoutan'][i] != '':
        ax4 = AddSignificance(ax4, Significance['orangoutan'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
for i in range(len(Significance['macaque'])):
    if Significance['macaque'][i] != '':
        ax5 = AddSignificance(ax5, Significance['macaque'][i], XPos[i][0], XPos[i][1], XPos[i][2], XPos[i][3], XPos[i][4])
 

# add legend relative to ax1 using ax1 coordinates
H = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'Hosts')
N = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 1, label= 'Nested')
U = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'Control')
ax1.legend(handles = [H, N, U], loc = (0.5, 1), fontsize = 8, frameon = False, ncol = 3)

# make sure subplots do not overlap
plt.tight_layout()

fig.savefig('truc.pdf', bbox_inches = 'tight')
#fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
