# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 18:38:17 2016

@author: Richard
"""

# use this script to plot expression divergence between gene pairs as a function
# of the distance between the 2 genes


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

# make lists of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
print('host-gene pairs in human', len(HumanHostNestedPairs))
ChimpHostNestedPairs = GetHostNestedPairs(ChimpHostGenes)
print('host-gene pairs in chimp', len(ChimpHostNestedPairs))
GorillaHostNestedPairs = GetHostNestedPairs(GorillaHostGenes)
print('host-gene pairs in gorilla', len(GorillaHostNestedPairs))
OrangOutanHostNestedPairs = GetHostNestedPairs(OrangOutanHostGenes)
print('host-gene pairs in orang-outan', len(OrangOutanHostNestedPairs))
MacaqueHostNestedPairs = GetHostNestedPairs(MacaqueHostGenes)
print('host-gene pairs in macaque', len(MacaqueHostNestedPairs))

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
PabGFF = 'Pongo_abelii.PPYG2.86.gff3'
MmlGFF = 'Macaca_mulatta.Mmul_8.0.1.86.gff3'    


# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF, PabGFF, MmlGFF]
species = ['human', 'chimp', 'gorilla', 'orangoutan', 'macaque']

# create lists to store the list od distances in each species [[human], [chimp], [gorilla], [orang-outan], [macaque]]
SpNestedDiv, SpProximalDiv, SpModerateDiv, SpIntermediateDiv, SpDistantDiv = [], [], [], [], []

# make a list of nested gene pairs lists
NestedPairs = [HumanHostNestedPairs, ChimpHostNestedPairs, GorillaHostNestedPairs, OrangOutanHostNestedPairs, MacaqueHostNestedPairs]

# loop over GFF files, find nested and intronic=nested genes in each species 
for i in range(len(GFFs)):
    print(GFFs[i][:GFFs[i].index('.')], species[i])
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
    SpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', species[i])
    # remove genes wuthout expression
    SpExpression = RemoveGenesLackingExpression(SpExpression)
    # get relative expression
    SpExpression = TransformRelativeExpression(SpExpression)
    # generate lists of gene pairs separated by distance 
    Proximal, Moderate, Intermediate, Distant = GenerateSetsGenePairsDistance(SpGeneCoord, SpOrderedGenes, SpExpression)
    # filter host-nested nested pairs based on expression
    NestedPairs[i] = FilterGenePairsWithoutExpression(NestedPairs[i], SpExpression)
    # compute expression divergence between pairs of genes
    NestedDiv = ComputeExpressionDivergenceGenePairs(NestedPairs[i], SpExpression)
    ProximalDiv = ComputeExpressionDivergenceGenePairs(Proximal, SpExpression)
    ModerateDiv = ComputeExpressionDivergenceGenePairs(Moderate, SpExpression)
    IntermediateDiv = ComputeExpressionDivergenceGenePairs(Intermediate, SpExpression)
    DistantDiv = ComputeExpressionDivergenceGenePairs(Distant, SpExpression)

    print('nested', len(NestedDiv), np.median(NestedDiv), np.mean(NestedDiv))
    print('proximal', len(ProximalDiv), np.median(ProximalDiv), np.mean(ProximalDiv))
    print('moderate', len(ModerateDiv), np.median(ModerateDiv), np.mean(ModerateDiv))
    print('intermediate', len(IntermediateDiv), np.median(IntermediateDiv), np.mean(IntermediateDiv))
    print('distant', len(DistantDiv), np.median(DistantDiv), np.mean(DistantDiv))
    
    # store lists of expression divergence in parallel lists
    SpNestedDiv.append(NestedDiv)
    SpProximalDiv.append(ProximalDiv)
    SpModerateDiv.append(ModerateDiv)
    SpIntermediateDiv.append(IntermediateDiv)
    SpDistantDiv.append(DistantDiv)
    
# create lists of lists with expression divergence for the different gene pairs for each species
HumanExpDiv = [SpNestedDiv[0], SpProximalDiv[0], SpModerateDiv[0], SpIntermediateDiv[0], SpDistantDiv[0]]
ChimpExpDiv = [SpNestedDiv[1], SpProximalDiv[1], SpModerateDiv[1], SpIntermediateDiv[1], SpDistantDiv[1]]
GorillaExpDiv = [SpNestedDiv[2], SpProximalDiv[2], SpModerateDiv[2], SpIntermediateDiv[2], SpDistantDiv[2]]
OrangOutanExpDiv = [SpNestedDiv[3], SpProximalDiv[3], SpModerateDiv[3], SpIntermediateDiv[3], SpDistantDiv[3]]    
MacaqueExpDiv = [SpNestedDiv[4], SpProximalDiv[4], SpModerateDiv[4], SpIntermediateDiv[4], SpDistantDiv[4]]


# create lists with means and with SEM
HumanMeans, HumanSEM = [], []
for i in range(len(HumanExpDiv)):
    HumanMeans.append(np.mean(HumanExpDiv[i]))
    HumanSEM.append(np.std(HumanExpDiv[i]) / math.sqrt(len(HumanExpDiv[i])))
ChimpMeans, ChimpSEM = [], []
for i in range(len(ChimpExpDiv)):
    ChimpMeans.append(np.mean(ChimpExpDiv[i]))
    ChimpSEM.append(np.std(ChimpExpDiv[i]) / math.sqrt(len(ChimpExpDiv[i])))
GorillaMeans, GorillaSEM = [], []
for i in range(len(GorillaExpDiv)):
    GorillaMeans.append(np.mean(GorillaExpDiv[i]))
    GorillaSEM.append(np.std(GorillaExpDiv[i]) / math.sqrt(len(GorillaExpDiv[i])))
OrangOutanMeans, OrangOutanSEM = [], []
for i in range(len(OrangOutanExpDiv)):
    OrangOutanMeans.append(np.mean(OrangOutanExpDiv[i]))
    OrangOutanSEM.append(np.std(OrangOutanExpDiv[i]) / math.sqrt(len(OrangOutanExpDiv[i])))
MacaqueMeans, MacaqueSEM = [], []
for i in range(len(MacaqueExpDiv)):
    MacaqueMeans.append(np.mean(MacaqueExpDiv[i]))
    MacaqueSEM.append(np.std(MacaqueExpDiv[i]) / math.sqrt(len(MacaqueExpDiv[i])))




## perform statistical tests between gene categories in all 3 species
## create dict to store results {species: [P_host-nested, P_host-unnested, P_nested-unnested]}
#PValues = {}
## loop over inner lists in data list
#for i in range(len(AllData)):
#    # initialize dict with empty list
#    PValues[SpeciesNames[i]] = []
#    # loop over inner list, compare gene categories
#    for j in range(0, len(AllData[i]) -1):
#        for k in range(j+1, len(AllData[i])):
#            P = stats.ranksums(AllData[i][j], AllData[i][k])[1]
#            PValues[SpeciesNames[i]].append(P)
## print p values
#for species in PValues:
#    print(species, PValues[species])


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
    colorscheme = ['#1f78b4', '#edf8fb','#b2e2e2','#66c2a4','#238b45']    
    # plot nucleotide divergence
    ax.bar([0, 0.4, 0.8, 1.2, 1.6], Means, 0.3, yerr = SEM, color = colorscheme,
           edgecolor = 'black', linewidth = 1,
           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    
    if YAxis == True:
        # write label for y and x axis
        ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    
    
    #ax.set_xlabel(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
        
    plt.title(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    
    
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
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 5))
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    if YAxis == True:
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
    elif YAxis == False:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='on', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 8,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
    
    
    
    
    # write label for x axis
    ax.set_xticks([0.15, 0.55, 0.95, 1.35, 1.75])
    ax.set_xticklabels(['Incl', 'Prox', 'Mod', 'Inter', 'Dist'], rotation = 0, ha = 'center', fontsize = 8, **FigFont)   
     
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    # create a margin around the x axis
    plt.margins(0.1)
    return ax      


# create figure
fig = plt.figure(1, figsize = (6, 3))

# plot data
ax1 = CreateAx(5, 1, 1, fig, HumanMeans, HumanSEM, 'Human', 'Expression\ndivergence', 0.8, True)
ax2 = CreateAx(5, 1, 2, fig, ChimpMeans, ChimpSEM, 'Chimp', 'Expression\ndivergence', 0.8, False)
ax3 = CreateAx(5, 1, 3, fig, GorillaMeans, GorillaSEM, 'Gorilla', 'Expression\ndivergence', 0.8, False)
ax4 = CreateAx(5, 1, 4, fig, OrangOutanMeans, OrangOutanSEM, 'Orangutan', 'Expression\ndivergence', 0.8, False)
ax5 = CreateAx(5, 1, 5, fig, MacaqueMeans, MacaqueSEM, 'Macaque', 'Expression\ndivergence', 0.8, False)

## annotate figure to add significance
## significant comparisons were already determined, just need to add letters to show significance
## get the x and y coordinates
#CelDiff = ['A', 'B', 'C', 'D', 'A']
#CbrDiff = ['A', 'B', 'C', 'A', 'D']
#CrmDiff = ['A', 'B', 'B', 'C', 'A']
#ypos = [19] * 5
#xpos =  [0.15, 0.55, 0.95, 1.35, 1.75]
#for i in range(len(CelDiff)):
#    ax1.text(xpos[i], ypos[i], CelDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#    ax2.text(xpos[i], ypos[i], CbrDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#    ax3.text(xpos[i], ypos[i], CrmDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)

# make sure subplots do not overlap
plt.tight_layout()

outputfile = 'truc'
# save figure
fig.savefig(outputfile + '.pdf', bbox_inches = 'tight')
#fig.savefig(outputfile + '.png', bbox_inches = 'tight')

