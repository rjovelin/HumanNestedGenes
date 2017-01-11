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


# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)
NestedPairs = GetHostNestedPairs(Nested)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)


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
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseGTEXExpression('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# generate lists of gene pairs separated by distance 
Proximal, Moderate, Intermediate, Distant = GenerateSetsGenePairsDistance(GeneCoord, OrderedGenes, ExpressionProfile)

# filter gene pairs lacking expression
AllPairs = [OverlappingPairs, NestedPairs, PiggybackPairs, ConvergentPairs, DivergentPairs]
for i in range(len(AllPairs)):
    AllPairs[i] = FilterGenePairsWithoutExpression(AllPairs[i], ExpressionProfile)



#NestedPairs[i] = FilterGenePairsWithoutExpression(NestedPairs[i], SpExpression)

# compute expression divergence between pairs of genes
OverlapDiv = ComputeExpressionDivergenceGenePairs(AllPairs[0], ExpressionProfile)
NestedDiv = ComputeExpressionDivergenceGenePairs(AllPairs[1], ExpressionProfile)
PiggybackDiv = ComputeExpressionDivergenceGenePairs(AllPairs[2], ExpressionProfile)
ConvergentDiv = ComputeExpressionDivergenceGenePairs(AllPairs[3], ExpressionProfile)
DivergentDiv = ComputeExpressionDivergenceGenePairs(AllPairs[4], ExpressionProfile)
ProximalDiv = ComputeExpressionDivergenceGenePairs(Proximal, ExpressionProfile)
ModerateDiv = ComputeExpressionDivergenceGenePairs(Moderate, ExpressionProfile)
IntermediateDiv = ComputeExpressionDivergenceGenePairs(Intermediate, ExpressionProfile)
DistantDiv = ComputeExpressionDivergenceGenePairs(Distant, ExpressionProfile)


print('overlapping', len(OverlapDiv), np.median(OverlapDiv), np.mean(OverlapDiv))


print('nested', len(NestedDiv), np.median(NestedDiv), np.mean(NestedDiv))


print('piggyback', len(PiggybackDiv), np.median(PiggybackDiv), np.mean(PiggybackDiv))

print('convergent', len(ConvergentDiv), np.median(ConvergentDiv), np.mean(ConvergentDiv))

print('divergent', len(DivergentDiv), np.median(DivergentDiv), np.mean(DivergentDiv))


print('proximal', len(ProximalDiv), np.median(ProximalDiv), np.mean(ProximalDiv))
print('moderate', len(ModerateDiv), np.median(ModerateDiv), np.mean(ModerateDiv))
print('intermediate', len(IntermediateDiv), np.median(IntermediateDiv), np.mean(IntermediateDiv))
print('distant', len(DistantDiv), np.median(DistantDiv), np.mean(DistantDiv))




#
#
#
#
#
#
#
#
#
#
#
## create lists to store the list od distances in each species [[human], [chimp], [gorilla], [orang-outan], [macaque]]
#SpNestedDiv, SpProximalDiv, SpModerateDiv, SpIntermediateDiv, SpDistantDiv = [], [], [], [], []
#
#
#
#
#  
## store lists of expression divergence in parallel lists
#SpNestedDiv.append(NestedDiv)
#SpProximalDiv.append(ProximalDiv)
#SpModerateDiv.append(ModerateDiv)
#SpIntermediateDiv.append(IntermediateDiv)
#SpDistantDiv.append(DistantDiv)
#
#
#
###################
#
#
#
#
#
#
#
## make a list of all gene sets
#AllGeneSets = [NonOverlappingGenes, NestedGenes, PiggyBackGenes, ConvergentGenes, DivergentGenes]
#
## make a list of tissues
#infile = open('GTEX_Median_Normalized_FPKM.txt')
#header = infile.readline().rstrip().split('\t')
#Tissues = header[1:]
## replace spaces in tissue names
#for i in range(len(Tissues)):
#    if ' ' in Tissues[i]:
#        Tissues[i] = Tissues[i].replace(' ', '_')
#
#
## make a list of genes with expression
#Expressed = list(ExpressionProfile.keys())
## check that all genes have the same number of tissues
#for gene in Expressed:
#    assert len(ExpressionProfile[gene]) == len(Tissues)
#
#
## make a list of gene category names parallel to the list of gene sets
#GeneCats = ['NoOv', 'Nst', 'Pbk', 'Conv', 'Div']
#
## create a dict to count the number of genes in each category with highest expression in each tissue
#HighestExpression = {}
## inititialize dict with list of 0s
#for i in GeneCats:
#    HighestExpression[i] = [0] * len(Tissues)
#
## loop over the gene sets
#for i in range(len(AllGeneSets)):
#    # loop over each gene in given set
#    for gene in AllGeneSets[i]:
#        # check if gene has expression profile
#        if gene in ExpressionProfile:
#            # get the index of the maximum expression value
#            pos = ExpressionProfile[gene].index(max(ExpressionProfile[gene]))
#            # check that there is a single maximum expression value
#            assert ExpressionProfile[gene].count(max(ExpressionProfile[gene])) == 1, '> 1 max value'
#            # update counter at pos index
#            HighestExpression[GeneCats[i]][pos] += 1
#
## divide by the total number of genes in each category to get proportions
#for i in range(len(GeneCats)):
#    for j in range(len(HighestExpression[GeneCats[i]])):
#        HighestExpression[GeneCats[i]][j] = round(HighestExpression[GeneCats[i]][j] / len(AllGeneSets[i]), 4)
#
## create a dictionary with tissue as key and a list of gene proportions for each gene category as value
#Proportions = {}
#for i in range(len(Tissues)):
#    Proportions[Tissues[i]] = []
#    for j in range(len(GeneCats)):
#        Proportions[Tissues[i]].append(HighestExpression[GeneCats[j]][i])
#
#
## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YMax, YLabel):
#    '''
#    (int, int, int, list, figure_object, str, int, list, list)
#    Take the number of a column, and rows in the figure object and the position of
#    the ax in figure, a list of data, a title, a maximum value for the Y axis,
#    a list with species names and list of X axis tick positions and return an
#    ax instance in the figure
#    '''    
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#    # set colors
#    colorscheme = ['#fb9a99', '#a6cee3','#1f78b4','#b2df8a','#33a02c']
#    # plot nucleotide divergence
#    ax.bar([0, 0.2, 0.4, 0.6, 0.8], Data, 0.2, color = colorscheme,
#           edgecolor = 'black', linewidth = 0.5,
#           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write y axis label
#    if YLabel == True:
#        ax.set_ylabel('Proportion', color = 'black',  size = 7, ha = 'center', **FigFont)
#    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
#    
#    # add a range for the Y axis
#    plt.ylim([0, YMax])
#    
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(False)    
#    ax.spines["right"].set_visible(False)
#    if YLabel == True:
#        ax.spines["left"].set_visible(True)  
#    elif YLabel == False:
#        ax.spines['left'].set_visible(False)
#    
#    # edit tick parameters    
#    if YLabel == True:
#        plt.tick_params(axis='both', which='both', bottom='off', top='off',
#                        right = 'off', left = 'on', labelbottom='off', 
#                        colors = 'black', labelsize = 7, direction = 'out')  
#    elif YLabel == False:
#        plt.tick_params(axis='both', which='both', bottom='off', top='off',
#                        right = 'off', left = 'off', labelleft='off', labelbottom='off', 
#                        colors = 'black', labelsize = 7, direction = 'out')  
#    # Set the tick labels font name
#    for label in ax.get_yticklabels():
#        label.set_fontname('Arial')   
#    return ax      
#
#
## create figure
#fig = plt.figure(1, figsize = (5, 3))
#
## plot data
#j = 1
#for i in range(len(Tissues)):
#    tissue = Tissues[i]
#    if j == 1 or j == 11 or j == 21:
#        YLabel = True
#    else:
#        YLabel = False
#    ax = CreateAx(10, 3, j, fig, Proportions[tissue], tissue.lower().replace('_', '\n'), 0.31, YLabel)
#    j += 1
#
#    # create legend
#    if j == 10:
#        No = mpatches.Patch(facecolor = '#fb9a99', edgecolor = 'black', linewidth = 0.5, label= 'NoOv')
#        Ns = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 0.5, label= 'Nst')
#        Pk = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 0.5, label= 'Pbk')
#        Co = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 0.5, label= 'Conv')
#        Dv = mpatches.Patch(facecolor = '#33a02c', edgecolor = 'black', linewidth = 0.5, label= 'Div')
#        ax.legend(handles = [No, Ns, Pk, Co, Dv], bbox_to_anchor=(-10, 0.6), loc = 3, fontsize = 6, frameon = False, ncol = 5)
#        
#        
## adjust padding between subplots
## pad controls the padding around the figure border
## hpad and wpad control the padding between subplots
#plt.tight_layout(pad=0.2, w_pad=0, h_pad=0.5)
#
#
#
#
#
#
#
#
#
#
#
###########################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#   
## create lists of lists with expression divergence for the different gene pairs for each species
#HumanExpDiv = [SpNestedDiv[0], SpProximalDiv[0], SpModerateDiv[0], SpIntermediateDiv[0], SpDistantDiv[0]]
#ChimpExpDiv = [SpNestedDiv[1], SpProximalDiv[1], SpModerateDiv[1], SpIntermediateDiv[1], SpDistantDiv[1]]
#GorillaExpDiv = [SpNestedDiv[2], SpProximalDiv[2], SpModerateDiv[2], SpIntermediateDiv[2], SpDistantDiv[2]]
#OrangOutanExpDiv = [SpNestedDiv[3], SpProximalDiv[3], SpModerateDiv[3], SpIntermediateDiv[3], SpDistantDiv[3]]    
#MacaqueExpDiv = [SpNestedDiv[4], SpProximalDiv[4], SpModerateDiv[4], SpIntermediateDiv[4], SpDistantDiv[4]]
#
#
## create lists with means and with SEM
#HumanMeans, HumanSEM = [], []
#for i in range(len(HumanExpDiv)):
#    HumanMeans.append(np.mean(HumanExpDiv[i]))
#    HumanSEM.append(np.std(HumanExpDiv[i]) / math.sqrt(len(HumanExpDiv[i])))
#ChimpMeans, ChimpSEM = [], []
#for i in range(len(ChimpExpDiv)):
#    ChimpMeans.append(np.mean(ChimpExpDiv[i]))
#    ChimpSEM.append(np.std(ChimpExpDiv[i]) / math.sqrt(len(ChimpExpDiv[i])))
#GorillaMeans, GorillaSEM = [], []
#for i in range(len(GorillaExpDiv)):
#    GorillaMeans.append(np.mean(GorillaExpDiv[i]))
#    GorillaSEM.append(np.std(GorillaExpDiv[i]) / math.sqrt(len(GorillaExpDiv[i])))
#OrangOutanMeans, OrangOutanSEM = [], []
#for i in range(len(OrangOutanExpDiv)):
#    OrangOutanMeans.append(np.mean(OrangOutanExpDiv[i]))
#    OrangOutanSEM.append(np.std(OrangOutanExpDiv[i]) / math.sqrt(len(OrangOutanExpDiv[i])))
#MacaqueMeans, MacaqueSEM = [], []
#for i in range(len(MacaqueExpDiv)):
#    MacaqueMeans.append(np.mean(MacaqueExpDiv[i]))
#    MacaqueSEM.append(np.std(MacaqueExpDiv[i]) / math.sqrt(len(MacaqueExpDiv[i])))
#
#
## perform statistical tests between gene categories in all species
## create dict to store results {species: [P_host-nested, P_host-unnested, P_nested-unnested]}
#AllData = [HumanExpDiv, ChimpExpDiv, GorillaExpDiv, OrangOutanExpDiv, MacaqueExpDiv]
#PValues = {}
## loop over inner lists in data list
#for i in range(len(AllData)):
#    # initialize dict with empty list
#    PValues[species[i]] = []
#    # loop over inner list, compare gene categories
#    for j in range(0, len(AllData[i]) -1):
#        for k in range(j+1, len(AllData[i])):
#            P = stats.ranksums(AllData[i][j], AllData[i][k])[1]
#            PValues[species[i]].append(P)
## print p values
#for sp in PValues:
#    print(sp, PValues[sp])
#
#
## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, figure, Means, SEM, XLabel, YLabel, YMax, YAxis):
#    '''
#    (int, int, int, list, figure_object, str, int, list, list)
#    Take the number of a column, and rows in the figure object and the position of
#    the ax in figure, a list of data, a title, a maximum value for the Y axis,
#    a list with species names and list of X axis tick positions and return an
#    ax instance in the figure
#    '''    
#    # create subplot in figure
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#    # set colors
#    colorscheme = ['#1f78b4', '#edf8fb','#b2e2e2','#66c2a4','#238b45']    
#    # plot nucleotide divergence
#    ax.bar([0, 0.4, 0.8, 1.2, 1.6], Means, 0.3, yerr = SEM, color = colorscheme,
#           edgecolor = 'black', linewidth = 1,
#           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write label for y and x axis    
#    if YAxis == True:
#        ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
#    ax.set_xlabel(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
#        
#    # add a range for the Y axis
#    plt.ylim([0, YMax])
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)
#    if YAxis == True:
#        ax.spines["left"].set_visible(True)  
#    elif YAxis == False:
#        ax.spines["left"].set_visible(False)
#        
#    if YAxis == True:
#        # do not show ticks
#        plt.tick_params(
#            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#            which='both',      # both major and minor ticks are affected
#            bottom='off',      # ticks along the bottom edge are off
#            top='off',         # ticks along the top edge are off
#            right = 'off',
#            left = 'on',          
#            labelbottom='off', # labels along the bottom edge are on
#            colors = 'black',
#            labelsize = 8,
#            direction = 'out') # ticks are outside the frame when bottom = 'on'  
#    elif YAxis == False:
#        # do not show ticks
#        plt.tick_params(
#            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#            which='both',      # both major and minor ticks are affected
#            bottom='off',      # ticks along the bottom edge are off
#            top='off',         # ticks along the top edge are off
#            right = 'off',
#            left = 'off',          
#            labelbottom='off', # labels along the bottom edge are on
#            colors = 'black',
#            labelsize = 8,
#            labelleft = 'off',
#            direction = 'out') # ticks are outside the frame when bottom = 'on'  
#    
#    ## write label for x axis
#    #ax.set_xticks([0.15, 0.55, 0.95, 1.35, 1.75])
#    #ax.set_xticklabels(['Incl', 'Prox', 'Mod', 'Inter', 'Dist'], rotation = 0, ha = 'center', fontsize = 8, **FigFont)   
#    
#    if YAxis == True:
#        # Set the tick labels font name
#        for label in ax.get_yticklabels():
#            label.set_fontname('Arial')
#    # create a margin around the x axis
#    plt.margins(0.1)
#    return ax      
#
#
## create figure
#fig = plt.figure(1, figsize = (7, 2.5))
#
## plot data
#ax1 = CreateAx(5, 1, 1, fig, HumanMeans, HumanSEM, 'Human', 'Expression divergence', 0.71, True)
#ax2 = CreateAx(5, 1, 2, fig, ChimpMeans, ChimpSEM, 'Chimp', 'Expression divergence', 0.71, False)
#ax3 = CreateAx(5, 1, 3, fig, GorillaMeans, GorillaSEM, 'Gorilla', 'Expression divergence', 0.71, False)
#ax4 = CreateAx(5, 1, 4, fig, OrangOutanMeans, OrangOutanSEM, 'Orangutan', 'Expression divergence', 0.71, False)
#ax5 = CreateAx(5, 1, 5, fig, MacaqueMeans, MacaqueSEM, 'Macaque', 'Expression divergence', 0.71, False)
#
## annotate figure to add significance
## significant comparisons were already determined, just need to add letters to show significance
## get the x and y coordinates
#HumanDiff = ['A', 'B', 'C', 'AD', 'AE']
#ChimpDiff = ['A', 'B', 'AD', 'AE', 'C']
#GorillaDiff = ['A', 'B', 'C', 'A', 'D']
#OrangutanDiff = ['A', 'B', 'C', 'AD', 'AE']
#MacaqueDiff = ['A', 'B', 'AD', 'AE', 'C']
#ypos = [0.62] * 5
#xpos =  [0.15, 0.55, 0.95, 1.35, 1.75]
#for i in range(len(HumanDiff)):
#    ax1.text(xpos[i], ypos[i], HumanDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#    ax2.text(xpos[i], ypos[i], ChimpDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#    ax3.text(xpos[i], ypos[i], GorillaDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#    ax4.text(xpos[i], ypos[i], OrangutanDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#    ax5.text(xpos[i], ypos[i], MacaqueDiff[i], horizontalalignment='center',
#             verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#
## add legend relative to ax1 using ax1 coordinates
#N = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 1, label= 'Nested')
#P = mpatches.Patch(facecolor = '#edf8fb', edgecolor = 'black', linewidth = 1, label= 'Proximal')
#M = mpatches.Patch(facecolor = '#b2e2e2', edgecolor = 'black', linewidth = 1, label= 'Moderate')
#I = mpatches.Patch(facecolor = '#66c2a4', edgecolor = 'black', linewidth = 1, label= 'Intermediate')
#D = mpatches.Patch(facecolor = '#238b45', edgecolor = 'black', linewidth = 1, label= 'Distant')
#ax1.legend(handles = [N, P, M, I, D], loc = (0.5, 1), fontsize = 8, frameon = False, ncol = 5)
#
## make sure subplots do not overlap
#plt.tight_layout()
#
## save figure
#fig.savefig('ExpressionDivergenceDistance.pdf', bbox_inches = 'tight')
#fig.savefig('ExpressionDivergenceDistance.eps', bbox_inches = 'tight')
#
