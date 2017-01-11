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

# add gene pairs to Allpairs list
AllPairs.append(Proximal)
AllPairs.append(Moderate)
AllPairs.append(Intermediate)
AllPairs.append(Distant)


# compute expression divergence between pairs of genes
Divergence = []
for i in range(len(AllPairs)):
    Div = ComputeExpressionDivergenceGenePairs(AllPairs[i], ExpressionProfile)
    Divergence.append(Div)

#
## compute expression divergence between pairs of genes
#OverlapDiv = ComputeExpressionDivergenceGenePairs(AllPairs[0], ExpressionProfile)
#NestedDiv = ComputeExpressionDivergenceGenePairs(AllPairs[1], ExpressionProfile)
#PiggybackDiv = ComputeExpressionDivergenceGenePairs(AllPairs[2], ExpressionProfile)
#ConvergentDiv = ComputeExpressionDivergenceGenePairs(AllPairs[3], ExpressionProfile)
#DivergentDiv = ComputeExpressionDivergenceGenePairs(AllPairs[4], ExpressionProfile)
#ProximalDiv = ComputeExpressionDivergenceGenePairs(Proximal, ExpressionProfile)
#ModerateDiv = ComputeExpressionDivergenceGenePairs(Moderate, ExpressionProfile)
#IntermediateDiv = ComputeExpressionDivergenceGenePairs(Intermediate, ExpressionProfile)
#DistantDiv = ComputeExpressionDivergenceGenePairs(Distant, ExpressionProfile)


# make a list of gene category names parallel to the list of gene pairs
GeneCats = ['NoOv', 'Nst', 'Pbk', 'Conv', 'Div', 'Prox', 'Mod', 'Int', 'Dist']


for i in range(len(Divergence)):
    print(GeneCats[i], len(Divergence[i]), np.median(Divergence[i]), np.mean(Divergence[i]))



#print('overlapping', len(OverlapDiv), np.median(OverlapDiv), np.mean(OverlapDiv))
#print('nested', len(NestedDiv), np.median(NestedDiv), np.mean(NestedDiv))
#print('piggyback', len(PiggybackDiv), np.median(PiggybackDiv), np.mean(PiggybackDiv))
#print('convergent', len(ConvergentDiv), np.median(ConvergentDiv), np.mean(ConvergentDiv))
#print('divergent', len(DivergentDiv), np.median(DivergentDiv), np.mean(DivergentDiv))
#print('proximal', len(ProximalDiv), np.median(ProximalDiv), np.mean(ProximalDiv))
#print('moderate', len(ModerateDiv), np.median(ModerateDiv), np.mean(ModerateDiv))
#print('intermediate', len(IntermediateDiv), np.median(IntermediateDiv), np.mean(IntermediateDiv))
#print('distant', len(DistantDiv), np.median(DistantDiv), np.mean(DistantDiv))



## store lists with divergence values in a single list
#Divergence = [OverlapDiv, NestedDiv, PiggybackDiv, ConvergentDiv, DivergentDiv,
#              ProximalDiv, ModerateDiv, IntermediateDiv, DistantDiv]

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

# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = GetMeanSEM(Divergence)







# create figure
fig = plt.figure(1, figsize = (3, 2))

# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# set colors
#colorscheme = ['#fb9a99', '#a6cee3','#1f78b4','#b2df8a','#33a02c', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
colorscheme = ['grey', 'grey','grey','grey','grey', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']


#No = mpatches.Patch(facecolor = '#fb9a99', edgecolor = 'black', linewidth = 0.5, label= 'NoOv')
#Ns = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 0.5, label= 'Nst')
#Pk = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 0.5, label= 'Pbk')
#Co = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 0.5, label= 'Conv')
#Dv = mpatches.Patch(facecolor = '#33a02c', edgecolor = 'black', linewidth = 0.5, label= 'Div')
#ax.legend(handles = [No, Ns, Pk, Co, Dv], bbox_to_anchor=(-10, 0.6), loc = 3, fontsize = 6, frameon = False, ncol = 5)


# plot nucleotide divergence
ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
       edgecolor = 'black', linewidth = 0.5,
       error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))

# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)

# add ticks and lebels
plt.xticks([0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9, 2.2, 2.5], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)


# add a range for the Y and X axes
plt.ylim([0, 0.6])
plt.xlim([0, 2.7])


  
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)  

   
# edit tick parameters    
plt.tick_params(axis='both', which='both', bottom='on', top='off',
                right = 'off', left = 'on', labelbottom='on',
                colors = 'black', labelsize = 7, direction = 'out')  

# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   
      






##################


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


# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

