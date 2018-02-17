# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 12:31:41 2017

@author: RJovelin
"""

# use this script to plot expression divergence between gene pairs as a function
# of the distance between the 2 genes in mouse and chimp


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


# load dictionaries of the 4 types of overlapping genes
JsonFiles = ['Nested', 'PiggyBack', 'Convergent', 'Divergent']
# make a list of dictionaries
Overlap = []
# loop over files
for species in ['Mouse', 'Chimp']:
    OverlappingGenes = []
    for i in range(len(JsonFiles)):
        # load dictionary of overlapping gene pairs
        json_data = open(species + JsonFiles[i] + 'Genes.json')
        overlapping = json.load(json_data)
        json_data.close()
        OverlappingGenes.append(overlapping)
    Overlap.append(OverlappingGenes)
    
# create a list of list of gene pairs
AllPairs = []
for i in range(len(Overlap)):
    sppairs = []
    # loop over dicts
    for j in range(len(Overlap[i])):
        sppairs.append(GetHostNestedPairs(Overlap[i][j]))
    AllPairs.append(sppairs)

# get GFF file
GFF = ['Mus_musculus.GRCm38.88.gff3', 'Pan_troglodytes.CHIMP2.1.4.88.gff3']
# get gene coordinates and oreded genes along chromosome
SpCoord, SpOrdered = [], []
for i in range(len(GFF)):
    # get the coordinates of genes on each chromo {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    SpCoord.append(GeneCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    SpOrdered.append(OrderedGenes)    
print('extracted gene coordinates')

    
# parse expression data
ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
MouseExpression = ParseExpressionFile('Mouse_Median_Normalized_FPKM.txt')
# remove genes without any expression
ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
MouseExpression = RemoveGenesLackingExpression(MouseExpression)
# transform absulte expression in relative expression
ChimpExpression = TransformRelativeExpression(ChimpExpression)
MouseExpression = TransformRelativeExpression(MouseExpression)
print('got gene expression')

# generate lists of gene pairs separated by distance in mouse
MmuProximal, MmuModerate, MmuIntermediate, MmuDistant = GenerateSetsGenePairsDistance(SpCoord[0], SpOrdered[0], MouseExpression)
# generate lists of gene pairs separated by distance in mouse
PtrProximal, PtrModerate, PtrIntermediate, PtrDistant = GenerateSetsGenePairsDistance(SpCoord[1], SpOrdered[1], ChimpExpression)
print('generated pairs of distant genes')


# remove gene pairs if any gene in the pair lacks expression
for i in range(len(AllPairs)):
    for j in range(len(AllPairs[i])):
        if i == 0:
            AllPairs[i][j] = FilterGenePairsWithoutExpression(AllPairs[i][j], MouseExpression, 'strict')
        elif i == 1:
            AllPairs[i][j] = FilterGenePairsWithoutExpression(AllPairs[i][j], ChimpExpression, 'strict')
print('removed pairs lacking expression') 
 
# make list of gene pairs for each gene categories and each species
MousePairs, ChimpPairs = AllPairs[0], AllPairs[1]            
# add gene pairs to overlapping gene pairs            
for L in [MmuProximal, MmuModerate, MmuIntermediate, MmuDistant]:
    MousePairs.append(L)
for L in [PtrProximal, PtrModerate, PtrIntermediate, PtrDistant]:
    ChimpPairs.append(L)
print('generate gene pairs in chimp and mouse')


# compute expression divergence between pairs of genes
MmuDivergence, PtrDivergence = [], []
for i in range(len(MousePairs)):
    Div = ComputeExpressionDivergenceGenePairs(MousePairs[i], MouseExpression)
    MmuDivergence.append(Div)
for i in range(len(ChimpPairs)):
    Div = ComputeExpressionDivergenceGenePairs(ChimpPairs[i], ChimpExpression)
    PtrDivergence.append(Div)
print('computed expression divergence')


# make a list of gene category names parallel to the list of gene pairs
GeneCats = ['Nst', 'Pbk', 'Conv', 'Div', 'Prox', 'Mod', 'Int', 'Dist']

# create lists with means and SEM for each gene category
MmuMeanExpDiv, MmuSEMExpDiv, PtrMeanExpDiv, PtrSEMExpDiv = [], [], [], []
# loop over lists in Divergence list
for i in range(len(MmuDivergence)):
    MmuMeanExpDiv.append(np.mean(MmuDivergence[i]))
    MmuSEMExpDiv.append(np.std(MmuDivergence[i]) / math.sqrt(len(MmuDivergence[i])))
for i in range(len(PtrDivergence)):
    PtrMeanExpDiv.append(np.mean(PtrDivergence[i]))
    PtrSEMExpDiv.append(np.std(PtrDivergence[i]) / math.sqrt(len(PtrDivergence[i])))
print('computed means and SEM')




# perform statistical tests between gene categories
# save P values to file
newfile = open('MouseChimpExpDivergDistancePVals.txt', 'w')
newfile.write('\t'.join(['Species', 'Genes1', 'Genes2', 'index1', 'index2', 'P']) + '\n')        
# loop over inner list, compare gene categories
for i in range(0, len(MmuDivergence) -1):
    for j in range(i+1, len(MmuDivergence)):
        P = PermutationResampling(MmuDivergence[i], MmuDivergence[j], 1000, statistic = np.mean)
        print('Mouse', GeneCats[i], GeneCats[j], i, j, P)
        newfile.write('\t'.join(['Mouse', GeneCats[i], GeneCats[j], str(i), str(j), str(P)]) + '\n')        
# loop over inner list, compare gene categories
for i in range(0, len(PtrDivergence) -1):
    for j in range(i+1, len(PtrDivergence)):
        P = PermutationResampling(PtrDivergence[i], PtrDivergence[j], 1000, statistic = np.mean)
        print('Chimp', GeneCats[i], GeneCats[j], i, j, P)
        newfile.write('\t'.join(['Chimp', GeneCats[i], GeneCats[j], str(i), str(j), str(P)]) + '\n')        
newfile.close()


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XTickLabel, Species):
    '''
    return an ax object part of figure
    '''

    # add a plot to figure (N row, N column, plot N)
    ax = fig.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['#f03b20', '#43a2ca', '#fee391', '#74c476', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
    # plot nucleotide divergence
    ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], Data[0], 0.2, yerr = Data[1], color = colorscheme,
           edgecolor = 'black', linewidth = 0.7,
           error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write axis label
    ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
    ax.set_xlabel(Species, color = 'black', size = 7, ha = 'center', **FigFont)
    # add ticks and lebels
    plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], XTickLabel, size = 7, color = 'black', ha = 'center', **FigFont)
    
    # add a range for the Y axis
    plt.yticks(np.arange(0, 0.8, 0.1))
    # add a range for the Y and X axes
    #plt.ylim([0, 0.8])
    plt.xlim([0, 2.45])
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
    return ax  

# create figure
fig = plt.figure(1, figsize = (5, 2))
# set colors
colorscheme = ['#f03b20', '#43a2ca', '#fee391', '#74c476', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
# plot data
ax1 = CreateAx(2, 1, 1, fig, [MmuMeanExpDiv, MmuSEMExpDiv], GeneCats, 'Mouse')
ax2 = CreateAx(2, 1, 2, fig, [PtrMeanExpDiv, PtrSEMExpDiv], GeneCats, 'Chimp')


# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
MmuDiff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'A']
PtrDiff = ['A', 'AB', 'BC', 'BC', 'BC', 'ABD', 'ABE', 'BF']

ypos = [0.65] * 8
xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
for i in range(len(MmuDiff)):
    ax1.text(xpos[i], ypos[i], MmuDiff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 6)
for i in range(len(PtrDiff)):
    ax2.text(xpos[i], ypos[i], PtrDiff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 6)

# add subplot labels
ax1.text(-0.5, 0.8, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)   
ax1.text(2.8, 0.8, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)

# make sure subplots do not overlap
plt.tight_layout() 

# save figure
outputfile = 'MouseChimpExpDivgDist'
for extension in ['.pdf', '.eps', '.png']:
    fig.savefig(outputfile + extension, bbox_inches = 'tight')
