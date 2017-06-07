# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:46:17 2017

@author: RJovelin
"""

# use this script to plot the % of orthologous gene pairs with same topology between human and mouse

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


# load dictionaries of overlapping genes
JsonFiles = ['Overlapping', 'Nested', 'PiggyBack', 'Convergent', 'Divergent']

# make a list of dictionaries
HsaAllOverlap, MmuAllOverlap = [], []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open('Human' + JsonFiles[i] + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    HsaAllOverlap.append(overlapping)
    json_data = open('Mouse' + JsonFiles[i] + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    MmuAllOverlap.append(overlapping)

# get GFF file
GFF = ['Homo_sapiens.GRCh38.88.gff3', 'Mus_musculus.GRCm38.88.gff3']

# make a list of gene coordinates       
AllCoordinates, AllOrdered = [], []
# loop over GFF files
for i in range(len(GFF)):
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    AllCoordinates.append(GeneCoord)
    AllOrdered.append(OrderedGenes)
HumanOrdered, MouseOrdered = AllOrdered[0], AllOrdered[1]
HumanCoord, MouseCoord = AllCoordinates[0], AllCoordinates[1]

# get orthologs between human and mouse
Orthos = MatchOrthologs('HumanMouseOrthologs.txt')
# reverse dictionary 
MouseOrthologs = {}
for gene in Orthos:
    for ortho in Orthos[gene]:
        if ortho not in MouseOrthologs:
            MouseOrthologs[ortho] = [gene]
        else:
            MouseOrthologs[ortho].append(gene)


# record adjacent gene pairs in human (remove order)
HumanAdjacentgenePairs = {}
for i in ['Proximal', 'Moderate', 'Intermediate', 'Distant']:
    HumanAdjacentgenePairs[i] = []
# loop over chromo in human
for chromo in HumanOrdered:
    # loop over the list of ordered genes
    for i in range(len(HumanOrdered[chromo]) -1):
        gene1, gene2 = HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]
        # only consider gene pairs with orthologs
        if gene1 in Orthos and gene2 in Orthos:
            # get the end position of gene 1
            EndGene1 = HumanCoord[gene1][2]
            # get the start position of adjacent gene 2
            StartGene2 = HumanCoord[gene2][1]
            # computance distance between genes
            D = StartGene2 - EndGene1
            if D >= 0 and D < 1000:
                HumanAdjacentgenePairs['Proximal'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))
            elif D >= 1000 and D < 10000:
                HumanAdjacentgenePairs['Moderate'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))
            elif D >= 10000 and D < 50000:
                HumanAdjacentgenePairs['Intermediate'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))            
            elif D >= 50000:
                HumanAdjacentgenePairs['Distant'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))

# make pairs of overlapping genes
for i in range(1, len(HsaAllOverlap)):
    pairs = GetHostNestedPairs(HsaAllOverlap[i])
    # remove pairs if any gene is lacking an ortholog
    to_remove = [L for L in pairs if L[0] not in Orthos or L[1] not in Orthos]    
    for L in to_remove:
        pairs.remove(L)
    HumanAdjacentgenePairs[JsonFiles[i]] = pairs    
    
# count gene pairs conserved in mouse    
ConservedPairs = {}
for GeneType in HumanAdjacentgenePairs:
    # initialize counters
    ConservedPairs[GeneType] = 0
    if GeneType in ['Proximal', 'Moderate', 'Intermediate', 'Distant']:
        # check if mouse orthologs are adjacent
        for pair in HumanAdjacentgenePairs[GeneType]:
            pair = list(pair)
            # make a list of orthologs in mouse
            L = []
            for ortho1 in Orthos[pair[0]]:
                for ortho2 in Orthos[pair[1]]:
                    L.append([ortho1, ortho2])
            for genes in L:
                gene1, gene2 = genes[0], genes[1]
                # check if genes are valid mouse genes
                if gene1 in MouseCoord and gene2 in MouseCoord:
                    # check if genes are on the same chromo and that genes are different
                    if (MouseCoord[gene1][0] == MouseCoord[gene2][0]) and (gene1 != gene2):
                        # get indices of gene1 and gene2
                        I1, I2 = MouseOrdered[MouseCoord[gene1][0]].index(gene1), MouseOrdered[MouseCoord[gene2][0]].index(gene2)
                        assert I1 != I2                        
                        # check if genes are adjacent
                        if I1 == I2 + 1 or I2 == I1 + 1:
                            ConservedPairs[GeneType] += 1
                        # record only 1 pair of orthologs if multiple co-ortholiogs exit, exit loop 
                        break
                
                
# make pairs of human orthologs for each mouse gene pair
PairsOrthos = []
for i in range(1, len(MmuAllOverlap)):
    # create a list to store the gene pairs
    GenePairs = []
    Pairs = GetHostNestedPairs(MmuAllOverlap[i])
    # loop over mouse gene pairs
    for pair in Pairs:
        # check that both genes have coordinates
        assert pair[0] in MouseCoord and pair[1] in MouseCoord
        # count only pairs in which both genes have orthologs in human
        if pair[0] in MouseOrthologs and pair[1] in MouseOrthologs:
            for ortho1 in MouseOrthologs[pair[0]]:
                for ortho2 in MouseOrthologs[pair[1]]:
                    # remove order
                    GenePairs.append(set([ortho1, ortho2]))
    PairsOrthos.append(GenePairs)
# loop over human overlapping gene classes    
for i in range(1, len(JsonFiles)):
    for pair in HumanAdjacentgenePairs[JsonFiles[i]]:
        if set(pair) in PairsOrthos[i-1]:
            ConservedPairs[JsonFiles[i]] += 1

# compute proportions
for i in ConservedPairs:
    print(i, ConservedPairs[i] / len(HumanAdjacentgenePairs[i]))                















      
        


## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, figure, Data, YLabel):
#    '''
#    Returns a ax instance in figure
#    '''    
#
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#
#    if Position == 1:
#        BarPos = [0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4, 2.7, 2.9, 3.2, 3.4]
#        XTickpos = [0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4]
#        Alignment = 'right'
#        XTicklabels = ['< 0', '0-1', '1-10', '10-50', '50-100', '100-150', '> 150']
#        # draw x axis line
#        ax.plot([0, 3.8], [0, 0], lw = 0.7, color = 'black')
#    else:
#        BarPos = [0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4]
#        XTickpos = [0.4, 0.9, 1.4, 1.9, 2.4]
#        Alignment = 'center'
#        XTicklabels = ['all', 'nst', 'pbk', 'conv', 'div']
#        # draw x axis line
#        ax.plot([0, 2.8], [0, 0], lw = 0.7, color = 'black')
#    
#    if Position == 4:
#        YTicksRange = np.arange(-20, 120, 20)
#        YMin, YMax = -20, 100
#    else:
#        YTicksRange = np.arange(0, 120, 20)
#        YMin, YMax = 0, 100
#        
#    Colors = ['black', 'lightgrey'] * 5
#
#    # plot data
#    ax.bar(BarPos, Data, width = 0.2, color = Colors, edgecolor = 'black', linewidth = 0.7)
#    
##    # draw x axis line
##    if Position == 4:
##        ax.plot([0, 2.8], [0, 0], lw = 0.7, color = 'black')
#    
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write y axis label
#    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
#    # add ticks and lebels
#    if Position == 1:
#        Rotation = 30
#    else:
#        Rotation = 0
#    plt.xticks(XTickpos, XTicklabels, rotation = Rotation, size = 7, color = 'black', ha = Alignment, **FigFont)
#    # edit y axis ticks
#    plt.yticks(YTicksRange)
#    # add a range for the Y and X axes
#    plt.ylim([YMin, YMax])    
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(False)    
#    ax.spines["right"].set_visible(False)
#    ax.spines["left"].set_visible(True)  
#    # make sure the y axis crosses the x axis at 0
#    ax.spines['left'].set_position('zero')
#
#    # edit tick parameters
#    if Position == 4:
#        BottomTicks = 'off'
#    else:
#        BottomTicks = 'on'
#    plt.tick_params(axis='both', which='both', bottom=BottomTicks, top='off',
#                    right = 'off', left = 'on', labelbottom='on',
#                    colors = 'black', labelsize = 7, direction = 'out') 
#    # Set the tick labels font name
#    for label in ax.get_yticklabels():
#        label.set_fontname('Arial')   
#    
#    ## add margins
#    #plt.margins(0.1)
#    return ax
#
#
#
## create figure
#fig = plt.figure(1, figsize = (5, 4))
#
#YLabels = ['% with orthologous\nadjacent gene pairs',
#           '% overlapping gene pairs\nwith conserved topology',
#           '% with orthologous\noverlapping gene pairs ',
#           '% excess of orthologous gene\npairs with conserved topology']
#
#
## 1) plot proportions of gene pairs with varying distance conserved in chimp and mouse
#ax1 = CreateAx(2, 2, 1, fig, NonOverlap, YLabels[0]) 
## 2) plot proportions of overlapping gene pairs with conserved topology in chimp and mouse
#ax2 = CreateAx(2, 2, 2, fig, OverlapCat, YLabels[1]) 
## 3) plot proportions of overlapping gene pairs that are overlapping in chimp and mouse
#ax3 = CreateAx(2, 2, 3, fig, OverlapAll, YLabels[2]) 
## 4) plot differences between conservation of human overalapping in chimop and mouse
##    and overlapping genes in chimp and mouse conserved in human human
#ax4 = CreateAx(2, 2, 4, fig, Differences, YLabels[3])
#
## add subplot labels
#ax1.text(-0.8, 113, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#ax2.text(-0.8, 113, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#ax3.text(-0.8, 113, 'C', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#ax4.text(-0.8, 113, 'D', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#
## add legend
#mouse = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'mouse')
#chimp = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'chimp')
#ax1.legend(handles = [chimp, mouse], loc = (0.17, 1), fontsize = 7, frameon = False, ncol = 2)
#
## make sure subplots do not overlap
#plt.tight_layout()
#
## save figure to file
#fig.savefig('ConservationTopology.pdf', bbox_inches = 'tight')
#fig.savefig('ConservationTopology.eps', bbox_inches = 'tight')
