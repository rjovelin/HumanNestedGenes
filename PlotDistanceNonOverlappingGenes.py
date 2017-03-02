# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:30:11 2017

@author: RJovelin
"""


# use this script to plot the distances between non-overlapping orthologs of human overlapping genes

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
jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json']

# make a list of dictionaries
AllOverlap = []
# loop over files
for i in range(len(jsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(jsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    AllOverlap.append(overlapping)

# get GFF file
GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3']

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
HumanOrdered, ChimpOrdered = AllOrdered[0], AllOrdered[1]
HumanCoord, ChimpCoord = AllCoordinates[0], AllCoordinates[1]

# get 1:1 orthologs between human and chimp
Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')

# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# get the gene pairs
HumanPairs = AllPairs[:2]
ChimpPairs = AllPairs[2:]

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos or pair[1] not in Orthos:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
# remove order for chimp gene pairs
for i in range(len(ChimpPairs)):
    for j in range(len(ChimpPairs[i])):
        ChimpPairs[i][j] = set(ChimpPairs[i][j])



# 1) plot the proportions of orthologous gene pairs that are adjacent, separated
# and on different chromosomes
# 2) plot a histogram of the intergenic distance between non-nested adjacent
# chimp orthologs of human nested genes

SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG = 0, 0, 0, 0, 0

# get the intergenic distance between non-nested adjacent orthologs of human nested genes
GeneDist = []  
 
# loop over human nested gene pairs
for i in range(len(HumanPairs[1])):
    # get the orthologs of the human genes
    ortho1, ortho2 = Orthos[HumanPairs[1][i][0]], Orthos[HumanPairs[1][i][1]]
    # check if human gene pairs is nested in chimp
    if set([ortho1, ortho2]) not in ChimpPairs[1]:
        # get chromosomes of chimp orthologs
        chromo1, chromo2 = ChimpCoord[ortho1][0], ChimpCoord[ortho2][0]        
        # check if chromosomes are different
        if chromo1 != chromo2:
            DiffLG += 1
        else:
            # get start positions of chimp orthologs
            S1, S2 = ChimpCoord[ortho1][1], ChimpCoord[ortho2][1]
            # get end positions of chimp orthologs
            E1, E2 = ChimpCoord[ortho1][2], ChimpCoord[ortho2][2]            
            # get indices of chimp orthologs in the ordred gene list
            P1, P2 = ChimpOrdered[chromo1].index(ortho1), ChimpOrdered[chromo2].index(ortho2)
            D = 'na'            
            # check if orthologs are overlapping in chimp
            if set([ortho1, ortho2]) in ChimpPairs[0]:
                # chimp genes are overlapping
                # check if they are adjcent
                if S1 < S2:
                    assert P1 < P2
                    if P2 != P1 + 1:
                        # overlapping non adjacent
                        SepOvlp += 1
                    elif P2 == P1 + 1:
                        # overlapping adjacent
                        AdjOvlp += 1
                        # get the distance between genes
                        D = S2 - E1
                elif S1 > S2:
                    assert P2 < P1
                    if P1 != P2 + 1:
                        # overlapping non adjacent
                        SepOvlp += 1
                    elif P1 == P2 + 1:
                        # overlapping adjacent
                        AdjOvlp += 1
                        # get the distance between genes              
                        D = S1 - E2
            elif set([ortho1, ortho2]) not in ChimpPairs[0]:
                # chimp genes are not overlapping
                # check if they are adjacent
                if S1 < S2:
                    assert P1 < P2
                    if P2 != P1 + 1:
                        # non-overlapping non-adjacent
                        SepNonOvlp += 1
                    elif P2 == P1 + 1:
                        # non-overlapping adjacent
                        AdjNonOvlp += 1
                        # get the distance between genes
                        D = S2 - E1
                elif S1 > S2:
                    assert P2 < P1
                    if P1 != P2 + 1:
                        # non-overlapping non adjacent
                        SepNonOvlp += 1
                    elif P1 == P2 + 1:
                        # non-overlapping adjacent
                        AdjNonOvlp += 1
                        # get the distance between genes              
                        D = S1 - E2
            if D != 'na':
                GeneDist.append(D)  
                
            
print(SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG)

negatif, positif, zero = [], [], []
for i in GeneDist:
    if i > 0:
        positif.append(i)
    elif i < 0:
        negatif.append(i)
    elif i == 0:
        zero.append(i)
negatif.sort()
positif.sort()
print(len(negatif), len(positif), len(zero))
print(min(negatif), max(negatif))
print(min(positif), max(positif))


print(positif)
print('\n\n')
print(negatif)


# 1) plot the proportions of orthologous gene pairs that are adjacent, separated and on different chromosomes






# 2) plot a histogram of the intergenic distance between non-nested adjacent chimp orthologs of human nested genes


###################



## create figure
#fig = plt.figure(1, figsize = (2.5, 4.5))
#
## plot distribution of intron number
#YMax1 = max(max([WithIntronCount.count(i) for i in WithIntronCount]),  max([IntronlessCount.count(i) for i in IntronlessCount])) + 10
#ax1 = CreateAx(1, 3, 1, fig, [WithIntronCount, IntronlessCount], 'N external genes', 'Number of introns', np.arange(0, YMax1 + 10, 10), YMax1, np.arange(0, 80, 10), ['orange', 'blue'], 'histo')
## plot distribution of intron position
#YMax2 = max(max([WithIntronPos.count(i) for i in WithIntronPos]), max([IntronlessPos.count(i) for i in IntronlessPos])) + 10
#ax2 = CreateAx(1, 3, 2, fig, [WithIntronPos, IntronlessPos], 'N external genes', 'Intron position', np.arange(0, YMax2 + 20, 20), YMax2, np.arange(0, 60, 10), ['orange', 'blue'], 'histo')
## plot distribution of intron position / intron number
#ax3 = CreateAx(1, 3, 3, fig, [WithIntronPosNorm, PWithPosNorm, IntronlessPosNorm, PNonePosNorm] , 'Probability', 'Intron position / intron number', np.arange(0, 1.1, 0.1), 1, np.arange(0, 1.1, 0.1), ['orange', 'blue'], 'cdf')
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
#    # create subplot in figure
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#    # plot data    
#    if GraphType == 'histo':
#        ax.hist(Data[0], bins = np.arange(0, max([max(Data[0]), max(Data[1])]) + 1, 1), linewidth = 0.7, histtype='step', fill = True, facecolor = Colors[0], edgecolor = Colors[0], alpha = 0.5, stacked = False)    
#        ax.hist(Data[1], bins = np.arange(0, max([max(Data[0]), max(Data[1])]) + 1, 1), linewidth = 0.7, histtype='step', fill = True, facecolor = Colors[1], edgecolor = Colors[1], alpha = 0.5, stacked = False)
#    elif GraphType == 'cdf':
#        ax.step(Data[0], Data[1], linewidth = 1, linestyle = '-', color = Colors[0], alpha = 0.5)
#        ax.step(Data[2], Data[3], linewidth = 1, linestyle = '-', color = Colors[1], alpha = 0.5)
#        
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write label axis
#    ax.set_ylabel(YLabel, color = 'black',  size = 6.5, ha = 'center', **FigFont)
#    ax.set_xlabel(XLabel, color = 'black',  size = 6.5, ha = 'center', **FigFont)
#    # set a limit to y axis
#    plt.ylim([0, YMax])
#    # add a range to the axis
#    plt.xticks(XRange, size = 6.5, color = 'black', ha = 'center', **FigFont)
#    plt.yticks(YRange, size = 6.5, color = 'black', ha = 'right', **FigFont)        
#    
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)    
#    ax.spines["left"].set_visible(True)  
#    # edit tick paramters
#    plt.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
#                    left = 'on', labelbottom='on', colors = 'black', labelsize = 7, direction = 'out')  
#    # add ticks on the x axis
#    #plt.xticks(TickPos, Ticklabel)    
#    # Set the tick labels font name
#    for label in ax.get_yticklabels():
#        label.set_fontname('Arial')   
#    # create a margin around the x axis
#    plt.margins(0.05)
#    return ax      
#
#
## make sure subplots do not overlap
#plt.tight_layout()
#
### annotate graphs with legends
## create patches 
#Title = mpatches.Patch(facecolor = 'none', edgecolor = 'none', linewidth = 0, label= 'Internal genes:', alpha = 0)
#Yes = mpatches.Patch(facecolor = 'orange', edgecolor = 'none', linewidth = 0.5, label= 'with introns', alpha = 0.5)
#No = mpatches.Patch(facecolor = 'blue', edgecolor = 'black', linewidth = 0.5, label= 'intronless', alpha = 0.5)
#ax1.legend(handles = [Title, Yes, No], bbox_to_anchor=(0.2, 0.4), loc = 3, fontsize = 5, frameon = False, ncol = 1)
#ax2.legend(handles = [Title, Yes, No], bbox_to_anchor=(0.2, 0.4), loc = 3, fontsize = 5, frameon = False, ncol = 1)
## create lines
#orange_line = mlines.Line2D([], [], color='orange', marker='', markersize=15, label='with introns', alpha = 0.5)
#blue_line = mlines.Line2D([], [], color='blue', marker='', markersize=15, label='intronless', alpha = 0.5)
#ax3.legend(handles = [Title, orange_line, blue_line], bbox_to_anchor=(0.05, 0.5), loc = 3, fontsize = 5, frameon = False, ncol = 1)
#
## annotate graphs with p values
#for i in range(len(PVals)):
#    if PVals[i] != '':
#        if i == 0:
#            ax1.text(25, 20, PVals[i], color = 'black', size = 5, fontname = 'Arial')
#        elif i == 1:
#            ax2.text(15, 60, PVals[i], color = 'black', size = 5, fontname = 'Arial')
#        elif i == 2:
#            ax3.text(0.4, 0.3, PVals[i], color = 'black', size = 5, fontname = 'Arial')
#
#fig.savefig('truc.pdf', bbox_inches = 'tight')
#
#
############################
#
#
#
