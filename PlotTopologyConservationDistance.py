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

# get 1:1 orthologs between human and mouse
Orthos = MatchOrthologs('HumanMouseOrthologs.txt')

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
        # check distance between orthologs in mouse
        for pair in HumanAdjacentgenePairs[GeneType]:
            pair = list(pair)
            # make a list of orthologs in mouse
            L = []
            for ortho1 in Orthos[pair[0]]:
                for ortho2 in Orthos[pair[1]]:
                    L.append([ortho1, ortho2])
            # compute divergence between genes in mouse
            for genes in L:
                gene1, gene2 = genes[0], genes[1]
                # check if genes are valid mouse genes
                if gene1 in MouseCoord and gene2 in MouseCoord:
                    # check if genes are on the same chromo and that genes are different
                    if (MouseCoord[gene1][0] == MouseCoord[gene2][0]) and (gene1 != gene2):
                        # get indices of gene1 and gene2
                        I1, I2 = MouseOrdered[MouseCoord[gene1][0]].index(gene1), MouseOrdered[MouseCoord[gene2][0]].index(gene2)
                        assert I1 != I2                        
                        if I1 < I2:
                            start1, end1 = MouseCoord[gene1][1], MouseCoord[gene1][2]
                            start2, end2 = MouseCoord[gene2][1], MouseCoord[gene2][2]
                        elif I1 > I2:
                            start1, end1 = MouseCoord[gene2][1], MouseCoord[gene2][2]
                            start2, end2 = MouseCoord[gene1][1], MouseCoord[gene1][2]
                        D = start2 - end1
                        if D >= 0 and D < 1000 and GeneType == 'Proximal':
                            ConservedPairs[GeneType] += 1
                        elif D >= 1000 and D < 10000 and GeneType == 'Moderate':
                            ConservedPairs[GeneType] += 1
                        elif D >= 10000 and D < 50000 and GeneType == 'Intermediate':
                            ConservedPairs[GeneType] += 1             
                        elif D >= 50000 and GeneType == 'Distant':
                            ConservedPairs[GeneType] += 1
                        # record only 1 pair of orthologs if multiple co-ortholiogs exit, exit loop 
                        break
                
                
for i in ConservedPairs:
    print(i, ConservedPairs[i] / len(HumanAdjacentgenePairs[i]))                
                
        
        
##### see below for checking conservation
## loop over human overlapping gene classes
#for i in range(len(AllPairs)):
#    # create a list with pair counts [both gene conserved, 1 conserved, non conserved]
#    counts = [0, 0, 0]
#    # make a set of orthologous overlapping gene
#    orthosovlp = set()
#    for pair in PairsOrthos[i]:
#        pair = list(pair)
#        for item in pair:
#            orthosovlp.add(item)
#    # loop over human gene pairs
#    for pair in AllPairs[i][0]:
#        if set(pair) in PairsOrthos[i]:
#            counts[0] += 1
#        else:
#            if pair[0] in orthosovlp or pair[1] in orthosovlp:
#                counts[1] += 1
#            elif pair[0] not in orthosovlp and pair[1] not in orthosovlp:
#                counts[2] += 1
#    # divide by number of gene pairs to get proportions    
#    assert sum(counts) == len(AllPairs[i][0])
#    for j in range(len(counts)):
#        counts[j] = counts[j] / len(AllPairs[i][0])
#    PairCounts.append(counts)
#    assert sum(counts) == 1















#
#
## 1) plot proportions of gene pairs with varying distance conserved in chimp and mouse
## 2) plot proportions of overlapping gene pairs with conserved topology in chimp and mouse
## 3) plot proportions of overlapping gene pairs that are overlapping in chimp and mouse
## 4) plot differences between conservation of human overalapping in chimop and mouse
##    and overlapping genes in chimp and mouse conserved in human human
#
## create a dict with lists of gene pairs
#AllGenePairs = {}
#AllGenePairs['hsa_ptr'] = copy.deepcopy(HumanPairs)
#AllGenePairs['hsa_mmu'] = copy.deepcopy(HsaPairs)
#AllGenePairs['ptr'] = copy.deepcopy(ChimpPairs)
#AllGenePairs['mmu'] = copy.deepcopy(MousePairs)
#
## replace human genes by their orthologs
#species = ['hsa_ptr', 'hsa_mmu']
#for i in range(len(species)):
#    for j in range(len(AllGenePairs[species[i]])):
#        for k in range(len(AllGenePairs[species[i]][j])):
#            AllGenePairs[species[i]][j][k][0] = Orthos[i][AllGenePairs[species[i]][j][k][0]]
#            AllGenePairs[species[i]][j][k][1] = Orthos[i][AllGenePairs[species[i]][j][k][1]]
#
## sort gene pairs
#for sp in AllGenePairs:
#    for i in range(len(AllGenePairs[sp])):
#        for j in range(len(AllGenePairs[sp][i])):
#            AllGenePairs[sp][i][j].sort()
#
## replace gene lists with strings
#for sp in AllGenePairs:
#    for i in range(len(AllGenePairs[sp])):
#        for j in range(len(AllGenePairs[sp][i])):
#            AllGenePairs[sp][i][j] = ':'.join(AllGenePairs[sp][i][j])
#
## do qc
#for sp in AllGenePairs:
#    for i in range(1, len(AllGenePairs[sp])):
#        for pair in AllGenePairs[sp][i]:
#            assert pair in AllGenePairs[sp][0]
#
## make lists of adjacent gene pairs with varying distance in human conserved in chimp and mouse
## [[gene1, gene2], ....[gene n, gene n+1]]
#for i in range(len(species)):
#    HsaPairsDist = []    
#    for j in range(7):
#        HsaPairsDist.append([])
#    # loop over chromosomes
#    for chromo in HumanOrdered:
#        # loop over the list of ordered genes
#        for k in range(len(HumanOrdered[chromo]) - 1):
#            # get the end position of gene 1
#            EndGene1 = HumanCoord[HumanOrdered[chromo][k]][2]                
#            # get the start position of adjacent gene 2
#            StartGene2 = HumanCoord[HumanOrdered[chromo][k+1]][1]
#            # computance distance between genes
#            D = StartGene2 - EndGene1
#            # assign infinity value to k
#            m = float('inf')
#            if D < 0:
#                m = 0
#            elif D >= 0 and D < 1000:
#                m = 1
#            elif D >= 1000 and D < 10000:
#                m = 2                
#            elif D >= 10000 and D < 50000:
#                m = 3
#            elif D >= 50000 and D < 100000:
#                m = 4
#            elif D >= 100000 and D < 150000:
#                m = 5
#            elif D >= 150000:
#                m = 6
#            if HumanOrdered[chromo][k] in Orthos[i] and HumanOrdered[chromo][k+1] in Orthos[i] and m in range(7):
#                # add the human orthologs
#                pair = [Orthos[i][HumanOrdered[chromo][k]], Orthos[i][HumanOrdered[chromo][k+1]]]
#                # sort gene pair
#                pair.sort()
#                # concert list to string
#                pair = ':'.join(pair)
#                HsaPairsDist[m].append(pair)
#    # add the pairs of non-overlapping genes to the lists of gene pairs
#    AllGenePairs[species[i]].extend(HsaPairsDist)
#
#
## make list of adjacent gene pairs regardless of distance in chimp and mouse [gene1:gene2, ....gene_n:gene_n+1]
#sp2 = ['ptr', 'mmu']
#Sp2Coord = [ChimpCoord, MouseCoord]
#Sp2Ordered = [ChimpOrdered, MouseOrdered]
#for i in range(len(sp2)):
#    Sp2PairsDist = []
#    # loop over chromosomes
#    for chromo in Sp2Ordered[i]:
#        # loop over the list of ordered genes
#        for j in range(len(Sp2Ordered[i][chromo]) - 1):
#            # get the end position of gene 1
#            EndGene1 = Sp2Coord[i][Sp2Ordered[i][chromo][j]][2]                
#            # get the start position of adjacent gene 2
#            StartGene2 = Sp2Coord[i][Sp2Ordered[i][chromo][j+1]][1]
#            # check that genes are not overlapping
#            D = StartGene2 - EndGene1
#            if D >= 0:
#                # get gene pair
#                pair = [Sp2Ordered[i][chromo][j], Sp2Ordered[i][chromo][j+1]]
#                # sort pair
#                pair.sort()
#                # convert list to string
#                pair = ':'.join(pair)
#                Sp2PairsDist.append(pair)
#    AllGenePairs[sp2[i]].append(Sp2PairsDist)
#
#
## convert lists to numpy arrays
#for sp in AllGenePairs:
#    for i in range(len(AllGenePairs[sp])):
#        AllGenePairs[sp][i] = np.array(AllGenePairs[sp][i])
#
#
## 1) plot proportions of gene pairs with varying distance conserved in chimp and mouse
## get the proportions of non-overlapping adjacent gene pairs with conserved topology
#NonOverlap = []
#for i in range(5, len(AllGenePairs['hsa_ptr'])):
#    # count the number of adjacent gene pairs gene pairs that are adjacent in chimp
#    total = sum(np.in1d(AllGenePairs['hsa_ptr'][i], AllGenePairs['ptr'][5], invert = False))
#    NonOverlap.append((total / len(AllGenePairs['hsa_ptr'][i])) * 100)
#    # count the number of adjacent gene pairs gene pairs that are adjacent in mouse
#    total = sum(np.in1d(AllGenePairs['hsa_mmu'][i], AllGenePairs['mmu'][5], invert = False))
#    NonOverlap.append((total / len(AllGenePairs['hsa_mmu'][i])) * 100)    
#
## 2) plot proportions of overlapping gene pairs with conserved topology in chimp and mouse
## get the proportions of overlapping gene pairs with conserved topology in each overlapping category
#OverlapCat = []
#for i in range(5):
#    # count the number of overlapping gene pairs that are also overlapping in the same category in chimp
#    total = sum(np.in1d(AllGenePairs['hsa_ptr'][i], AllGenePairs['ptr'][i], invert = False))
#    OverlapCat.append((total / len(AllGenePairs['hsa_ptr'][i])) * 100)
#    # count the number of overlapping gene pairs that are also overlapping in the same category in mouse
#    total = sum(np.in1d(AllGenePairs['hsa_mmu'][i], AllGenePairs['mmu'][i], invert = False))
#    OverlapCat.append((total / len(AllGenePairs['hsa_mmu'][i])) * 100)
#
## 3) plot proportions of overlapping gene pairs that are overlapping in chimp and mouse
## get the proportions of overlapping gene pairs that are also overlapping (regadless of category)
#OverlapAll = []
#for i in range(5):
#    # count the number of overlapping gene pairs that are also overlapping in chimp
#    total = sum(np.in1d(AllGenePairs['hsa_ptr'][i], AllGenePairs['ptr'][0], invert = False))
#    OverlapAll.append((total / len(AllGenePairs['hsa_ptr'][i])) * 100)
#    # count the number of overlapping gene pairs that are also overlapping in mouse
#    total = sum(np.in1d(AllGenePairs['hsa_mmu'][i], AllGenePairs['mmu'][0], invert = False))
#    OverlapAll.append((total / len(AllGenePairs['hsa_mmu'][i])) * 100)
#    
#
## 4) plot differences between conservation of human overalapping in chimop and mouse
##    and overlapping genes in chimp and mouse conserved in human human
#
## make a reverse dictionary of orthologs 
#ChimpOrthos = {}
#for gene in HsaPtrOrthos:
#    assert HsaPtrOrthos[gene] not in ChimpOrthos
#    ChimpOrthos[HsaPtrOrthos[gene]] = gene
#MouseOrthos = {}
#for gene in HsaMmuOrthos:
#    assert HsaMmuOrthos[gene] not in MouseOrthos
#    MouseOrthos[HsaMmuOrthos[gene]] = gene
#
## remove gene pairs lacking orthologs
#for i in range(len(ChimpPairs)):
#    to_remove = []
#    for pair in ChimpPairs[i]:
#        if pair[0] not in ChimpOrthos or pair[1] not in ChimpOrthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        ChimpPairs[i].remove(pair)
#for i in range(len(MousePairs)):
#    to_remove = []
#    for pair in MousePairs[i]:
#        if pair[0] not in MouseOrthos or pair[1] not in MouseOrthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        MousePairs[i].remove(pair)
#
#
## add dictionaries of orthologs to list
#Orthos.append(ChimpOrthos)
#Orthos.append(MouseOrthos)
#
## make a list of gene pair lists
#GenePairs = [HumanPairs, HsaPairs, ChimpPairs, MousePairs]
## make a list of sets of gene pairs
#GeneSets = []
#for i in range(len(GenePairs)):
#    spsets = []
#    for j in range(len(GenePairs[i])):
#        spsets.append([set(k) for k in GenePairs[i][j]])
#    GeneSets.append(spsets)
#
## make a list of counts of conserved and non-conserved human overlapping genes in chimp and mouse
#Conserved = []
#for i in range(len(GenePairs)):
#    if i < 2:
#        j = i + 2
#    else:
#        j = i - 2
#    conservation = []
#    for k in range(len(GenePairs[i])):
#        conserved, divergent = 0, 0
#        for pair in GenePairs[i][k]:
#            if set([Orthos[i][pair[0]], Orthos[i][pair[1]]]) in GeneSets[j][k]:
#                conserved += 1
#            else:
#                divergent += 1
#        conservation.append([conserved, divergent])
#    Conserved.append(conservation) 
#
#
## create lists with proportions of overlapping genes with conserved topologies
## [[prop human overlapping genes conserved, prop chimp overlapping genes conserved in human, etc]]
#Proportions = []
#for i in range(2):
#    prop = []
#    for j in range(len(Conserved[i])):
#        man = Conserved[i][j][0] / sum(Conserved[i][j])
#        sp2 = Conserved[i+2][j][0] / sum(Conserved[i+2][j])
#        diff = (sp2 - man) * 100
#        prop.append(diff)
#    Proportions.append(prop)
#
#  
## create a single list with differences between human and chimp and between
## human and mouse for each overlapping gene category
#Differences = []
#for i in range(len(Proportions[0])):
#    Differences.append(Proportions[0][i])
#    Differences.append(Proportions[1][i])
#
#
#
#
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
