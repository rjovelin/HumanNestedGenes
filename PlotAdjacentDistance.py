# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:05:32 2017

@author: Richard
"""


# use this script to plot conservation of overlapping genes as a function of the physical distance

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



# 1) plot the % of orthologous gene pairs for which orthologs are adjacent in mouse

# load dictionaries of overlapping genes
JsonFiles = ['Nested', 'PiggyBack', 'Convergent', 'Divergent']
# make a list of dictionaries
HsaAllOverlap, MmuAllOverlap = [], []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open('Human' + JsonFiles[i] + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    HsaAllOverlap.append(overlapping)
    data = open('Mouse' + JsonFiles[i] + 'Genes.json')
    overlapping = json.load(data)
    data.close()
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
for i in range(len(HsaAllOverlap)):
    pairs = GetHostNestedPairs(HsaAllOverlap[i])
    # remove pairs if any gene is lacking an ortholog
    to_remove = [L for L in pairs if L[0] not in Orthos or L[1] not in Orthos]    
    for L in to_remove:
        pairs.remove(L)
    HumanAdjacentgenePairs[JsonFiles[i]] = pairs    
    
    
# make pairs of human orthologs for each mouse gene pair
PairsOrthos = []
for i in range(len(MmuAllOverlap)):
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


    
# count gene pairs adjacent in mouse    
ConservedPairs = {}
for GeneType in ['Proximal', 'Moderate', 'Intermediate', 'Distant']:
    # initialize counters
    ConservedPairs[GeneType] = 0
    # check if mouse orthologs are adjacent
    for pair in HumanAdjacentgenePairs[GeneType]:
        pair = list(pair)
        # make a list of orthologs in mouse
        L = []
        for ortho1 in Orthos[pair[0]]:
            for ortho2 in Orthos[pair[1]]:
                # check that both genes have coordinates in mouse
                if ortho1 in MouseCoord and ortho2 in MouseCoord:
                    # check that genes are different and on the same chromo 
                    if ortho1 != ortho2 and MouseCoord[ortho1][0] == MouseCoord[ortho2][0]:
                        L.append([ortho1, ortho2])
        # check if a pair of orthologs exist
        if len(L) != 0:
            # set up boolean to be changed if an adjacent pair is found
            FoundAdjacent = False
            # loop over the mouse gene pairs
            for i in range(len(L)):
                gene1, gene2 = L[i][0], L[i][1]
                # get indices of gene1 and gene2
                I1, I2 = MouseOrdered[MouseCoord[gene1][0]].index(gene1), MouseOrdered[MouseCoord[gene2][0]].index(gene2)
                # check if genes are adjacent
                if I1 == I2 + 1 or I2 == I1 + 1:
                    # update boolean
                    FoundAdjacent = True
            # update counter if at least 1 pair of adjacent orthologs exist in mouse
            if FoundAdjacent == True:
                ConservedPairs[GeneType] += 1


for i in range(len(JsonFiles)):
    ConservedPairs[JsonFiles[i]] = 0
    # count conserved pairs and adjacent pairs
    for pair in HumanAdjacentgenePairs[JsonFiles[i]]:
        if set(pair) in PairsOrthos[i]:
            # count conserved pairs
            ConservedPairs[JsonFiles[i]] += 1
        else:
            # check if orthologs are adjacent
            pair = list(pair)
            # make a list of orthologs in mouse
            L = []
            for ortho1 in Orthos[pair[0]]:
                for ortho2 in Orthos[pair[1]]:
                    # check that both genes have coordinates in mouse
                    if ortho1 in MouseCoord and ortho2 in MouseCoord:
                        # check that genes are different and on the same chromo 
                        if ortho1 != ortho2 and MouseCoord[ortho1][0] == MouseCoord[ortho2][0]:
                            L.append([ortho1, ortho2])
            # check if a pair of orthologs exist
            if len(L) != 0:
                # set up boolean to be changed if an adjacent pair is found
                FoundAdjacent = False
                # loop over the mouse gene pairs
                for j in range(len(L)):
                    gene1, gene2 = L[j][0], L[j][1]
                    # get indices of gene1 and gene2
                    I1, I2 = MouseOrdered[MouseCoord[gene1][0]].index(gene1), MouseOrdered[MouseCoord[gene2][0]].index(gene2)
                    # check if genes are adjacent
                    if I1 == I2 + 1 or I2 == I1 + 1:
                        # update boolean
                        FoundAdjacent = True
                # update counter if at least 1 pair of adjacent orthologs exist in mouse
                if FoundAdjacent == True:
                    ConservedPairs[JsonFiles[i]] += 1

              
   
# test differences among gene categories
# write P values to file
newfile = open('AdjacentPairsPvalues.txt', 'w')
GeneTypes = JsonFiles + ['Proximal', 'Moderate', 'Intermediate', 'Distant']
Header = '\t'.join(['\t']+GeneTypes)
newfile.write(Header + '\n')
for i in range(len(GeneTypes)):
    # create a list with pvalues for each line
    Line = [GeneTypes[i]]
    for j in range(len(GeneTypes)):
        # add p values for each pairwise comparison
        Line.append(str(stats.fisher_exact([[ConservedPairs[GeneTypes[i]], len(HumanAdjacentgenePairs[GeneTypes[i]])-ConservedPairs[GeneTypes[i]]], [ConservedPairs[GeneTypes[j]], len(HumanAdjacentgenePairs[GeneTypes[j]])-ConservedPairs[GeneTypes[j]]]])[1]))
    newfile.write('\t'.join(Line) + '\n')
newfile.close()
   
   
# compute proportions
for i in ConservedPairs:
    ConservedPairs[i] = ConservedPairs[i] / len(HumanAdjacentgenePairs[i])    
   
   
GeneTypes = JsonFiles + ['Proximal', 'Moderate', 'Intermediate', 'Distant']   
   

# create figure
figure = plt.figure(1, figsize = (2.5, 2))
# add a plot to figure (N row, N column, plot N)
ax = figure.add_subplot(1, 1, 1)

# make a list of gene category names parallel to the list of gene pairs
GeneCats = ['Nst', 'Pbk', 'Con', 'Div', 'Prx', 'Mod', 'Int', 'Dst']
# set colors
colorscheme = ['#f03b20', '#43a2ca', '#fee391', '#74c476', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
# plot proportions of gene pairs
ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], [ConservedPairs[i] for i in GeneTypes], 0.2, color = colorscheme,
        edgecolor = 'black', linewidth = 0.7)
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
## add a range for the Y and X axes
plt.ylim([0, 1])
# edit y axis ticks
plt.yticks(np.arange(0, 1.2, 0.2)) 
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
      
# save figure
for extension in ['.pdf', '.eps', '.png']:
    figure.savefig('AdjacentPairs' + extension, bbox_inches = 'tight')


