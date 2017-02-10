# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:46:17 2017

@author: RJovelin
"""

# use this script to plot the % of orthologous gene pairs with same topology between human and chimp

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

species = sys.argv[1]

# load dictionaries of overlapping genes

if species == 'chimp':
    jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json', 'HumanPiggyBackGenes.json',
                 'HumanConvergentGenes.json', 'HumanDivergentGenes.json', 
                 'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json', 'ChimpPiggyBackGenes.json',
                 'ChimpConvergentGenes.json', 'ChimpDivergentGenes.json']
elif species == 'mouse':
    jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json', 'HumanPiggyBackGenes.json',
                 'HumanConvergentGenes.json', 'HumanDivergentGenes.json', 
                 'MouseOverlappingGenes.json', 'MouseNestedGenes.json', 'MousePiggyBackGenes.json',
                 'MouseConvergentGenes.json', 'MouseDivergentGenes.json']


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
if species == 'chimp':
    GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3']
elif species == 'mouse':
    GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Mus_musculus.GRCm38.86.gff3']


# make a list of gene coordinates       
AllCoordinates = []
AllOrdered = []

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
HumanOrdered = AllOrdered[0]
Sp2Ordered = AllOrdered[1]
HumanCoord = AllCoordinates[0]
Sp2Coord = AllCoordinates[1]


# make a list of gene sets
GeneSets = []
for i in range(len(AllOverlap)):
    overlapping = MakeFullPartialOverlapGeneSet(AllOverlap[i])
    GeneSets.append(overlapping)

# get 1:1 orthologs between human and chimp
if species == 'chimp':
    Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
elif species == 'mouse':
    Orthos = MatchOrthologPairs('HumanMouseOrthologs.txt')

#Counts = []
#for i in range(len(GeneSets[:5])):
#    total = [0,0,0]
#    for j in GeneSets[i]:
#        if j in Orthos:
#            total[0] += 1
#            if Orthos[j] in GeneSets[i+5]:
#                total[1] += 1
#            else:
#                if Orthos[j] in GeneSets[5]:
#                    total[2] += 1
#    Counts.append(total)
#  
#for i in range(len(GeneSets[:5])):
#    print(i, len(GeneSets[i]), Counts[i])


# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
for i in AllPairs:
    print(len(i), end = ' ')
print('\n')


# remove human genes lacking orthologs
for i in range(len(AllPairs[:5])):
    to_remove = []
    for pair in AllPairs[i]:
        if pair[0] not in Orthos or pair[1] not in Orthos:
            to_remove.append(pair)
    for pair in to_remove:
        AllPairs[i].remove(pair)
for i in AllPairs:
    print(len(i), end = ' ')
print('\n')


# count the number of genes of non-overlapping orthologs of human overlapping genes
AllCounts = []
HumanPairs, Sp2Pairs = [], []
for i in range(len(AllPairs)):
    if i < 5:
        pairs = [j for j in AllPairs[i]]
        HumanPairs.append(pairs)
    else:
        pairs = [set(j) for j in AllPairs[i]]
        Sp2Pairs.append(pairs)

NotOverlap = []
Missing = []
for i in range(len(HumanPairs)):
    L = []
    total = 0
    missing = 0
    for pair in HumanPairs[i]:
        # get the orthologs
        orthos = set()
        orthos.add(Orthos[pair[0]])
        orthos.add(Orthos[pair[1]])
        if orthos not in Sp2Pairs[0]:
            total += 1
            if Orthos[pair[0]] in Sp2Coord and Orthos[pair[1]] in Sp2Coord:
                # get the coordinates of the non-overlapping orthologs
                chromo1 = Sp2Coord[Orthos[pair[0]]][0]
                chromo2 = Sp2Coord[Orthos[pair[1]]][0]
                if chromo1 != chromo2:
                    L.append(1000)
                else:
                    print(Orthos[pair[0]], Orthos[pair[1]])
                    # get the index of both genes on chromo
                    j = Sp2Ordered[chromo1].index(Orthos[pair[0]])
                    k = Sp2Ordered[chromo1].index(Orthos[pair[1]])
                    if k > j:
                        L.append(k - j - 1)
                    elif k < j:
                        L.append(j - k - 1)
            else:
                missing += 1
    AllCounts.append(L)
    NotOverlap.append(total)
    Missing.append(missing)
for i in range(len(AllCounts)):
    print(i, len(AllCounts[i]), min(AllCounts[i]), AllCounts[i].count(min(AllCounts[i])), max(AllCounts[i]), AllCounts[i].count(max(AllCounts[i])), np.mean(AllCounts[i]), np.median(AllCounts[i]))

print(NotOverlap)
print(Missing)


#############


## create lists of gene pairs
#OverlappingPairs = GetHostNestedPairs(Overlapping)
#NestedPairs = GetHostNestedPairs(Nested)
#PiggybackPairs = GetHostNestedPairs(Piggyback)
#ConvergentPairs = GetHostNestedPairs(Convergent)
#DivergentPairs = GetHostNestedPairs(Divergent)
#
#
## use this function to generate sets of gene pairs separated by a given distance
#def GenerateSetsGenePairsDistance(GeneCoord, OrderedGenes, ExpressionProfile):
#    '''
#    (dict, dict, dict) -> tuple
#    Take the dictionary of gene coordinates, the dictionary of ordered genes
#    along each chromosome, the dictionary of gene: expression pairs and return a
#    tuple with lists of gene pairs separated by a given distance (in bp):
#    proximal, intermediate and distant
#    '''
#        
#    # make lists of gene pairs [[gene1, gene2], ....[gene n, gene n+1]]
#    Proximal, Moderate, Intermediate, Distant = [], [], [], []
#    # loop over chromosomes
#    for chromo in OrderedGenes:
#        # loop over the list of ordered genes
#        for i in range(len(OrderedGenes[chromo])):
#            # check that gene is not host or nested, has expression
#            if OrderedGenes[chromo][i] in ExpressionProfile:
#                # get the end position of gene 1
#                EndGene1 = GeneCoord[OrderedGenes[chromo][i]][2]                
#                # grab 2nd gene to form a pair                
#                for j in range(i+1, len(OrderedGenes[chromo])):
#                    # check that gene is not host or nested and has expression
#                    if OrderedGenes[chromo][j] in ExpressionProfile:
#                        # get the start position of gene 2
#                        StartGene2 = GeneCoord[OrderedGenes[chromo][j]][1]
#                        # check if distance is less that 500 bp
#                        D = StartGene2 - EndGene1
#                        if D >= 0 and D < 1000:
#                            # add gene pair to Proximal
#                            Proximal.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
#                        elif D >= 1000 and D < 10000:
#                            # add gene pair to Intermediate
#                            Moderate.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
#                        elif D >= 10000 and D < 50000:
#                            # add gene pair to Intermediate
#                            Intermediate.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
#                        elif D >= 50000:
#                            # add gene pair to Distant
#                            Distant.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
#    return Proximal, Moderate, Intermediate, Distant
#    
    
