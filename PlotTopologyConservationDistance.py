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
HumanOrdered, Sp2Ordered = AllOrdered[0], AllOrdered[1]
HumanCoord, Sp2Coord = AllCoordinates[0], AllCoordinates[1]


# get 1:1 orthologs between human and other species
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
HumanPairs = AllPairs[:5]
Sp2Pairs = AllPairs[5:]

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos or pair[1] not in Orthos:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)





def MakeGenePairsDistance(GeneCoord, OrderedGenes, Orthos, species = 'human'):
    '''
    (dict, dict) -> tuple
    Take the dictionary of gene coordinates, the dictionary of ordered genes
    along each chromosome, and return a tuple with lists of gene pairs
    separated by a given distance (in bp): proximal, intermediate and distant
    '''
        
    # make lists of gene pairs [[gene1, gene2], ....[gene n, gene n+1]]
    Proximal, Moderate, Intermediate, Distant = [], [], [], []
    # loop over chromosomes
    for chromo in OrderedGenes:
        # loop over the list of ordered genes
        for i in range(len(OrderedGenes[chromo]) - 1):
            # get the end position of gene 1
            EndGene1 = GeneCoord[OrderedGenes[chromo][i]][2]                
            # grab 2nd gene to form a pair                
            for j in range(i+1, len(OrderedGenes[chromo])):
                # get the start position of gene 2
                StartGene2 = GeneCoord[OrderedGenes[chromo][j]][1]
                # check if distance is less that 500 bp
                D = StartGene2 - EndGene1
                if D >= 0 and D < 1000:
                    # add gene pair to Proximal
                    Proximal.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
                elif D >= 1000 and D < 10000:
                    # add gene pair to Intermediate
                    Moderate.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
                elif D >= 10000 and D < 50000:
                    # add gene pair to Intermediate
                    Intermediate.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
                elif D >= 50000:
                    # add gene pair to Distant
                    Distant.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
    return Proximal, Moderate, Intermediate, Distant


#########################################

# make lists of gene pairs [[gene1, gene2], ....[gene n, gene n+1]]
#HsaProximal, HsaModerate, HsaIntermediate, HsaDistant = [], [], [], []
HsaPairsDist = [[], [], [], []]
# loop over chromosomes
for chromo in HumanOrdered:
    # loop over the list of ordered genes
    for i in range(len(HumanOrdered[chromo]) - 1):
        # get the end position of gene 1
        EndGene1 = HumanCoord[HumanOrdered[chromo][i]][2]                
        # grab 2nd gene to form a pair                
        for j in range(i+1, len(HumanOrdered[chromo])):
            # get the start position of gene 2
            StartGene2 = HumanCoord[HumanOrdered[chromo][j]][1]
            # check if distance is less that 500 bp
            D = StartGene2 - EndGene1
            if D >= 0 and D < 1000:
                # add gene pair to Proximal
                k = 0
            elif D >= 1000 and D < 10000:
                # add gene pair to Intermediate
                k = 1                
            elif D >= 10000 and D < 50000:
                # add gene pair to Intermediate
                k = 2
            elif D >= 50000:
                # add gene pair to Distant
                k = 3
            # check that human genes in pair have orthologs    
            if HumanOrdered[chromo][i] in Orthos and HumanOrdered[chromo][j] in Orthos:
                HsaPairsDist[k].append([HumanOrdered[chromo][i], HumanOrdered[chromo][j]])
            












########################################

DistantPairs = MakeGenePairsDistance(HumanCoord, HumanOrdered)
for i in DistantPairs:
    HumanPairs.append(i)
DistantPairs = MakeGenePairsDistance(Sp2Coord, Sp2Ordered)
for i in DistantPairs:
    Sp2Pairs.append(i)



# create a list of overlapping gene categories parallel to the list of overlapping pairs
GeneCats = ['overlapping', 'nested', 'piggyback', 'convergent', 'divergent']
# create dictionary of gene pairs for each overlapping categories
HsaGenes = {}
for i in range(len(GeneCats)):
    HsaGenes[GeneCats[i]] = HumanPairs[i]
# create a dictionary of gene pairs without any order, for each overlapping gene categories
Sp2Genes = {}
for i in range(len(GeneCats)):
    Sp2Genes[GeneCats[i]] = [set(j) for j in Sp2Pairs[i]]
     
# do qc
for i in range(1, len(GeneCats)):
    for pair in Sp2Genes[GeneCats[i]]:
        assert pair in Sp2Genes['overlapping']

