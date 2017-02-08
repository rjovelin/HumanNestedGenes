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



# load dictionaries of overlapping genes
jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json', 'HumanPiggyBackGenes.json',
             'HumanConvergentGenes.json', 'HumanDivergentGenes.json', 
             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json', 'ChimpPiggyBackGenes.json',
             'ChimpConvergentGenes.json', 'ChimpDivergentGenes.json']






# make a list of dictionaries
AllOverlap = []
# loop over files
for i in range(len(jsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(jsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    AllOverlap.append(overlapping)

# assign dicts to variables
HsaOverlapping, HsaNested, HsaPiggyback = AllOverlap[0], AllOverlap[1], AllOverlap[2]
HsaConvergent, HsaDivergent = AllOverlap[3], AllOverlap[4]

PtrOverlapping, PtrNested, PtrPiggyback = AllOverlap[5], AllOverlap[6], AllOverlap[7]
PtrConvergent, PtrDivergent = AllOverlap[8], AllOverlap[9]

# get GFF file
GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3']



# make a list of gene coordinates       
AllCoordinates = []
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
    AllCoordinates.append(GeneCoord)

HsaGeneCoord, PtrGeneCoord = AllCoordinates[0], AllCoordinates[1]

# make a list of gene sets
GeneSets = []
for i in range(len(AllOverlap)):
    overlapping = MakeFullPartialOverlapGeneSet(AllOverlap[i])
    GeneSets.append(overlapping)

HsaOverlappingGenes, HsaNestedGenes, HsaPiggyBackGenes, HsaConvergentGenes, HsaDivergentGenes = GeneSets[0], GeneSets[1], GeneSets[2], GeneSets[3], GeneSets[4]
PtrOverlappingGenes, PtrNestedGenes, PtrPiggyBackGenes, PtrConvergentGenes, PtrDivergentGenes = GeneSets[5], GeneSets[6], GeneSets[7], GeneSets[8], GeneSets[9]

# get 1:1 orthologs between human and chimp
Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')



Counts = []
for i in range(len(GeneSets[:5])):
    total = [0,0]
    for j in GeneSets[i]:
        if j in Orthos:
            if Orthos[j] in GeneSets[i+5]:
                total[0] += 1
            else:
                if Orthos[j] in GeneSets[5]:
                    total[1] += 1
    Counts.append(total)
  
for i in range(len(GeneSets[:5])):
    print(i, len(GeneSets[i]), Counts[i])




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
    
