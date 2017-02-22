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

print('ordered genes')
print('got gene coordinates')


# get 1:1 orthologs between human and other species
if species == 'chimp':
    Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
elif species == 'mouse':
    Orthos = MatchOrthologPairs('HumanMouseOrthologs.txt')

print('mapped orthologs')


# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# get the gene pairs
HumanPairs = AllPairs[:5]
Sp2Pairs = AllPairs[5:]
print('made lists of gene pairs')

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos or pair[1] not in Orthos:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
print('removed human gene pairs lacking orthologs')

# replace human genes by their orthologs
for i in range(len(HumanPairs)):
    for j in range(len(HumanPairs[i])):
        HumanPairs[i][j][0] = Orthos[HumanPairs[i][j][0]]
        HumanPairs[i][j][1] = Orthos[HumanPairs[i][j][1]]

# sort gene pairs
for i in range(len(HumanPairs)):
    for j in range(len(HumanPairs[i])):
        HumanPairs[i][j].sort()
for i in range(len(Sp2Pairs)):
    for j in range(len(Sp2Pairs[i])):
        Sp2Pairs[i][j].sort()

# replace gene lists with strings
for i in range(len(HumanPairs)):
    for j in range(len(HumanPairs[i])):
        HumanPairs[i][j] = ':'.join(HumanPairs[i][j])
for i in range(len(Sp2Pairs)):
    for j in range(len(Sp2Pairs[i])):
        Sp2Pairs[i][j] = ':'.join(Sp2Pairs[i][j])

       
# do qc
for i in range(1, len(Sp2Pairs)):
    for pair in Sp2Pairs[i]:
        assert pair in Sp2Pairs[0]
print('done with QC')


# make lists of gene pairs in human [[gene1, gene2], ....[gene n, gene n+1]]
HsaPairsDist = [[], [], [], []]
# loop over chromosomes
for chromo in HumanOrdered:
    # loop over the list of ordered genes
    for i in range(len(HumanOrdered[chromo]) - 1):
        # get the end position of gene 1
        EndGene1 = HumanCoord[HumanOrdered[chromo][i]][2]                
        # get the start position of adjacent gene 2
        StartGene2 = HumanCoord[HumanOrdered[chromo][i+1]][1]
        # check if distance is less that 500 bp
        D = StartGene2 - EndGene1
        # assign infinity value to k
        k = float('inf')
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
        if HumanOrdered[chromo][i] in Orthos and HumanOrdered[chromo][i+1] in Orthos and k in range(4):
            # add the human orthologs
            pair = [Orthos[HumanOrdered[chromo][i]], Orthos[HumanOrdered[chromo][i+1]]]
            # sort gene pair
            pair.sort()
            # concert list to string
            pair = ':'.join(pair)
            HsaPairsDist[k].append(pair)

print('generated human gene pairs by distance')


# make lists of sets of gene pairs in species 2 [{gene1, gene2}, ....{gene n, gene n+1}]
Sp2PairsDist = [[], [], [], []]
# loop over chromosomes
for chromo in Sp2Ordered:
    # loop over the list of ordered genes
    for i in range(len(Sp2Ordered[chromo]) - 1):
        # get the end position of gene 1
        EndGene1 = Sp2Coord[Sp2Ordered[chromo][i]][2]                
        # get the start position of adjacent gene 2
        StartGene2 = Sp2Coord[Sp2Ordered[chromo][i+1]][1]
        # check if distance is less that 500 bp
        D = StartGene2 - EndGene1
        # assign infinity value to k
        k = float('inf')
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
        # populate lists with sets of gene pairs    
        if k in range(4):
            # get gene pair
            pair = [Sp2Ordered[chromo][i], Sp2Ordered[chromo][i+1]]
            # sort pair
            pair.sort()
            # convert list to string
            pair = ':'.join(pair)
            Sp2PairsDist[k].append(pair)

print('generated species 2 gene pairs by distance')


# add the pairs of non-overlapping genes to the lists of gene pairs
HumanPairs.extend(HsaPairsDist)
Sp2Pairs.extend(Sp2PairsDist)


# convert lists to numpy arrays
for i in range(len(HumanPairs)):
    HumanPairs[i] = np.array(HumanPairs[i])
for i in range(len(Sp2Pairs)):
    Sp2Pairs[i] = np.array(Sp2Pairs[i])


for i in range(len(HumanPairs)):
    print(i, len(HumanPairs[i]), len(Sp2Pairs[i]))





# count the number of pairs with conserved topology
CountPairs = []
for i in range(len(HumanPairs)):
    total = sum(np.in1d(HumanPairs[i], Sp2Pairs[i], invert = False))    
    CountPairs.append([total, len(HumanPairs[i])])

# create a list of overlapping gene categories parallel to the list of overlapping pairs
GeneCats = ['overlapping', 'nested', 'piggyback', 'convergent', 'divergent',
            'proximal', 'moderate', 'intermediate', 'distant']

for i in range(len(GeneCats)):
    print(GeneCats[i], CountPairs[i][0] / CountPairs[i][1])


#newfile = open('PairsCounts.txt', 'a')
#header = '\t'.join(['species', 'genes', 'conserved', 'total', 'ratio'])
#newfile.write(header + '\n')
#
#for i in range(len(GeneCats)):
#    line = '\t'.join([species, GeneCats[i], str(PairCounts[GeneCats[i]][0]), str(PairCounts[GeneCats[i]][1]), str(PairCounts[GeneCats[i]][0] / PairCounts[GeneCats[i]][1])])
#    newfile.write(line + '\n')
#newfile.close()


#newfile = open('test.txt', 'w')
#newfile.write('a\n')
#newfile.write(str(sum(np.in1d(a, a))) + '\n')
#newfile.write(str(sum(np.in1d(a, a, invert = False)) / len(a)) + '\n')
#newfile.write('\n\n')
#
#
#truc = '\t'.join([str(sum(np.in1d(b, b))), str(sum(np.in1d(b, b[:15]))), str(sum(np.in1d(b, b[:20]))), str(sum(np.in1d(b, b[:30]))), str(sum(np.in1d(b, b[:50])))])
#newfile.write(truc + '\n')
#newfile.close()

