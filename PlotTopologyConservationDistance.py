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
# get the human gene pairs
HumanPairs = AllPairs[:5]

print('made lists of gene pairs')


# make a list of sets of gene pairs for the 2nd species
Sp2Pairs = []
for i in range(5, len(AllPairs)):
    Sp2Pairs.append([set(j) for j in AllPairs[i]])
print('made sets of gene pairs in species 2')        
        
        
# do qc
for i in range(1, len(Sp2Pairs)):
    for pair in Sp2Pairs[i]:
        assert pair in Sp2Pairs[0]
print('done with QC')

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos or pair[1] not in Orthos:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
print('removed human gene pairs lacking orthologs')
for i in HumanPairs:
    print(len(i))

# make lists of gene pairs in human [[gene1, gene2], ....[gene n, gene n+1]]
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
            if HumanOrdered[chromo][i] in Orthos and HumanOrdered[chromo][j] in Orthos and k in range(4):
                #HsaPairsDist[k].append([HumanOrdered[chromo][i], HumanOrdered[chromo][j]])
                
                HsaPairsDist[k].append({Orthos[HumanOrdered[chromo][i]], Orthos[HumanOrdered[chromo][j]]})
   


for i in HsaPairsDist:
    print(len(i))

         
print('generated human gene pairs by distance')












# make lists of sets of gene pairs in species 2 [{gene1, gene2}, ....{gene n, gene n+1}]
Sp2PairsDist = [[], [], [], []]
# loop over chromosomes
for chromo in Sp2Ordered:
    # loop over the list of ordered genes
    for i in range(len(Sp2Ordered[chromo]) - 1):
        # get the end position of gene 1
        EndGene1 = Sp2Coord[Sp2Ordered[chromo][i]][2]                
        # grab 2nd gene to form a pair                
        for j in range(i+1, len(Sp2Ordered[chromo])):
            # get the start position of gene 2
            StartGene2 = Sp2Coord[Sp2Ordered[chromo][j]][1]
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
                Sp2PairsDist[k].append({Sp2Ordered[chromo][i], Sp2Ordered[chromo][j]})

    

print('generated species 2 gene pairs by distance')


# add the pairs of non-overlapping genes to the lists of gene pairs
#HumanPairs.extend(HsaPairsDist)

for i in range(len(HumanPairs)):
    for j in range(len(HumanPairs[i])):
        HumanPairs[i][j][0] = Orthos[HumanPairs[i][j][0]]
        HumanPairs[i][j][1] = Orthos[HumanPairs[i][j][1]]        
        HumanPairs[i][j] = set(HumanPairs[i][j])
    print(HumanPairs[i][:10])
        
HumanPairs.extend(HsaPairsDist)
    toadd = []    
    for pair in HumanPairs[i]:
        toadd.append({Orthos[pair[0]], Orthos[pair[1]]})
    HumanOrthos.append(toadd)
HumanOrthos.extend(HsaPairsDist)



Sp2Pairs.extend(Sp2PairsDist)
print(len(HumanOrthos), len(Sp2Pairs))

HumanGenes, Sp2Genes = [], []

for i in range(len(HumanOrthos)):
    HumanGenes.append(np.array(HumanOrthos[i]))
    Sp2Genes.append(np.array(Sp2Pairs[i]))


# create a list of overlapping gene categories parallel to the list of overlapping pairs
GeneCats = ['overlapping', 'nested', 'piggyback', 'convergent', 'divergent',
            'proximal', 'moderate', 'intermediate', 'distant']

# create dictionary of gene pairs for each overlapping categories
HsaGenes = {}
for i in range(len(GeneCats)):
    HsaGenes[GeneCats[i]] = HumanPairs[i]

# create a dictionary of gene pairs without any order, for each overlapping gene categories
Sp2Genes = {}
for i in range(len(GeneCats)):
    Sp2Genes[GeneCats[i]] = Sp2Pairs[i]
print('generated dictionaries')


for i in range(len(GeneCats)):
    print(GeneCats[i], len(HsaGenes[GeneCats[i]]), len(Sp2Genes[GeneCats[i]]))

for i in range(len(HumanGenes)):
    print(i, 'hsa', sum(np.in1d(HumanGenes[i], HumanGenes[i], invert = False)))
    print(i, 'sp2', sum(np.in1d(HumanGenes[i], Sp2Genes[i], invert = False)))









#
#for i in HsaPairsDist:
#    print(len(i))
#for i in Sp2PairsDist:
#    print(len(i))
#
#
#
#
#
## create a list of overlapping gene categories parallel to the list of overlapping pairs
#GeneCats = ['overlapping', 'nested', 'piggyback', 'convergent', 'divergent',
#            'proximal', 'moderate', 'intermediate', 'distant']
#
## create dictionary of gene pairs for each overlapping categories
#HsaGenes = {}
#for i in range(len(GeneCats)):
#    #HsaGenes[GeneCats[i]] = HumanPairs[i]
#    HsaGenes[GeneCats[i]] = HumanOrthos[i]
#
#
#
## create a dictionary of gene pairs without any order, for each overlapping gene categories
#Sp2Genes = {}
#for i in range(len(GeneCats)):
#    Sp2Genes[GeneCats[i]] = Sp2Pairs[i]
#print('generated dictionaries')
#
#
#for i in range(len(GeneCats)):
#    print(GeneCats[i], len(HsaGenes[GeneCats[i]]), len(Sp2Genes[GeneCats[i]]))
#
#
#
#
### count the number of gene pairs for which orthologs in are the same topology
##PairCounts = {}
##for i in range(len(GeneCats) -1):
##    # initialize counter
##    total = 0
##    # loop over gene pairs for the given gene category
##    for pair in HsaGenes[GeneCats[i]]:
##        # check if pair has same topology
##        if set([Orthos[pair[0]], Orthos[pair[1]]]) in Sp2Genes[GeneCats[i]]:
##            total += 1
##    # populate dict
##    PairCounts[GeneCats[i]] = [total, len(HsaGenes[GeneCats[i]])]
#
#
#
#
#
#
## count the number of gene pairs for which orthologs in are the same topology
#PairCounts = {}
#for i in range(len(GeneCats) -1):
#    # initialize counter
#    total = 0
#    # loop over gene pairs for the given gene category
#    for pair in HsaGenes[GeneCats[i]]:
#        # check if pair has same topology
#        if pair in Sp2Genes[GeneCats[i]]:
#            total += 1
#    # populate dict
#    PairCounts[GeneCats[i]] = total
#
#
#for i in HsaGenes:
#    HsaGenes[i] = np.array(HsaGenes[i])
#for i in Sp2Genes:
#    Sp2Genes[i] = np.array(Sp2Genes[i])
#
#
#for i in HsaGenes:
#    print(len(HsaGenes[i]))
#for i in Sp2Genes:
#    print(len(Sp2Genes[i]))
#
#
#CountPairs = {}
#for i in range(len(GeneCats) -1):
#    truc = 0
#    for j in HsaGenes[GeneCats[i]]:
#        if j in Sp2Genes[GeneCats[i]]:
#            truc += 1
#    #total = np.in1d(HsaGenes[GeneCats[i]], Sp2Genes[GeneCats[i]], invert = False)
#    print(GeneCats[i], HsaGenes[GeneCats[i]][:10], Sp2Genes[GeneCats[i]][:10])
#    total = np.in1d(HsaGenes[GeneCats[i]], HsaGenes[GeneCats[i]], invert = False)
#    
#    print(GeneCats[i], truc, sum(total))    
#    CountPairs[GeneCats[i]] = sum(total)
#
#for i in range(len(GeneCats) -1):
#    print(GeneCats[i], PairCounts[GeneCats[i]], CountPairs[GeneCats[i]])




## count the number of gene pairs for which orthologs in are the same topology
#PairCounts = {}
#for i in range(len(GeneCats) -1):
#    # initialize counter
#    total = 0
#    # loop over gene pairs for the given gene category
#    for pair in HsaGenes[GeneCats[i]]:
#        # check if pair has same topology
#        if set([Orthos[pair[0]], Orthos[pair[1]]]) in Sp2Genes[GeneCats[i]]:
#            total += 1
#    # populate dict
#    PairCounts[GeneCats[i]] = [total, len(HsaGenes[GeneCats[i]])]



# count the number of gene pairs for which orthologs in are the same topology
PairCounts = {}
for i in range(len(GeneCats) -1):
    # initialize counter
    total = 0
    # loop over gene pairs for the given gene category
    for pair in HsaGenes[GeneCats[i]]:
        # check if pair has same topology
        if pair in Sp2Genes[GeneCats[i]]:
            total += 1
    # populate dict
    PairCounts[GeneCats[i]] = [total, len(HsaGenes[GeneCats[i]])]


for i in HsaGenes:
    HsaGenes[i] = np.array(HsaGenes[i])
for i in Sp2Genes:
    Sp2Genes[i] = np.array(Sp2Genes[i])



CountPairs = {}
for i in range(len(GeneCats) - 1):
    print('array', i)
    CountPairs[GeneCats[i]] = sum(np.in1d(HsaGenes[GeneCats[i]], Sp2Genes[GeneCats[i]], invert = False))


for i in range(len(GeneCats) - 1):
    print(i, PairCounts[GeneCats[i]], CountPairs[GeneCats[i]], PairCounts[GeneCats[i]] == CountPairs[GeneCats[i]])



#newfile = open('PairsCounts.txt', 'a')
#header = '\t'.join(['species', 'genes', 'conserved', 'total', 'ratio'])
#newfile.write(header + '\n')
#
#for i in range(len(GeneCats)):
#    line = '\t'.join([species, GeneCats[i], str(PairCounts[GeneCats[i]][0]), str(PairCounts[GeneCats[i]][1]), str(PairCounts[GeneCats[i]][0] / PairCounts[GeneCats[i]][1])])
#    newfile.write(line + '\n')
#newfile.close()


a = np.array([{'MouseNestedGenes.json', 'ChimpOverlappingGenes.json'},
       {'MouseDivergentGenes.json', 'HumanConvergentGenes.json'},
       {'ChimpConvergentGenes.json', 'HumanPiggyBackGenes.json'},
       {'MouseOverlappingGenes.json', 'HumanNestedGenes.json'},
       {'ChimpDivergentGenes.json', 'ChimpNestedGenes.json'},
       {'MouseConvergentGenes.json', 'HumanDivergentGenes.json'},
       {'MousePiggyBackGenes.json', 'ChimpPiggyBackGenes.json'}])
print(sum(np.in1d(a, a)))
print(sum(np.in1d(a, a, invert = False)) / len(a))




