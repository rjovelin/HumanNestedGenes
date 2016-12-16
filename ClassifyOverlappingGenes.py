# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:53:18 2016

@author: RJovelin
"""

# use this script to classify overlapping genes into 4 groups

# usage ClassifyOverlappingGenes.py [options]
# [human/chimp/gorilla]: species to consider

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
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


# get the species to consider from the command line
CurrSpecies = sys.argv[1]
assert CurrSpecies in ['human', 'chimp', 'gorilla'], 'species name is not valid'

# load dictionary of overlapping genes
if CurrSpecies == 'human':
    json_data = open('HumanOverlappingGenes.json')
elif CurrSpecies == 'chimp':
    json_data = open('ChimpOverlappingGenes.json')
elif CurrSpecies == 'gorilla':
    json_data = open('GorillaOverlappingGenes.json')
# load dictionary
OverlappingGenes = json.load(json_data)
json_data.close()

# make pairs of overlapping genes
OverlappingPairs = GetHostNestedPairs(OverlappingGenes)

# get GFF file
if CurrSpecies == 'human':
    GFF = 'Homo_sapiens.GRCh38.86.gff3'
elif CurrSpecies == 'chimp':
    GFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
elif CurrSpecies == 'gorilla':
    GFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
    
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
# Map Transcript names to gene names {transcript: gene}
MapTranscriptGene = TranscriptToGene(GFF)
# get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
ExonCoord = GeneExonCoord(GFF)
ExonCoord = CleanGeneFeatureCoord(ExonCoord, MapTranscriptGene)
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
IntronCoord = GeneIntronCoord(ExonCoord)
IntronCoord = CleanGeneFeatureCoord(IntronCoord, MapTranscriptGene)
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
CombinedIntronCoord = CombineAllGeneRegions(IntronCoord, MapTranscriptGene)



# define 4 categories of overlapping genes:
# 1) nested gene pairs: one gene is fully contained within the intron of another gene
#    both genes can be on the same strand or on opposite strands
# 2) piggyback gene pairs: both gene have the same orientation
#    overlap can be partial or complete
# 3) convergent gene pairs: both genes have different orientation
#    first gene on chromosome is +, second gene on chromo is -
# 4) divergent gene pairs: both genes have different orientation
#    first gene on chromosome is -, second gene on chromo is +

# create lists for groups of overlapping genes
Nested, Piggyback, Convergent, Divergent = [], [], [], []

# find intronic nested genes first
ContainedGenes = FindContainedGenePairs(GeneCoord, OverlappingGenes)
# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
HostGenes = FindIntronicNestedGenePairs(ContainedGenes, CombinedIntronCoord, GeneCoord)
# make a list of Host, nested gene pairs
Nested = GetHostNestedPairs(HostGenes)

# make a list of set of genes to remove nested relationships
NestedSets = []
for pair in Nested:
    NestedSets.append(set(pair))

# find piggyback genes, convergent and divergent genes
# loop over overlapping gene pairs
# first gene in the pair comes first on chromosome
# but sometimes both members of the pair have same starting position
for pair in OverlappingPairs:
    # check that pair is not nested
    if set(pair) not in NestedSets:
        # check orientation of both genes
        if GeneCoord[pair[0]][-1] == GeneCoord[pair[1]][-1]:
            # same orientation, piggyback gene pairs
            Piggyback.append(pair)
        elif GeneCoord[pair[0]][-1] != GeneCoord[pair[1]][-1]:
            # check orientation of first and second gene in the pair
            assert GeneCoord[pair[0]][1] <= GeneCoord[pair[1]][1]
            if GeneCoord[pair[0]][-1] == '+':
                assert GeneCoord[pair[1]][-1] == '-'
                # convergent gene pair
                Convergent.append(pair)
            elif GeneCoord[pair[0]][-1] == '-':
                assert GeneCoord[pair[1]][-1] == '+'
                # divergent gene pairs
                Divergent.append(pair)
            

print('nested', len(Nested))
print('piggyback', len(Piggyback))
print('convergent', len(Convergent))
print('divergent', len(Divergent))

assert len(OverlappingPairs) == len(Nested) + len(Piggyback) + len(Convergent) + len(Divergent)


# save overlapping relationships to json files

# save nested genes as json file
newfile = open(CurrSpecies.title() + 'NestedGenes.json', 'w')
json.dump(HostGenes, newfile, sort_keys = True, indent = 4)
newfile.close()

# create a dictionary with Piggyback genes
PiggyBackGenes = {}
for pair in Piggyback:
    # use first gene in pair as key
    if pair[0] in PiggyBackGenes:
        PiggyBackGenes[pair[0]].append(pair[1])
    else:
        PiggyBackGenes[pair[0]] = [pair[1]]
newfile = open(CurrSpecies.title() + 'PiggyBackGenes.json', 'w')
json.dump(PiggyBackGenes, newfile, sort_keys = True, indent = 4)
newfile.close()

# create a dictionary with convergent genes
ConvergentGenes = {}
for pair in Convergent:
    # use first gene in pair as key
    if pair[0] in ConvergentGenes:
        ConvergentGenes[pair[0]].append(pair[1])
    else:
        ConvergentGenes[pair[0]] = [pair[1]]
newfile = open(CurrSpecies.title() + 'ConvergentGenes.json', 'w')
json.dump(ConvergentGenes, newfile, sort_keys = True, indent = 4)
newfile.close()

# create a dictionary with divergent genes
DivergentGenes = {}
for pair in Divergent:
    # use first gene in pair as key
    if pair[0] in DivergentGenes:
        DivergentGenes[pair[0]].append(pair[1])
    else:
        DivergentGenes[pair[0]] = [pair[1]]
newfile = open(CurrSpecies.title() + 'DivergentGenes.json', 'w')
json.dump(DivergentGenes, newfile, sort_keys = True, indent = 4)
newfile.close()

