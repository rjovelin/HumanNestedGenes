# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 21:24:45 2017

@author: Richard
"""

# use this script to plot differences in functional similarities between overlapping genes


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



# load the overlapping gene dictionaries


# make a set of non-overlapping genes

# load dictionaries of overlapping genes
JsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',
             'HumanPiggybackGenes.json', 'HumanConvergentGenes.json',
             'HumanDivergentGenes.json']
# get GFF file
GFF = 'Homo_sapiens.GRCh38.88.gff3'

# make a list of dictionaries
Overlap = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    Overlap.append(overlapping)

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

# make pairs of overlapping genes
OverlappingPairs = []
for i in range(len(Overlap)):
    pairs = GetHostNestedPairs(Overlap[i])
    OverlappingPairs.append(pairs)

# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

# map gene names with gene ID {gene ID: name}
Names = MapNametoID(GFF)
# parse file with GO annotations 
GOAnnotations = ParseGOFile('goa_human.gaf')
# map gene IDs to GO annotations            
GeneOntology = MapEnsemblGenesToGOTerms(GOAnnotations, Names)


# create a list of lists with functional similarities between genes of the same pair
FunctionalSimilarity = []
for i in range(len(OverlappingPairs)):
    # compute jaccard similarity index between each pair
    GOOverlap = []
    for pair in OverlappingPairs[i]:
        JI = JaccardIndex(pair[0], pair[1])
        GOOverlap.append(JI)
    FunctionalSimilarity.append(GOOverlap)
        
        























#######################################





# make sets of host and nested nested genes
NestedSets = []
for i in range(1, len(AllOverlap), 2):
    nestedset = MakeFullPartialOverlapGeneSet(AllOverlap[i])
    NestedSets.append(nestedset)

# make sets of overlapping genes
OverlapSets = []
for i in range(0, len(AllOverlap), 2):
    overlap = MakeFullPartialOverlapGeneSet(AllOverlap[i])
    OverlapSets.append(overlap)







if Analysis == 'pairs':
    # compare expression divergence between human host and nested genes and their un-nested orthologs in sister-species   
    # remove human pairs if orthologs are nested in sister-species
    to_remove = [pair for pair in HumanYoung if OrthoPairs[pair[0]] in NestedSets[1] or OrthoPairs[pair[1]] in NestedSets[1]]
    for pair in to_remove:
        HumanYoung.remove(pair)
    # remove sister species pairs if orthologs are nested in human
    to_remove = [pair for pair in SisterSpYoung if SisterOrthos[pair[0]] in NestedSets[0] or SisterOrthos[pair[1]] in NestedSets[0]]
    for pair in to_remove:
        SisterSpYoung.remove(pair)



