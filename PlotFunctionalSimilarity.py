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
             'HumanPiggyBackGenes.json', 'HumanConvergentGenes.json',
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
        JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
        GOOverlap.append(JI)
    FunctionalSimilarity.append(GOOverlap)
        
        
# generate random pairs of non-overlapping genes as a control
# convert set of non-overlapping genes to list
NonOverlappingGenes = list(NonOverlappingGenes)
ControlPairs = []
# draw 10000 random pairs
for i in range(10000):
    # pick 2 random genes
    j = random.randint(0, len(NonOverlappingGenes) -1)
    k = random.randint(0, len(NonOverlappingGenes) -1)
    gene1, gene2 = NonOverlappingGenes[j], NonOverlappingGenes[k]
    JI = JaccardIndex(GeneOntology[gene1], GeneOntology[gene2])
    ControlPairs.append(JI)

# insert pairs of JI in list
FunctionalSimilarity.insert(0, ControlPairs)

for i in range(len(FunctionalSimilarity)):
    print(np.mean(FunctionalSimilarity[i]))
















