# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:21:31 2017

@author: RJovelin
"""

# use this script to compare the proportions of genes with and without homologs

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


# load dictionary of overlapping gene pairs
json_data = open('HumanOverlappingGenes.json')
Overlapping = json.load(json_data)
json_data.close()
# load dictionary of nested gene pairs
json_data = open('HumanNestedGenes.json')
Nested = json.load(json_data)
json_data.close()
# load dictionary of pibbyback gene pairs
json_data = open('HumanPiggyBackGenes.json')
Piggyback = json.load(json_data)
json_data.close()
# load dictionary of convergent gene pairs
json_data = open('HumanConvergentGenes.json')
Convergent = json.load(json_data)
json_data.close()
# load dictionary of divergent gene pairs
json_data = open('HumanDivergentGenes.json')
Divergent = json.load(json_data)
json_data.close()

# get GFF file
GFF = 'Homo_sapiens.GRCh38.86.gff3'
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)

# generate gene sets
NestedGenes  = MakeFullPartialOverlapGeneSet(Nested)
OverlappingGenes = MakeFullPartialOverlapGeneSet(Overlapping)
ConvergentGenes = MakeFullPartialOverlapGeneSet(Convergent)
DivergentGenes = MakeFullPartialOverlapGeneSet(Divergent)
PiggyBackGenes = MakeFullPartialOverlapGeneSet(Piggyback)

# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)

# create sets of internal and external nested gene pairs
NestedPairs = GetHostNestedPairs(Nested)
InternalGenes, ExternalGenes = set(), set()
for pair in NestedPairs:
    ExternalGenes.add(pair[0])
    InternalGenes.add(pair[1])
    

# create a dictionary of homologs
Homologs = {}
infile = open('HumanChimpOrthologs.txt')
infile.readline()
# loop over file, match all orthologs, including 1:1 and 1: to many
for line in infile:
    if 'ortholog' in line:
        line = line.rstrip().split('\t')
        gene1, gene2 = line[0], line[2]
        assert 'ENS' in gene1, 'gene id is not valid'
        assert 'ortholog' in line[4], 'ortholog should be in homology type'
        if gene1 in Homologs:
            Homologs[gene1].add(gene2)
        else:
            Homologs[gene1] = set()
            Homologs[gene1].add(gene2)
infile.close()            


# create lists of orthologous pairs for each gene category 
GeneCats = ['NoOv', 'Nst', 'Int', 'Ext', 'Pbk', 'Conv', 'Div']
AllGenes = [NonOverlappingGenes, NestedGenes, InternalGenes, ExternalGenes,
            PiggyBackGenes, ConvergentGenes, DivergentGenes] 


# count genes with and without homologs
GeneCounts = []
# loop over gene sets
for i in range(len(AllGenes)):
    # initialize counters
    homo, nohomo = 0, 0
    # loop over genes in given set
    for gene in AllGenes[i]:
        if gene in Homologs:
            homo += 1
        else:
            nohomo += 1
    GeneCounts.append([homo, nohomo])    

# compare the proportions of gene with and without homologs
# create a list to store the P-values
PValues = ['']
for i in range(1, len(GeneCounts)):
    p = stats.fisher_exact([GeneCounts[0], GeneCounts[i]])[1]
    PValues.append(p)
PVals = []
for p in PValues:
    if p == '':
        PVals.append(p)
    elif p >= 0.05:
        PVals.append(str(round(p, 3)))
    elif p < 0.05 and p >= 0.01:
        PVals.append('< 0.05')
    elif p < 0.01 and p >= 0.001:
        PVals.append('< 0.01')
    elif p < 0.001:
        PVals.append('< 0.001')

# save proportions to file
newfile = open('ProportionHomologs.txt', 'w')    
# create contingency table table
newfile.write('Table 1. Proportion of genes with chimp homologs\n')
newfile.write('\t'.join(['Genes', 'Homologs', 'No-homologs', 'Proportion', 'P']) + '\n')
GeneCats = ['Non-Overlapping', 'Nested', 'Internal', 'External', 'Piggyback', 'Convergent', 'Divergent']
for i in range(len(GeneCats)):
    newfile.write('\t'.join([GeneCats[i], str(GeneCounts[i][0]), str(GeneCounts[i][1]), str(round(GeneCounts[i][0]/sum(GeneCounts[i]), 4)), str(PVals[i])]) + '\n')
newfile.close()          
          
