# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:51:38 2017

@author: RJovelin
"""

# use this script to make a table with correlation between expression divergence and functional similarity among gene pairs

# usage python3 MakeTableCorrelationGOExpDiv.py 

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
L = len(GeneOntology)

# get the expression profile of each gene
Expression = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without expression
Expression = RemoveGenesLackingExpression(Expression)
# get relative expression
Expression = TransformRelativeExpression(Expression)

# make a dict of correlation coefficients for each GO class {all: [baseline, nsted, pbk, conv, div], ...'biolproc':[baseline, nsted, pbk, conv, div]}
Correlations = {}
for GOClass in ['all', 'molfunc', 'celcomp', 'biolproc']:
    # make a set of GO ontology terms
    if GOClass == 'molfunc':
        GOSubSet = FilterGOTerms('GOMolecularFunction.txt')
    elif GOClass == 'celcomp':
        GOSubSet = FilterGOTerms('GOCellularCompartments.txt')
    elif GOClass == 'biolproc':
        GOSubSet = FilterGOTerms('GOBiologicalProcesses.txt')
    # copy GeneOntology dict
    GeneGOTerms = copy.deepcopy(GeneOntology)
    assert GeneGOTerms == GeneOntology
    assert L == len(GeneGOTerms) 
    # keep only terms for specific GO class
    if GOClass != 'all':
        for gene in GeneGOTerms:
            GeneGOTerms[gene] = set(filter(lambda x: x in GOSubSet, GeneGOTerms[gene]))
    # remove genes if gene has no GO term
    to_remove = [gene for gene in GeneGOTerms if len(GeneGOTerms[gene]) == 0]
    if len(to_remove) != 0:
        for gene in to_remove:
            del GeneGOTerms[gene]
    # make a list of non-overlapping genes that have GO and that are expressed            
    NonOverlapping = [gene for gene in NonOverlappingGenes if gene in GeneGOTerms and gene in Expression]
        
    # generate random pairs of non-overlapping genes and compute JI and expression divergence
    # store JI and expdiv in parrallel lists [[JI], [exodiv]]
    BaseLine = [[], []]
    # draw 10000 random pairs
    replicates = 10000
    while replicates != 0:
        # pick 2 random genes
        j = random.randint(0, len(NonOverlapping) -1)
        k = random.randint(0, len(NonOverlapping) -1)
        gene1, gene2 = NonOverlapping[j], NonOverlapping[k]
        # compute functional similarity index
        JI = JaccardIndex(GeneGOTerms[gene1], GeneGOTerms[gene2])
        # compute expression divergence        
        expdiv = EuclidianDistance(Expression[gene1], Expression[gene2])
        BaseLine[0].append(JI)
        BaseLine[1].append(expdiv)    
        replicates -= 1
    
    # initialize dict with correlation between expression divergence and functional similarity for the random genes to list
    a, b = stats.spearmanr(BaseLine[0], BaseLine[1])
    Correlations[GOClass] = [[a, b]]      
    
    # include only the 4 overlapping types    
    for i in range(1, len(OverlappingPairs)):
        # make a list of list of overlapping pairs that have GO terms and that are expressed
        OvlPairs = [pair for pair in OverlappingPairs[i] if (pair[0] in GeneGOTerms and pair[1] in GeneGOTerms) and (pair[0] in Expression and pair[1] in Expression)]
        # compute jaccard similarity index  and expression divergence between each pair    
        GOOverlap, ExpDiv = [], []
        for pair in OvlPairs:
            JI = JaccardIndex(GeneGOTerms[pair[0]], GeneGOTerms[pair[1]])
            expdiv = EuclidianDistance(Expression[pair[0]], Expression[pair[1]])
            GOOverlap.append(JI)
            ExpDiv.append(expdiv)
        # compute correlation between functional similarity and expression divergence
        a, b = stats.spearmanr(GOOverlap, ExpDiv)
        Correlations[GOClass].append([a, b])        
        
# open file to store the results of the correlation in table
newfile = open('CorrelationsGOExpDiv.txt', 'w')

# create table
newfile.write('Table 1. Correlation between expression divergence and Gene Ontology similarity\n')
newfile.write('\t'.join(['GO Class', 'Random', 'Nested', 'PiggyBack', 'Convergent', 'DivergentGenes']) + '\n')

GOTypes = ['all', 'molfunc', 'celcomp', 'biolproc']
GONames = ['All', 'Molecular Function', 'Cellular Compartments', 'Biological Processes']
for i in range(len(GOTypes)):
    line = [GONames[i]]
    for correl in Correlations[GOTypes[i]]:
        line.append(str(correl[0]) + ' (' + str(correl[1]) + ')')
    line = '\t'.join(line)
    newfile.write(line + '\n')
# close file after writing
newfile.close()

