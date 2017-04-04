# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 21:24:45 2017

@author: Richard
"""

# use this script to plot differences in functional similarities between overlapping genes

# usage PlotFunctionalSimilarity
#- [all/molfunc/celcomp/biolproc]: use all GO terns, molecular function, cellular compartments or biological processes

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



#get option from command
GOClass = sys.argv[1]
assert GOClass in ['all', 'molfunc', 'celcomp', 'biolproc']



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

# check option to see GO ids are filtered or not
if GOClass == 'molfunc':
    GOSubSet = FilterGOTerms('GOMolecularFunction.txt')
elif GOClass == 'celcomp':
    GOSubSet = FilterGOTerms('GOCellularCompartments.txt')
elif GOClass == 'biolproc':
    GOSubSet = FilterGOTerms('GOBiologicalProcesses.txt')
# keep only terms for specific GO class
if GOClass != 'all':
    for gene in GeneOntology:
         GeneOntology[gene] = set(filter(lambda x: x in GOSubSet, GeneOntology[gene]))

# remove genes if gene has no GO term
to_remove = [gene for gene in GeneOntology if GeneOntology[gene] == 0]
if len(to_remove) != 0:
    for gene in to_remove:
        del GeneOntology[gene]

# create a list of lists with functional similarities between genes of the same pair
FunctionalSimilarity = []
# include only the 4 overlapping types
for i in range(1, len(OverlappingPairs)):
    # compute jaccard similarity index between each pair
    GOOverlap = []
    for pair in OverlappingPairs[i]:
        if pair[0] in GeneOntology and pair[1] in GeneOntology:
            JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
            GOOverlap.append(JI)
    FunctionalSimilarity.append(GOOverlap)
 
# generate random pairs of non-overlapping genes as a control
# make a list of non-overlapping genes convert set of non-overlapping genes to list
NonOverlappingGenes = list(NonOverlappingGenes)
ControlPairs = []
# draw 10000 random pairs
replicates = 10000
while replicates != 0:
    # pick 2 random genes
    j = random.randint(0, len(NonOverlappingGenes) -1)
    k = random.randint(0, len(NonOverlappingGenes) -1)
    gene1, gene2 = NonOverlappingGenes[j], NonOverlappingGenes[k]
    if gene1 in GeneOntology and gene2 in GeneOntology:
        JI = JaccardIndex(GeneOntology[gene1], GeneOntology[gene2])
        ControlPairs.append(JI)
        replicates -= 1
  
# insert pairs of JI in list
FunctionalSimilarity.insert(0, ControlPairs)

DataSets = ['Control', 'Nested', 'PiggyBack', 'Convergent', 'Divergent']

for i in range(1, len(FunctionalSimilarity)):
    P = PermutationResampling(FunctionalSimilarity[0], FunctionalSimilarity[i], 1000, statistic = np.mean)
    print(DataSets[0], DataSets[i], len(FunctionalSimilarity[0]), len(FunctionalSimilarity[i]), np.mean(FunctionalSimilarity[0]), np.mean(FunctionalSimilarity[i]), P)


###############################

    

# make sets of host and nested nested genes
OverlapSets = MakeFullPartialOverlapGeneSet(Overlap[0])
HumanExpression = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
# compute expression specificity
HumanSpecificity = ExpressionSpecificity(HumanExpression)
# generate a dict to draw genes in human    
HumanRandomGenes = GenerateAllUnNestedGenes(OverlapSets, OrderedGenes, HumanExpression)
    
    
NestedPairs = copy.deepcopy(OverlappingPairs[1])
# remove human pairs if genes are not expressed
to_remove = [pair for pair in NestedPairs if pair[0] not in HumanExpression or pair[1] not in HumanExpression]
for pair in to_remove:
    NestedPairs.remove(pair)
print('expressed nested', len(NestedPairs), len(OverlappingPairs[1]))
  



  
# make a list of control un-nested pairs in sister species
HumanControlPairs = []    
for pair in NestedPairs:
    # make a list of matching gene pairs (orientation, chromosome, distance)
    PairPool = GenerateMatchingPoolPairs(pair, HumanRandomGenes, GeneCoord, 2000)
    # draw a matching gene pair at random
    i = random.randint(0, len(PairPool) -1)
    assert PairPool[i][0] not in OverlapSets and PairPool[i][1] not in OverlapSets
    HumanControlPairs.append(PairPool[i])


NestedJI = []
for pair in NestedPairs:
    if pair[0] in GeneOntology and pair[1] in GeneOntology:
        JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
        NestedJI.append(JI)
ControlJI = []
for pair in HumanControlPairs:
    if pair[0] in GeneOntology and pair[1] in GeneOntology:
        JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
        ControlJI.append(JI)    

P = PermutationResampling(NestedJI, ControlJI, 1000, statistic = np.mean)
print(len(NestedJI), len(ControlJI), np.mean(NestedJI), np.mean(ControlJI), P)