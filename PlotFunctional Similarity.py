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









for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# make pairs of overlapping genes
HumanPairs = AllPairs[:2]
SisterPairs = AllPairs[2:4]
OutGroupPairs = AllPairs[4:]

# make list with sets of non-overlapping genes
NonOverlappingSets = []
for i in range(3):
    j = i * 2
    # make a set of non-overlapping gene
    nonoverlap = MakeNonOverlappingGeneSet(AllOverlap[j], AllCoordinates[i])
    NonOverlappingSets.append(nonoverlap)    

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

# get 1:1 orthologs between human and sister-species
# get 1:1 orthologs between human, sister-species and outgroup {human:[sistersp,outgroup]}
if SisterSp == 'chimp':
    OrthoPairs = MatchOrthologPairs('HumanChimpOrthologs.txt')
    OrthoTrios = MatchOrthologTrios('HumanChimpGorillaOrthologs.txt')
elif SisterSp == 'mouse':
    OrthoPairs = MatchOrthologPairs('HumanMouseOrthologs.txt')
    OrthoTrios = MatchOrthologTrios('HumanMouseDogOrthologs.txt')

# reverse dict with human and sister-species orthologs
SisterOrthos = {}
for gene in OrthoPairs:
    SisterOrthos[OrthoPairs[gene]] = gene
# make a dict of ortho trios with sister-species genes as key
SisterOrthoTrios = {}
for gene in OrthoTrios:
    sistergene, outgroupgene = OrthoTrios[gene][0], OrthoTrios[gene][1]
    SisterOrthoTrios[sistergene] = [gene, outgroupgene]


# infer young and old nesting events 
HumanOld, HumanYoung = InferYoungOldNestingEvents(OrthoPairs, OrthoTrios, SisterPairs[1], OutGroupPairs[1], HumanPairs[1])
SisterSpOld, SisterSpYoung = InferYoungOldNestingEvents(SisterOrthos, SisterOrthoTrios, HumanPairs[1], OutGroupPairs[1], SisterPairs[1])

# do some QC   
a = [set([OrthoPairs[pair[0]], OrthoPairs[pair[1]]]) for pair in HumanOld]
b = [set(pair) for pair in SisterSpOld]
assert len(a) == len(b)
for i in a:
    assert i in b    

# get expression profiles of human and sister-species
if SisterSp == 'chimp':
    HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
    SisterSpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
elif SisterSp == 'mouse':
    HumanExpression = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
    SisterSpExpression = ParseExpressionFile('Mouse_Median_Normalized_FPKM.txt')
    # match expression profiles between mouse and human 
    HumanExpression = MatchHumanToMouseExpressionProfiles(HumanExpression)

# remove genes without expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
SisterSpExpression = RemoveGenesLackingExpression(SisterSpExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
SisterSpExpression = TransformRelativeExpression(SisterSpExpression)

# compute expression specificity
HumanSpecificity = ExpressionSpecificity(HumanExpression)
SisterSpSpecificity = ExpressionSpecificity(SisterSpExpression)

HumanInferredPairs = [HumanOld, HumanYoung]
SisterSpInferredPairs = [SisterSpOld, SisterSpYoung]



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



