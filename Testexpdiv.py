# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:10:15 2016

@author: RJovelin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""


# use this script to plot expression divergence between external their
# un-nested orthologs and between internal and their un-nested orthologs

# usage python3 PlotExpDivergYoungOld.py 


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

# load dictionaries of host and nested genes 
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('ChimpHostNestedGenes.json') as chimp_json_data:
    ChimpHostGenes = json.load(chimp_json_data)
with open('GorillaHostNestedGenes.json') as gorilla_json_data:
    GorillaHostGenes = json.load(gorilla_json_data)

# load dictionaries of host and contained genes
with open('HumanContainedGenes.json') as human_json_data:
    HumanContainedGenes = json.load(human_json_data)
with open('ChimpContainedGenes.json') as chimp_json_data:
    ChimpContainedGenes = json.load(chimp_json_data)
with open('GorillaContainedGenes.json') as gorilla_json_data:
    GorillaContainedGenes = json.load(gorilla_json_data)

# load dictionaries with overlapping genes
with open('HumanOverlappingGenes.json') as human_json_data:
    HumanOverlappingGenes = json.load(human_json_data)
with open('ChimpOverlappingGenes.json') as chimp_json_data:
    ChimpOverlappingGenes = json.load(chimp_json_data)
with open('GorillaOverlappingGenes.json') as gorilla_json_data:
    GorillaOverlappingGenes = json.load(gorilla_json_data)


# make a list of host:nested genes dictionaries
HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes]

# get human GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'

# get the coordinates of human genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
HumanGeneChromoCoord = ChromoGenesCoord(HsaGFF)
# map each gene to its mRNA transcripts
HumanMapGeneTranscript = GeneToTranscripts(HsaGFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
HumanGeneChromoCoord = FilterOutGenesWithoutValidTranscript(HumanGeneChromoCoord, HumanMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
HumanGeneCoord = FromChromoCoordToGeneCoord(HumanGeneChromoCoord)
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
HumanOrderedGenes = OrderGenesAlongChromo(HumanGeneChromoCoord)
# get expression profile of the species genes
HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
# remove genes wuthout expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
# get expression specificity
HumanExpSpecificity = ExpressionSpecificity(HumanExpression)
# make a list of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HostGenes[0])
# remove gene pairs with genes lacking expression
HumanHostNestedPairs = FilterGenePairsWithoutExpression(HumanHostNestedPairs, HumanExpression)


# get expression profile of the species genes
ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
# remove genes wuthout expression
ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
# get relative expression
ChimpExpression = TransformRelativeExpression(ChimpExpression)
# make a list of host-nested gene pairs
ChimpHostNestedPairs = GetHostNestedPairs(HostGenes[1])
# remove gene pairs with genes lacking expression
ChimpHostNestedPairs = FilterGenePairsWithoutExpression(ChimpHostNestedPairs, ChimpExpression)


# get expression profile of the species genes
GorillaExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'gorilla')
# remove genes wuthout expression
GorillaExpression = RemoveGenesLackingExpression(GorillaExpression)
# get relative expression
GorillaExpression = TransformRelativeExpression(GorillaExpression)
# make a list of host-nested gene pairs
GorillaHostNestedPairs = GetHostNestedPairs(HostGenes[2])
# remove gene pairs with genes lacking expression
GorillaHostNestedPairs = FilterGenePairsWithoutExpression(GorillaHostNestedPairs, GorillaExpression)

# get chimp and gorilla orthologs of human genes
HumanOrthologs = ParseOrthologFile('HumanChimpGorillaOrthologs.txt')
# create a dict with chimp genes as key {chimp_gene: [human_orto, gorilla_ortho]}
ChimpOrthologs = {}
for gene in HumanOrthologs:
    ChimpOrthologs[HumanOrthologs[gene][0]] = [gene, HumanOrthologs[gene][1]]


# make lists of overlapping gene pairs
HumanOverlappingPairs = GetHostNestedPairs(HumanOverlappingGenes)
ChimpOverlappingPairs = GetHostNestedPairs(ChimpOverlappingGenes)
GorillaOverlappingPairs = GetHostNestedPairs(GorillaOverlappingGenes)
# make lists of contained gene pairs
HumanContainedPairs = GetHostNestedPairs(HumanContainedGenes)
ChimpContainedPairs = GetHostNestedPairs(ChimpContainedGenes)
GorillaContainedPairs = GetHostNestedPairs(GorillaContainedGenes)


# make a set of human host and nested genes
HumanNestedGeneSet = MakeFullPartialOverlapGeneSet(HostGenes[0])
# make a set of chimp host and nested genes
ChimpNestedGeneSet = MakeFullPartialOverlapGeneSet(HostGenes[1])

# make a set of human genes excluding host and nested genes
HumanNonNestedGeneSet = MakeNonOverlappingGeneSet(HostGenes[0], HumanGeneCoord)


# create a dict to randomly draw genes
ToDrawFrom = GenerateAllUnNestedGenes(HumanNonNestedGeneSet, HumanOrderedGenes)



# infer young and old nesting events in human and chimp
HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, ChimpHostNestedPairs, GorillaHostNestedPairs, HumanHostNestedPairs)
ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HumanHostNestedPairs, GorillaHostNestedPairs, ChimpHostNestedPairs)


# compare expression divergence between young nested genes and their un-nested orthologs
# compare expression divergence between young host genes and their un-nested orthologs


ExternalReplicates, InternalReplicates = [], []


for num in range(100):


    # populate lists with young genes and their un-nested orthologs 
    YoungInternal, YoungExternal = [], []
    # get pairs of human young external and internal genes and their un-nested chimp orthologs
    for pair in HumanYoung:
        # get the ortholog of the host and nested genes
        extortho, internortho = HumanOrthologs[pair[0]][0], HumanOrthologs[pair[1]][0]
        # check that these genes are not nested
        if extortho not in ChimpNestedGeneSet and internortho not in ChimpNestedGeneSet:
            if extortho in ChimpExpression:
                assert pair[0] in HumanExpression
                YoungExternal.append([pair[0], extortho])
            if internortho in ChimpExpression:
                assert pair[1] in HumanExpression
                YoungInternal.append([pair[1], internortho])
    # get pairs of chimp young external and internal genes and their un-nested human orthologs      
    for pair in ChimpYoung:
        # get the ortholog of the host and nested genes
        extortho, internortho = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
         # check that these genes are not nested
        if extortho not in HumanNestedGeneSet and internortho not in HumanNestedGeneSet:
            if extortho in HumanExpression:
                assert pair[0] in ChimpExpression
                YoungExternal.append([extortho, pair[0]])
            if internortho in HumanExpression:
                assert pair[1] in ChimpExpression
                YoungInternal.append([internortho, pair[1]])
    

## compute divergence between young nested pairs and their un-nested orthologs 
#YoungInternalDiv = ComputeExpressionDivergenceOrthologs(YoungInternal, HumanExpression, ChimpExpression)
#YoungExternalDiv = ComputeExpressionDivergenceOrthologs(YoungExternal, HumanExpression, ChimpExpression)
#print(len(YoungInternalDiv), np.mean(YoungInternalDiv), len(YoungExternalDiv), np.mean(YoungExternalDiv), stats.ranksums(YoungInternalDiv, YoungExternalDiv)[1])



    

    

    # create pairs of random internal-like and external-like genes in human and their un-nested orthologs in chimp
    InternalLike, ExternalLike = [], []

    # create a list of pairs to remove when genes have no match
    to_remove = []

    # for each human gene, match a random un-nested gene with similar characterisitics
    for pair in YoungInternal:
        # get the chromosome of the human gene
        chromo = HumanGeneCoord[pair[0]][0]
        # create a list of genes corresponding to all possible genes on chromo to draw from
        PossibleGenes = list(ToDrawFrom[chromo].keys())
        # draw a random gene on that chromo with matching characteristics
        NotFound = True
        while len(PossibleGenes) != 0 and NotFound == True:
            i = random.choice(PossibleGenes)
            gene = ToDrawFrom[chromo][i]
            # assert gene not nested
            assert gene not in HumanNestedGeneSet
            # match gene by expression specificity (+- 0.01)
            if gene in HumanExpSpecificity:
                assert gene in HumanExpression
                # match by expression specificity
                if HumanExpSpecificity[pair[0]] - 0.02 <= HumanExpSpecificity[gene] <= HumanExpSpecificity[pair[0]] + 0.02:
                    # match by tissue with highest expression
                    if HumanExpression[gene].index(max(HumanExpression[gene])) == HumanExpression[pair[0]].index(max(HumanExpression[pair[0]])):
                        # check that matching gene has a un-nested chimp ortholog
                        if gene in HumanOrthologs:
                            # check that ortholog is not nested and is expressed
                            if HumanOrthologs[gene][0] not in ChimpNestedGeneSet and HumanOrthologs[gene][0] in ChimpExpression:
                                # found internal-like gene, populate list
                                InternalLike.append([gene, HumanOrthologs[gene][0]])
                                # update boolean
                                NotFound = False
                            else:
                                PossibleGenes.remove(i)
                        else:
                            PossibleGenes.remove(i)
                    else:
                        PossibleGenes.remove(i)
                else:
                    PossibleGenes.remove(i)
            else:
                PossibleGenes.remove(i)
        if len(PossibleGenes) == 0:
            to_remove.append(pair)                
    print(len(YoungInternal), len(to_remove), len(InternalLike))

    if len(to_remove) != 0:
        for pair in to_remove:
            YoungInternal.remove(pair)

    # create a list of pairs to remove when genes have no match
    to_remove = []

    # for each human gene, match a random un-nested gene with similar characterisitics
    for pair in YoungExternal:
        # get the chromosome of the human gene
        chromo = HumanGeneCoord[pair[0]][0]
        # create a list of genes corresponding to all possible genes on chromo to draw from
        PossibleGenes = list(ToDrawFrom[chromo].keys())
        # draw a random gene on that chromo with matching characteristics
        NotFound = True
        while len(PossibleGenes) != 0 and NotFound == True:
            i = random.choice(PossibleGenes)
            gene = ToDrawFrom[chromo][i]
            # assert gene not nested
            assert gene not in HumanNestedGeneSet
            # match gene by expression specificity (+- 0.01)
            if gene in HumanExpSpecificity:
                assert gene in HumanExpression
                # match by expression specificity
                if HumanExpSpecificity[pair[0]] - 0.02 <= HumanExpSpecificity[gene] <= HumanExpSpecificity[pair[0]] + 0.02:
                    # match by tissue with highest expression
                    if HumanExpression[gene].index(max(HumanExpression[gene])) == HumanExpression[pair[0]].index(max(HumanExpression[pair[0]])):
                        # check that matching gene has a un-nested chimp ortholog
                        if gene in HumanOrthologs:
                            # check that ortholog is not nested and is expressed
                            if HumanOrthologs[gene][0] not in ChimpNestedGeneSet and HumanOrthologs[gene][0] in ChimpExpression:
                                # found internal-like gene, populate list
                                ExternalLike.append([gene, HumanOrthologs[gene][0]])
                                # update boolean
                                NotFound = False
                            else:
                                PossibleGenes.remove(i)
                        else:
                            PossibleGenes.remove(i)
                    else:
                        PossibleGenes.remove(i)
                else:
                    PossibleGenes.remove(i)
            else:
                PossibleGenes.remove(i)
        if len(PossibleGenes) == 0:
            to_remove.append(pair)                
    print(len(YoungExternal), len(to_remove), len(ExternalLike))

    if len(to_remove) != 0:
        for pair in to_remove:
            YoungExternal.remove(pair)

    # compute divergence between young nested pairs and their un-nested orthologs 
    YoungInternalDiv = ComputeExpressionDivergenceOrthologs(YoungInternal, HumanExpression, ChimpExpression)
    YoungExternalDiv = ComputeExpressionDivergenceOrthologs(YoungExternal, HumanExpression, ChimpExpression)
    print(len(YoungInternalDiv), np.mean(YoungInternalDiv), len(YoungExternalDiv), np.mean(YoungExternalDiv), stats.ranksums(YoungInternalDiv, YoungExternalDiv)[1])
    # compute divergence between external like and their un-nested orthologs
    ExternalDiv = ComputeExpressionDivergenceOrthologs(ExternalLike, HumanExpression, ChimpExpression)
    InternalDiv = ComputeExpressionDivergenceOrthologs(InternalLike, HumanExpression, ChimpExpression)
    # compare expression divergence between young external and external-like and between young internal and internal-like
    print(len(YoungExternalDiv), np.mean(YoungExternalDiv), len(ExternalDiv), np.mean(ExternalDiv), stats.ranksums(YoungExternalDiv, ExternalDiv)[1])
    print(len(YoungInternalDiv), np.mean(YoungInternalDiv), len(InternalDiv), np.mean(InternalDiv), stats.ranksums(YoungInternalDiv, InternalDiv)[1])
    ExternalReplicates.append(stats.ranksums(YoungExternalDiv, ExternalDiv)[1])
    InternalReplicates.append(stats.ranksums(YoungInternalDiv, InternalDiv)[1])
    

print(min(ExternalReplicates), max(ExternalReplicates))
total = 0
for i in ExternalReplicates:
    if i < 0.05:
        total += 1
print('N replicates with p 0.05 young external', total)
print(min(InternalReplicates), max(InternalReplicates))
total = 0
for i in InternalReplicates:
    if i < 0.05:
        total += 1
print('N replicates with p 0.05 young internal', total)




