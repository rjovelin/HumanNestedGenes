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


# load dictionaries of overlapping genes
jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json',
             'GorillaOverlappingGenes.json', 'GorillaNestedGenes.json']
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
GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3', 'Gorilla_gorilla.gorGor3.1.86.gff3']

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


# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# make pairs of overlapping genes
HumanPairs = AllPairs[:2]
ChimpPairs = AllPairs[2:4]
GorillaPairs = AllPairs[4:]

# make list with sets of non-overlapping genes
NonOverlapping = []
for i in range(3):
    # make a set of overlapping genes
    if i == 0:
        overlap = MakeFullPartialOverlapGeneSet(HumanPairs[0])
    elif i == 1:
        overlap = MakeFullPartialOverlapGeneSet(ChimpPairs[0])
    elif i == 2:
        overlap = MakeFullPartialOverlapGeneSet(GorillaPairs[0])
    # make a set of non-overlapping gene
    nonoverlap = MakeNonOverlappingGeneSet(overlap, AllCoordinates[i])    

# get 1:1 orthologs between human, chimp and gorilla {human:[chimp,gorilla]}
Orthos = MatchOrthologTrios('HumanChimpGorillaOrthologs.txt')

print(len(HumanPairs[0]), len(HumanPairs[1]), len(ChimpPairs[0]), len(ChimpPairs[1]), len(GorillaPairs[0]), len(GorillaPairs[1]))

# make sets of chimp and gorilla genes  with orthologs
chimporthos, gorillaorthos = set(), set()
for gene in Orthos:
    chimporthos.add(Orthos[gene][0])
    gorillaorthos.add(Orthos[gene][1])

# remove gene pairs lacking orthos
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos or pair[1] not in orthos:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
# remove chimp pairs lacking orthologs
for i in range(len(ChimpPairs)):
    to_remove = []
    for pair in ChimpPairs[i]:
        if pair[0] not in chimporthos or pair[1] not in chimporthos:
            to_remove.append(pair)
    for pair in to_remove:
        ChimpPairs[i].remove(pair)
# remove gorilla pairs lacking orthologs
for i in range(len(GorillaPairs)):
    to_remove = []
    for pair in GorillaPairs[i]:
        if pair[0] not in gorillaorthos or pair[1] not in gorillaorthos:
            to_remove.append(pair)
    for pair in to_remove:
        Gorillapairs[i].remove(pair)

print(len(HumanPairs[0]), len(HumanPairs[1]), len(ChimpPairs[0]), len(ChimpPairs[1]), len(GorillaPairs[0]), len(GorillaPairs[1]))













## get expression profile in human and chimp
#HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
#ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
#
## remove genes without expression
#HumanExpression = RemoveGenesLackingExpression(HumanExpression)
#ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
#
## get relative expression
#HumanExpression = TransformRelativeExpression(HumanExpression)
#ChimpExpression = TransformRelativeExpression(ChimpExpression)
#
## remove gene pairs with genes lacking expression
#for i in range(len(HumanPairs)):
#    HumanPairs[i] = FilterGenePairsWithoutExpression(HumanPairs[i], HumanExpression)
#for i in range(len(ChimpPairs)):
#    ChimpPairs[i] = FilterGenePairsWithoutExpression(ChimpPairs[i], ChimpExpression)
#
#
#print(len(HumanPairs[0]), len(HumanPairs[1]), len(ChimpPairs[0]), len(ChimpPairs[1]), len(GorillaPairs[0]), len(GorillaPairs[1]))
#
#
#
#
#
## create a dict with chimp genes as key {chimp_gene: [human_orto, gorilla_ortho]}
#ChimpOrthologs = {}
#for gene in HumanOrthologs:
#    ChimpOrthologs[HumanOrthologs[gene][0]] = [gene, HumanOrthologs[gene][1]]
#
#
#
## make a set of human host and nested genes
#HumanNestedGeneSet = MakeFullPartialOverlapGeneSet(HostGenes[0])
## make a set of chimp host and nested genes
#ChimpNestedGeneSet = MakeFullPartialOverlapGeneSet(HostGenes[1])
#
## make a set of human genes excluding host and nested genes
#HumanNonNestedGeneSet = MakeNonOverlappingGeneSet(HostGenes[0], HumanGeneCoord)
#
#
#
###########################
#
#
## use this function to sort young and ancestral nesting events
#def InferYoungOldNestingEvents(FirstSpOrthologs, SecondSpOverlapPairs, OutGroupOverlapPairs, FirstSpHostNestedPairs):
#    '''
#    (dict, list, list, list) -> (list, list)
#    Take a dictionary with ortholog gene trios, the lists of overlapping gene
#    pairs in the sister species and in the outgroup and the list of host:nested
#    pairs in the focal species and return a tuple with list of host: nested pairs
#    that are infered to be old and young (before the divergence of the 2 species or after)
#    '''   
#    # create lists of sets of gene pairs to remove the order between genes
#    SecondOverlap, OutGroupOverlap = [], []
#    for pair in SecondSpOverlapPairs:
#        SecondOverlap.append(set(pair))
#    for pair in OutGroupOverlapPairs:
#        OutGroupOverlap.append(set(pair))
#    # create lists of host-nested gene pairs that are old (present in second
#    # species, or second species and outgroup) or young (not found in outgroup
#    # and not found in the second species, ie species-specific)    
#    old, young = [], []
#    # loop over the gene pairs    
#    for pair in FirstSpHostNestedPairs:
#        # check that both genes have orthologs
#        if pair[0] in FirstSpOrthologs and pair[1] in FirstSpOrthologs:
#            # check if gene pair is overlapping in second species and outgroup
#            # get ortholog of host in 2nd species
#            hostorthosp2 = FirstSpOrthologs[pair[0]][0]
#            # get ortholog of nested in 2nd species
#            nestedorthosp2 = FirstSpOrthologs[pair[1]][0]
#            # get ortholog of host in outgroup
#            hostorthoout = FirstSpOrthologs[pair[0]][1]
#            # get ortholog of nested gene in outgroup
#            nestedorthoout = FirstSpOrthologs[pair[1]][1]
#            # check if orthologs are overlapping in 2nd species
#            if {hostorthosp2, nestedorthosp2} not in SecondOverlap and {hostorthoout, nestedorthoout} not in OutGroupOverlap:
#                # nested is species-specific
#                young.append(pair)
#            elif {hostorthosp2, nestedorthosp2} in SecondOverlap and {hostorthoout, nestedorthoout} not in OutGroupOverlap:
#                # gene pair is overlapping in 2nd species, suffcient to be called old nested pair
#                old.append(pair)
#            elif {hostorthosp2, nestedorthosp2} in SecondOverlap and {hostorthoout, nestedorthoout} in OutGroupOverlap:
#                # gene pair is overlapping in 2nd species and outgroup, infer old nested pair
#                old.append(pair)
#
#    return old, young
#
#
#
#
#
##############################
#
#
#
## infer young and old nesting events in human and chimp
#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, ChimpHostNestedPairs, GorillaHostNestedPairs, HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HumanHostNestedPairs, GorillaHostNestedPairs, ChimpHostNestedPairs)
#
#
## compare expression divergence between young nested genes and their un-nested orthologs
## compare expression divergence between young host genes and their un-nested orthologs
#ExternalReplicates, InternalReplicates = [], []
#
#
#Nreplicates = 500
#
#
#for num in range(Nreplicates):
#
#
#    # populate lists with young genes and their un-nested orthologs 
#    YoungInternal, YoungExternal = [], []
#    # get pairs of human young external and internal genes and their un-nested chimp orthologs
#    for pair in HumanYoung:
#        # get the ortholog of the host and nested genes
#        extortho, internortho = HumanOrthologs[pair[0]][0], HumanOrthologs[pair[1]][0]
#        # check that these genes are not nested
#        if extortho not in ChimpNestedGeneSet and internortho not in ChimpNestedGeneSet:
#            if extortho in ChimpExpression:
#                assert pair[0] in HumanExpression
#                YoungExternal.append([pair[0], extortho])
#            if internortho in ChimpExpression:
#                assert pair[1] in HumanExpression
#                YoungInternal.append([pair[1], internortho])
#    # get pairs of chimp young external and internal genes and their un-nested human orthologs      
#    for pair in ChimpYoung:
#        # get the ortholog of the host and nested genes
#        extortho, internortho = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
#         # check that these genes are not nested
#        if extortho not in HumanNestedGeneSet and internortho not in HumanNestedGeneSet:
#            if extortho in HumanExpression:
#                assert pair[0] in ChimpExpression
#                YoungExternal.append([extortho, pair[0]])
#            if internortho in HumanExpression:
#                assert pair[1] in ChimpExpression
#                YoungInternal.append([internortho, pair[1]])
#    
#
### compute divergence between young nested pairs and their un-nested orthologs 
##YoungInternalDiv = ComputeExpressionDivergenceOrthologs(YoungInternal, HumanExpression, ChimpExpression)
##YoungExternalDiv = ComputeExpressionDivergenceOrthologs(YoungExternal, HumanExpression, ChimpExpression)
##print(len(YoungInternalDiv), np.mean(YoungInternalDiv), len(YoungExternalDiv), np.mean(YoungExternalDiv), stats.ranksums(YoungInternalDiv, YoungExternalDiv)[1])
#
#
#
#    
#
#    
#
#    # create pairs of random internal-like and external-like genes in human and their un-nested orthologs in chimp
#    InternalLike, ExternalLike = [], []
#
#    # create a list of pairs to remove when genes have no match
#    to_remove = []
#
#    # for each human gene, match a random un-nested gene with similar characterisitics
#    for pair in YoungInternal:
#        # get the chromosome of the human gene
#        chromo = HumanGeneCoord[pair[0]][0]
#        # create a list of genes corresponding to all possible genes on chromo to draw from
#        PossibleGenes = list(ToDrawFrom[chromo].keys())
#        # draw a random gene on that chromo with matching characteristics
#        NotFound = True
#        while len(PossibleGenes) != 0 and NotFound == True:
#            i = random.choice(PossibleGenes)
#            gene = ToDrawFrom[chromo][i]
#            # assert gene not nested
#            assert gene not in HumanNestedGeneSet
#            # match gene by expression specificity (+- 0.01)
#            if gene in HumanExpSpecificity:
#                assert gene in HumanExpression
#                # match by expression specificity
#                if HumanExpSpecificity[pair[0]] - 0.02 <= HumanExpSpecificity[gene] <= HumanExpSpecificity[pair[0]] + 0.02:
#                    # match by tissue with highest expression
#                    if HumanExpression[gene].index(max(HumanExpression[gene])) == HumanExpression[pair[0]].index(max(HumanExpression[pair[0]])):
#                        # check that matching gene has a un-nested chimp ortholog
#                        if gene in HumanOrthologs:
#                            # check that ortholog is not nested and is expressed
#                            if HumanOrthologs[gene][0] not in ChimpNestedGeneSet and HumanOrthologs[gene][0] in ChimpExpression:
#                                # found internal-like gene, populate list
#                                InternalLike.append([gene, HumanOrthologs[gene][0]])
#                                # update boolean
#                                NotFound = False
#                            else:
#                                PossibleGenes.remove(i)
#                        else:
#                            PossibleGenes.remove(i)
#                    else:
#                        PossibleGenes.remove(i)
#                else:
#                    PossibleGenes.remove(i)
#            else:
#                PossibleGenes.remove(i)
#        if len(PossibleGenes) == 0:
#            to_remove.append(pair)                
#    print(len(YoungInternal), len(to_remove), len(InternalLike))
#
#    if len(to_remove) != 0:
#        for pair in to_remove:
#            YoungInternal.remove(pair)
#
#    # create a list of pairs to remove when genes have no match
#    to_remove = []
#
#    # for each human gene, match a random un-nested gene with similar characterisitics
#    for pair in YoungExternal:
#        # get the chromosome of the human gene
#        chromo = HumanGeneCoord[pair[0]][0]
#        # create a list of genes corresponding to all possible genes on chromo to draw from
#        PossibleGenes = list(ToDrawFrom[chromo].keys())
#        # draw a random gene on that chromo with matching characteristics
#        NotFound = True
#        while len(PossibleGenes) != 0 and NotFound == True:
#            i = random.choice(PossibleGenes)
#            gene = ToDrawFrom[chromo][i]
#            # assert gene not nested
#            assert gene not in HumanNestedGeneSet
#            # match gene by expression specificity (+- 0.01)
#            if gene in HumanExpSpecificity:
#                assert gene in HumanExpression
#                # match by expression specificity
#                if HumanExpSpecificity[pair[0]] - 0.02 <= HumanExpSpecificity[gene] <= HumanExpSpecificity[pair[0]] + 0.02:
#                    # match by tissue with highest expression
#                    if HumanExpression[gene].index(max(HumanExpression[gene])) == HumanExpression[pair[0]].index(max(HumanExpression[pair[0]])):
#                        # check that matching gene has a un-nested chimp ortholog
#                        if gene in HumanOrthologs:
#                            # check that ortholog is not nested and is expressed
#                            if HumanOrthologs[gene][0] not in ChimpNestedGeneSet and HumanOrthologs[gene][0] in ChimpExpression:
#                                # found internal-like gene, populate list
#                                ExternalLike.append([gene, HumanOrthologs[gene][0]])
#                                # update boolean
#                                NotFound = False
#                            else:
#                                PossibleGenes.remove(i)
#                        else:
#                            PossibleGenes.remove(i)
#                    else:
#                        PossibleGenes.remove(i)
#                else:
#                    PossibleGenes.remove(i)
#            else:
#                PossibleGenes.remove(i)
#        if len(PossibleGenes) == 0:
#            to_remove.append(pair)                
#    print(len(YoungExternal), len(to_remove), len(ExternalLike))
#
#    if len(to_remove) != 0:
#        for pair in to_remove:
#            YoungExternal.remove(pair)
#
#    # compute divergence between young nested pairs and their un-nested orthologs 
#    YoungInternalDiv = ComputeExpressionDivergenceOrthologs(YoungInternal, HumanExpression, ChimpExpression)
#    YoungExternalDiv = ComputeExpressionDivergenceOrthologs(YoungExternal, HumanExpression, ChimpExpression)
#    print(len(YoungInternalDiv), np.mean(YoungInternalDiv), len(YoungExternalDiv), np.mean(YoungExternalDiv), stats.ranksums(YoungInternalDiv, YoungExternalDiv)[1])
#    # compute divergence between external like and their un-nested orthologs
#    ExternalDiv = ComputeExpressionDivergenceOrthologs(ExternalLike, HumanExpression, ChimpExpression)
#    InternalDiv = ComputeExpressionDivergenceOrthologs(InternalLike, HumanExpression, ChimpExpression)
#    # compare expression divergence between young external and external-like and between young internal and internal-like
#    print(len(YoungExternalDiv), np.mean(YoungExternalDiv), len(ExternalDiv), np.mean(ExternalDiv), stats.ranksums(YoungExternalDiv, ExternalDiv)[1])
#    print(len(YoungInternalDiv), np.mean(YoungInternalDiv), len(InternalDiv), np.mean(InternalDiv), stats.ranksums(YoungInternalDiv, InternalDiv)[1])
#    ExternalReplicates.append(stats.ranksums(YoungExternalDiv, ExternalDiv)[1])
#    InternalReplicates.append(stats.ranksums(YoungInternalDiv, InternalDiv)[1])
#    
#
#print(min(ExternalReplicates), max(ExternalReplicates))
#total = 0
#for i in ExternalReplicates:
#    if i < 0.05:
#        total += 1
#print('% replicates with p 0.05 young external', (total / Nreplicates) * 100)
#print(min(InternalReplicates), max(InternalReplicates))
#total = 0
#for i in InternalReplicates:
#    if i < 0.05:
#        total += 1
#print('% replicates with p 0.05 young internal', (total / Nreplicates) * 100)
#
#
#
#
### use this function to sort young and ancestral nesting events
##def InferYoungOldNestingEvents(FirstSpOrthologs, SecondSpOverlapPairs, OutGroupOverlapPairs, FirstSpHostNestedPairs):
##    '''
#
#############################################
#
#
### compare expression divergence between old nested genes and their un-nested orthologs
### compare expression divergence between old host genes and their un-nested orthologs
##
### populate lists with old genes and their un-nested orthologs 
##OldInternal, OldExternal = [], []
### get pairs of human old external and internal genes and their nested chimp orthologs
##for pair in HumanOld:
##    # get the ortholog of the host and nested genes
##    extortho, internortho = HumanOrthologs[pair[0]][0], HumanOrthologs[pair[1]][0]
##    # check that these genes are nested
##    if extortho in ChimpNestedGeneSet and internortho in ChimpNestedGeneSet:
##        if extortho in ChimpExpression:
##            assert pair[0] in HumanExpression
##            OldExternal.append([pair[0], extortho])
##        if internortho in ChimpExpression:
##            assert pair[1] in HumanExpression
##            OldInternal.append([pair[1], internortho])
### get pairs of chimp young external and internal genes and their un-nested human orthologs      
##for pair in ChimpOld:
##    # get the ortholog of the host and nested genes
##    extortho, internortho = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
##    # check that these genes are nested
##    if extortho in HumanNestedGeneSet and internortho in HumanNestedGeneSet:
##        if extortho in HumanExpression:
##            assert pair[0] in ChimpExpression
##            OldExternal.append([extortho, pair[0]])
##        if internortho in HumanExpression:
##            assert pair[1] in ChimpExpression
##            OldInternal.append([internortho, pair[1]])
##print(len(HumanOld), len(ChimpOld), len(OldInternal), len(OldExternal), len(HumanOld) + len(ChimpOld), len(OldInternal) + len(OldExternal))
##
##
### compute divergence between old nested pairs and their nested orthologs 
##OldInternalDiv = ComputeExpressionDivergenceOrthologs(OldInternal, HumanExpression, ChimpExpression)
##OldExternalDiv = ComputeExpressionDivergenceOrthologs(OldExternal, HumanExpression, ChimpExpression)
##print(len(OldInternalDiv), np.mean(OldInternalDiv), len(OldExternalDiv), np.mean(OldExternalDiv), stats.ranksums(OldInternalDiv, OldExternalDiv)[1])
##
### create pairs of random internal-like and external-like genes in human and their un-nested orthologs in chimp
##InternalLike, ExternalLike = [], []
##
### create a list of pairs to remove when genes have no match
##to_remove = []
##
### for each human gene, match a random un-nested gene with similar characterisitics
##for pair in OldInternal:
##    # get the chromosome of the human gene
##    chromo = HumanGeneCoord[pair[0]][0]
##    # create a list of genes corresponding to all possible genes on chromo to draw from
##    PossibleGenes = list(ToDrawFrom[chromo].keys())
##    # draw a random gene on that chromo with matching characteristics
##    NotFound = True
##    while len(PossibleGenes) != 0 and NotFound == True:
##        i = random.choice(PossibleGenes)
##        gene = ToDrawFrom[chromo][i]
##        # assert gene not nested
##        assert gene not in HumanNestedGeneSet
##        # match gene by expression specificity (+- 0.01)
##        if gene in HumanExpSpecificity:
##            assert gene in HumanExpression
##            # match by expression specificity
##            if HumanExpSpecificity[pair[0]] - 0.05 <= HumanExpSpecificity[gene] <= HumanExpSpecificity[pair[0]] + 0.05:
##                # match by tissue with highest expression
##                if HumanExpression[gene].index(max(HumanExpression[gene])) == HumanExpression[pair[0]].index(max(HumanExpression[pair[0]])):
##                    # check that matching gene has a un-nested chimp ortholog
##                    if gene in HumanOrthologs:
##                        # check that ortholog is not nested and is expressed
##                        if HumanOrthologs[gene][0] not in ChimpNestedGeneSet and HumanOrthologs[gene][0] in ChimpExpression:
##                            # found internal-like gene, populate list
##                            InternalLike.append([gene, HumanOrthologs[gene][0]])
##                            # update boolean
##                            NotFound = False
##                        else:
##                            PossibleGenes.remove(i)
##                    else:
##                        PossibleGenes.remove(i)
##                else:
##                    PossibleGenes.remove(i)
##            else:
##                PossibleGenes.remove(i)
##        else:
##            PossibleGenes.remove(i)
##    if len(PossibleGenes) == 0:
##        to_remove.append(pair)                
##print(len(OldInternal), len(to_remove), len(InternalLike))
##
##if len(to_remove) != 0:
##    for pair in to_remove:
##        OldInternal.remove(pair)
##
##
### create a list of pairs to remove when genes have no match
##to_remove = []
##
### for each human gene, match a random un-nested gene with similar characterisitics
##for pair in OldExternal:
##    # get the chromosome of the human gene
##    chromo = HumanGeneCoord[pair[0]][0]
##    # create a list of genes corresponding to all possible genes on chromo to draw from
##    PossibleGenes = list(ToDrawFrom[chromo].keys())
##    # draw a random gene on that chromo with matching characteristics
##    NotFound = True
##    while len(PossibleGenes) != 0 and NotFound == True:
##        i = random.choice(PossibleGenes)
##        gene = ToDrawFrom[chromo][i]
##        # assert gene not nested
##        assert gene not in HumanNestedGeneSet
##        # match gene by expression specificity (+- 0.01)
##        if gene in HumanExpSpecificity:
##            assert gene in HumanExpression
##            # match by expression specificity
##            if HumanExpSpecificity[pair[0]] - 0.05 <= HumanExpSpecificity[gene] <= HumanExpSpecificity[pair[0]] + 0.05:
##                # match by tissue with highest expression
##                if HumanExpression[gene].index(max(HumanExpression[gene])) == HumanExpression[pair[0]].index(max(HumanExpression[pair[0]])):
##                    # check that matching gene has a un-nested chimp ortholog
##                    if gene in HumanOrthologs:
##                        # check that ortholog is not nested and is expressed
##                        if HumanOrthologs[gene][0] not in ChimpNestedGeneSet and HumanOrthologs[gene][0] in ChimpExpression:
##                            # found internal-like gene, populate list
##                            ExternalLike.append([gene, HumanOrthologs[gene][0]])
##                            # update boolean
##                            NotFound = False
##                        else:
##                            PossibleGenes.remove(i)
##                    else:
##                        PossibleGenes.remove(i)
##                else:
##                    PossibleGenes.remove(i)
##            else:
##                PossibleGenes.remove(i)
##        else:
##            PossibleGenes.remove(i)
##    if len(PossibleGenes) == 0:
##        to_remove.append(pair)                
##print(len(OldExternal), len(to_remove), len(ExternalLike))
##
##if len(to_remove) != 0:
##    for pair in to_remove:
##        OldExternal.remove(pair)
##
### compute divergence between young nested pairs and their un-nested orthologs 
##OldInternalDiv = ComputeExpressionDivergenceOrthologs(OldInternal, HumanExpression, ChimpExpression)
##OldExternalDiv = ComputeExpressionDivergenceOrthologs(OldExternal, HumanExpression, ChimpExpression)
##print(len(OldInternalDiv), np.mean(OldInternalDiv), len(OldExternalDiv), np.mean(OldExternalDiv), stats.ranksums(OldInternalDiv, OldExternalDiv)[1])
##
### compute divergence between external like and their un-nested orthologs
##ExternalDiv = ComputeExpressionDivergenceOrthologs(ExternalLike, HumanExpression, ChimpExpression)
##InternalDiv = ComputeExpressionDivergenceOrthologs(InternalLike, HumanExpression, ChimpExpression)
### compare expression divergence between young external and external-like and between young internal and internal-like
##print(len(OldExternalDiv), np.mean(OldExternalDiv), len(ExternalDiv), np.mean(ExternalDiv), stats.ranksums(OldExternalDiv, ExternalDiv)[1])
##print(len(OldInternalDiv), np.mean(OldInternalDiv), len(InternalDiv), np.mean(InternalDiv), stats.ranksums(OldInternalDiv, InternalDiv)[1])
##
##
