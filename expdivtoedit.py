# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 15:23:55 2016

@author: RJovelin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""


# use this script to plot expression divergence before young and old overlapped/nested genes

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
        

## get GFF file
#HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
#PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
#GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
#
## make a list of primate GFF files
#GFFs = [HsaGFF, PtrGFF, GgoGFF]
## make a list of species names
#SpeciesNames = ['human', 'chimp', 'gorilla']
## make a list of host:nested genes dictionaries
#HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes]
## make a list of host:contained genes dictionaries
#HostContained = [HumanContainedGenes, ChimpContainedGenes, GorillaContainedGenes]
## make a list of overlapping genes dictionaries
#Overlapping = [HumanOverlappingGenes, ChimpOverlappingGenes, GorillaOverlappingGenes]
#
#
## loop over GFF files
#for i in range(len(GFFs)):
#    print(GFFs[i][:GFFs[i].index('.')], SpeciesNames[i])
#    # get expression profile of the species genes
#    SpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[i])
#    # remove genes wuthout expression
#    SpExpression = RemoveGenesLackingExpression(SpExpression)
#    # get relative expression
#    SpExpression = TransformRelativeExpression(SpExpression)
#    # make a list of host-nested gene pairs
#    SpHostNestedPairs = GetHostNestedPairs(HostGenes[i])
#    # make a list of host-contained gene pairs
#    SpHostContainedPairs = GetHostNestedPairs(HostContained[i])
#    # make a list of overlapping gene pairs
#    SpOverlappingPairs = GetHostNestedPairs(Overlapping[i])
#    # assert each gene pair is found only once in each level of genomic organization
#    for pair in SpHostNestedPairs:
#        assert SpHostNestedPairs.count(pair) == 1
#    for pair in SpHostContainedPairs:
#        assert SpHostContainedPairs.count(pair) == 1
#    for pair in SpOverlappingPairs:
#        assert SpOverlappingPairs.count(pair) == 1
#    
#    print('total number of host-nested pairs', SpeciesNames[i], len(SpHostNestedPairs))
#    print('total number of host-contained pairs', SpeciesNames[i], len(SpHostContainedPairs))
#    print('total number of overlapping pairs', SpeciesNames[i], len(SpOverlappingPairs))
#    # remove gene pairs from higher hierarchical level present in lower hierarchical level    
#    # count gene pairs before filtering gene pairs
#    a, b = len(SpOverlappingPairs), len(SpHostContainedPairs)
#    # make a list of overlaping genes that are not contained
#    SpOverlappingPairs = RemoveGenePairsFromHigherLevel(SpOverlappingPairs, SpHostContainedPairs)
#    assert a == len(SpHostContainedPairs) + len(SpOverlappingPairs)    
#    # make a list of contained genes that are not host:intronic nested genes
#    SpHostContainedPairs = RemoveGenePairsFromHigherLevel(SpHostContainedPairs, SpHostNestedPairs)
#    assert b == len(SpHostContainedPairs) + len(SpHostNestedPairs)  
#    print('number of overlapping-not contained gene pairs', SpeciesNames[i], len(SpOverlappingPairs))
#    print('number of contained-not nested gene pairs', SpeciesNames[i], len(SpHostContainedPairs))
#    # remove gene pairs with genes lacking expression
#    SpHostNestedPairs = FilterGenePairsWithoutExpression(SpHostNestedPairs, SpExpression)
#    SpHostContainedPairs = FilterGenePairsWithoutExpression(SpHostContainedPairs, SpExpression)
#    SpOverlappingPairs = FilterGenePairsWithoutExpression(SpOverlappingPairs, SpExpression)
#    print('number of host-nested pairs with expression', SpeciesNames[i], len(SpHostNestedPairs))
#    print('number of host-contained pairs with expression', SpeciesNames[i], len(SpHostContainedPairs))
#    print('number of overlapping pairs with expression', SpeciesNames[i], len(SpOverlappingPairs))



# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
 

# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF]
# make a list of species names
SpeciesNames = ['human', 'chimp', 'gorilla']



# make a list of host:nested genes dictionaries
HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes]


#HostGenes = [HumanContainedGenes, ChimpContainedGenes, GorillaContainedGenes]

#HostGenes = [HumanOverlappingGenes, ChimpOverlappingGenes, GorillaOverlappingGenes]




# get expression profile of the species genes
HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[0])
# remove genes wuthout expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
# make a set of host and nested genes (include all host and nested even if not expressed)    
HumanNestedConformation = MakeFullPartialOverlapGeneSet(HostGenes[0])    
# make a list of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HostGenes[0])
# remove gene pairs with genes lacking expression
HumanHostNestedPairs = FilterGenePairsWithoutExpression(HumanHostNestedPairs, HumanExpression)


# get expression profile of the species genes
ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[1])
# remove genes wuthout expression
ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
# get relative expression
ChimpExpression = TransformRelativeExpression(ChimpExpression)
# make a set of host and nested genes (include all host and nested even if not expressed)    
ChimpNestedConformation = MakeFullPartialOverlapGeneSet(HostGenes[1])    
# make a list of host-nested gene pairs
ChimpHostNestedPairs = GetHostNestedPairs(HostGenes[1])
# remove gene pairs with genes lacking expression
ChimpHostNestedPairs = FilterGenePairsWithoutExpression(ChimpHostNestedPairs, ChimpExpression)


# get expression profile of the species genes
GorillaExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[2])
# remove genes wuthout expression
GorillaExpression = RemoveGenesLackingExpression(GorillaExpression)
# get relative expression
GorillaExpression = TransformRelativeExpression(GorillaExpression)
# make a set of host and nested genes (include all host and nested even if not expressed)    
GorillaNestedConformation = MakeFullPartialOverlapGeneSet(HostGenes[2])    
# make a list of host-nested gene pairs
GorillaHostNestedPairs = GetHostNestedPairs(HostGenes[2])
# remove gene pairs with genes lacking expression
GorillaHostNestedPairs = FilterGenePairsWithoutExpression(GorillaHostNestedPairs, GorillaExpression)


HumanOrthologs = ParseOrthologFile('HumanChimpGorillaOrthologs.txt')
# create a dicti with chimp genes as key {chimp_gene: [human_orto, gorilla_ortho]}
ChimpOrthologs = {}
for gene in HumanOrthologs:
    ChimpOrthologs[HumanOrthologs[gene][0]] = [gene, HumanOrthologs[gene][1]]


# make lists of overlapping gene pairs
HumanOverlappingPairs = GetHostNestedPairs(HumanOverlappingGenes)
ChimpOverlappingPairs = GetHostNestedPairs(ChimpOverlappingGenes)
GorillaOverlappingPairs = GetHostNestedPairs(GorillaOverlappingGenes)

HumanContainedPairs = GetHostNestedPairs(HumanContainedGenes)
ChimpContainedPairs = GetHostNestedPairs(ChimpContainedGenes)
GorillaContainedPairs = GetHostNestedPairs(GorillaContainedGenes)


# make lists of old and yound host:nested pairs
#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HostGenes[0], HostGenes[1], HostGenes[2], HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HostGenes[1], HostGenes[0], HostGenes[2], ChimpHostNestedPairs)
 


#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HumanContainedPairs, ChimpContainedPairs, GorillaContainedPairs, HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, ChimpContainedPairs, HumanContainedPairs, GorillaContainedPairs, ChimpHostNestedPairs)

HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, ChimpOverlappingPairs, GorillaOverlappingPairs, HumanHostNestedPairs)
ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HumanOverlappingPairs, GorillaOverlappingPairs, ChimpHostNestedPairs)

print('human nested pairs', len(HumanHostNestedPairs))
print('human old pairs', len(HumanOld))
print('human young', len(HumanYoung))

# make list of ancestral un-nested gene pairs
humanancestral, chimpancestral = [], []

to_remove = []
for pair in HumanYoung:
    ortho1, ortho2 = HumanOrthologs[pair[0]][0], HumanOrthologs[pair[1]][0]
    if ortho1 not in ChimpExpression or ortho2 not in ChimpExpression:
        to_remove.append(pair)
for pair in to_remove:
    HumanYoung.remove(pair)
        
for pair in HumanYoung:
    ortho1, ortho2 = HumanOrthologs[pair[0]][0], HumanOrthologs[pair[1]][0]
    assert ortho1 in ChimpExpression
    assert ortho2 in ChimpExpression
    chimpancestral.append([ortho1, ortho2])


to_remove = []
for pair in ChimpYoung:
    ortho1, ortho2 = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
    if ortho1 not in HumanExpression or ortho2 not in HumanExpression:
        to_remove.append(pair)
for pair in to_remove:
    ChimpYoung.remove(pair)
    
for pair in ChimpYoung:
    ortho1, ortho2 = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
    assert ortho1 in HumanExpression
    assert ortho2 in HumanExpression
    humanancestral.append([ortho1, ortho2])


print(len(HumanOld), len(HumanYoung), len(humanancestral))
print(len(ChimpOld), len(ChimpYoung), len(chimpancestral))

# compute divergence between young human nested pairs and their un-nested orthologs in chimp
HumanYoungDiv = ComputeExpressionDivergenceGenePairs(HumanYoung, HumanExpression)
ChimpUnNestedDiv = ComputeExpressionDivergenceGenePairs(chimpancestral, ChimpExpression)

# compute divergence between young chimp nested pairs and their un-nested orthologs in human
ChimpYoungDiv = ComputeExpressionDivergenceGenePairs(ChimpYoung, ChimpExpression)
HumanUnNestedDiv = ComputeExpressionDivergenceGenePairs(humanancestral, HumanExpression)

print(len(HumanYoungDiv), np.mean(HumanYoungDiv), len(ChimpUnNestedDiv), np.mean(ChimpUnNestedDiv), stats.ranksums(HumanYoungDiv, ChimpUnNestedDiv)[1])
print(len(ChimpYoungDiv), np.mean(ChimpYoungDiv), len(HumanUnNestedDiv), np.mean(HumanUnNestedDiv), stats.ranksums(ChimpYoungDiv, HumanUnNestedDiv)[1])







########### compare expression divergence between young nested genes and their un-nested orthologs
########### compare expression divergence between youn host genes and their un-nested orthologs

YoungInternal, YoungExternal = [], []

# make lists of old and yound host:nested pairs
#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HostGenes[0], HostGenes[1], HostGenes[2], HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HostGenes[1], HostGenes[0], HostGenes[2], ChimpHostNestedPairs)



#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HumanOverlappingPairs, ChimpOverlappingPairs, GorillaOverlappingPairs, HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, ChimpOverlappingPairs, HumanOverlappingPairs, GorillaOverlappingPairs, ChimpHostNestedPairs)

#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HumanContainedPairs, ChimpContainedPairs, GorillaContainedPairs, HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, ChimpContainedPairs, HumanContainedPairs, GorillaContainedPairs, ChimpHostNestedPairs)

HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, ChimpOverlappingPairs, GorillaOverlappingPairs, HumanHostNestedPairs)
ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HumanOverlappingPairs, GorillaOverlappingPairs, ChimpHostNestedPairs)


#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, ChimpHostNestedPairs, GorillaHostNestedPairs, HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HumanHostNestedPairs, GorillaHostNestedPairs, ChimpHostNestedPairs)




for pair in HumanYoung:
    # get the ortholog of the host and nested genes
    extortho, internortho = HumanOrthologs[pair[0]][0], HumanOrthologs[pair[1]][0]
    if extortho in ChimpExpression:
        assert pair[0] in HumanExpression
        YoungExternal.append([pair[0], extortho])
    if internortho in ChimpExpression:
        assert pair[1] in HumanExpression
        YoungInternal.append([pair[1], internortho])
        
for pair in ChimpYoung:
    # get the ortholog of the host and nested genes
    extortho, internortho = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
    if extortho in HumanExpression:
        assert pair[0] in ChimpExpression
        YoungExternal.append([extortho, pair[0]])
    if internortho in HumanExpression:
        assert pair[1] in ChimpExpression
        YoungInternal.append([internortho, pair[1]])
        
print(len(YoungInternal), len(YoungExternal))


# compute divergence between young human nested pairs and their un-nested orthologs in chimp
YoungInternalDiv = ComputeExpressionDivergenceOrthologs(YoungInternal, HumanExpression, ChimpExpression)
YoungExternalDiv = ComputeExpressionDivergenceOrthologs(YoungExternal, HumanExpression, ChimpExpression)

print(len(YoungInternalDiv), np.mean(YoungInternalDiv), len(YoungExternalDiv), np.mean(YoungExternalDiv), stats.ranksums(YoungInternalDiv, YoungExternalDiv)[1])

######################################################
######################################################
###### compare sequence divergence between young nested genes and their un-nested orthologs
###### compare sequence divergence between young host genes and their un-nested orthologs



## make a list of host-nested gene pairs
#HumanHostNestedPairs = GetHostNestedPairs(HostGenes[0])
#ChimpHostNestedPairs = GetHostNestedPairs(HostGenes[1])
#GorillaHostNestedPairs = GetHostNestedPairs(HostGenes[2])





# parse the se divergence file to extract dN
HumandN, ChimpdN = {}, {}
infile = open('HumanChimpSeqDiverg.txt')
infile.readline()
for line in infile:
    if line.startswith('ENS'):
        line = line.rstrip().split('\t')
        HumandN[line[0]] = float(line[2])
        ChimpdN[line[1]] = float(line[2])
infile.close()


## parse the se divergence file to extract dN
#HumandN, ChimpdN = {}, {}
#infile = open('HumanChimpSeqDiverg.txt')
#infile.readline()
#for line in infile:
#    if line.startswith('ENS'):
#        line = line.rstrip().split('\t')
#        if line[-1] != 'NA':
#            HumandN[line[0]] = float(line[-1])
#            ChimpdN[line[1]] = float(line[-1])
#infile.close()








YoungInternalSeq, YoungExternalSeq = [], []

# make lists of old and yound host:nested pairs
#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HostGenes[0], HostGenes[1], HostGenes[2], HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HostGenes[1], HostGenes[0], HostGenes[2], ChimpHostNestedPairs)



#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HumanOverlappingPairs, ChimpOverlappingPairs, GorillaOverlappingPairs, HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, ChimpOverlappingPairs, HumanOverlappingPairs, GorillaOverlappingPairs, ChimpHostNestedPairs)

#HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, HumanContainedPairs, ChimpContainedPairs, GorillaContainedPairs, HumanHostNestedPairs)
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, ChimpContainedPairs, HumanContainedPairs, GorillaContainedPairs, ChimpHostNestedPairs)

HumanOld, HumanYoung = InferYoungOldNestingEvents(HumanOrthologs, ChimpOverlappingPairs, GorillaOverlappingPairs, HumanHostNestedPairs)
ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HumanOverlappingPairs, GorillaOverlappingPairs, ChimpHostNestedPairs)


 
for pair in HumanYoung:
    # get dN for the host and nested gene
    if pair[0] in HumandN:
        YoungExternalSeq.append(HumandN[pair[0]])
    if pair[1] in HumandN:
        YoungInternalSeq.append(HumandN[pair[1]])
        
        
for pair in ChimpYoung:
    # get dN for the host and nested genes
    if pair[0] in ChimpdN:
        YoungExternalSeq.append(ChimpdN[pair[0]])
    if pair[1] in ChimpdN:
        YoungInternalSeq.append(ChimpdN[pair[1]])


print(len(YoungInternalSeq), np.mean(YoungInternalSeq), len(YoungExternalSeq), np.mean(YoungExternalSeq), stats.ranksums(YoungInternalSeq, YoungExternalSeq)[1])






















