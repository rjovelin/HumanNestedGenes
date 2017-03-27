# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""


# use this script to plot expression divergence between host-nested pairs and
# their un-nested orthologs and expression divergence between external genes and
# their un-nested orthologs and between internal and their un-nested orthologs

# usage python3 PlotExpDivergYoungOld.py [options]
# [mouse/chimp]: sister-species of interest
# -[pairs/orthos]: expression divergence within gene pairs or between orthologs

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


# get option from command
SisterSp = sys.argv[1]
Analysis = sys.argv[2]
assert SisterSp in ['mouse', 'chimp']
assert Analysis in ['pairs', 'orthos']

# load dictionaries of overlapping genes
if SisterSp == 'chimp':
    # sister species is chimp and outgroup is gorilla
    jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
                 'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json',
                 'GorillaOverlappingGenes.json', 'GorillaNestedGenes.json']
    # get GFF file
    GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3', 'Gorilla_gorilla.gorGor3.1.86.gff3']
elif SisterSp == 'mouse':
    # sister species is mouse and outgroup is dog
    jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
                 'MouseOverlappingGenes.json', 'MouseNestedGenes.json',
                 'DogOverlappingGenes.json', 'DogNestedGenes.json']
    # get GFF file
    GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Mus_musculus.GRCm38.86.gff3', 'Canis_familiaris.CanFam3.1.87.gff3']
    

# make a list of dictionaries
AllOverlap = []
# loop over files
for i in range(len(jsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(jsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    AllOverlap.append(overlapping)

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

HumanInferredPairs = [HumanOld, HumanYoung]
SisterSpInferredPairs = [SisterSpOld, SisterSpYoung]


if Analysis == 'pairs':
   
    # compare expression divergence between human host and nested genes and their un-nested orthologs in sister-species   
    # remove human pairs if orthologs are nested in sister-species
    print(len(HumanYoung))    
    to_remove = []
    for pair in HumanYoung:
        if OrthoPairs[pair[0]] in NestedSets[1] or OrthoPairs[pair[1]] in NestedSets[1]:
            to_remove.append(pair)
    for pair in to_remove:
        HumanYoung.remove(pair)
    print(len(HumanYoung)) 

#    # make a list of mouse genes if genes not nested in mouse but nested in human and outgroup
#    SecondOverlap, OutGroupOverlap = [], []
#    for pair in SisterPairs[1]:
#        SecondOverlap.append(set(pair))
#    for pair in OutGroupPairs[1]:
#        OutGroupOverlap.append(set(pair))
#    to_remove = set()
#    for pair in HumanPairs[1]:
#        if pair[0] in OrthoPairs and pair[1] in OrthoPairs:
#            hostorthosp2, nestedorthosp2 = OrthoPairs[pair[0]], OrthoPairs[pair[1]]        
#            if {hostorthosp2, nestedorthosp2} not in SecondOverlap:
#                if pair[0] not in OrthoTrios or pair[1] not in OrthoTrios:
#                    to_remove.add(pair[0])
#                    to_remove.add(pair[1])
#                elif pair[0] in OrthoTrios and pair[1] in OrthoTrios:
#                    hostorthoout, nestedorthoout = OrthoTrios[pair[0]][1], OrthoTrios[pair[1]][1]                 
#                    if {hostorthoout, nestedorthoout} in OutGroupOverlap:
#                        to_remove.add(pair[0])
#                        to_remove.add(pair[1])
#    to_delete = []    
#    for pair in HumanYoung:
#        if pair[0] in to_remove or pair[1] in to_remove:
#            to_delete.append(pair)
#    for pair in to_delete:
#        HumanYoung.remove(pair)            
#    print(len(HumanYoung))
    
    
    
    
    # remove human pairs if genes are not expressed
    to_remove = [pair for pair in HumanYoung if pair[0] not in HumanExpression or pair[1] not in HumanExpression]
    for pair in to_remove:
        HumanYoung.remove(pair)
    print(len(HumanYoung))
    # get the sister-species un-nested orthologs
    SisterSpUnested = []
    for pair in HumanYoung:
        # check if orthologs are expressed in sister species
        if OrthoPairs[pair[0]] in SisterSpExpression and OrthoPairs[pair[1]] in SisterSpExpression:
            SisterSpUnested.append([OrthoPairs[pair[0]], OrthoPairs[pair[1]]])
    # generate a dict to draw random genes in sister-species
    SisterRandomGenes = GenerateAllUnNestedGenes(NestedSets[1], AllOrdered[1], SisterSpExpression)
    # make a list of control un-nested pairs
    ControlPairs = []
    for pair in SisterSpUnested:
        # make a list of matching gene pairs (orientation, chromosome, distance)
        PairPool = GenerateMatchingPoolPairs(pair, SisterRandomGenes, AllCoordinates[1], 2000)
        # draw a matching gene pair at random
        i = random.randint(0, len(PairPool) -1)
        ControlPairs.append(PairPool[i])
                        
    # compute expression divergence betwen human nested gene pairs
    HumanDiv = ComputeExpressionDivergenceGenePairs(HumanYoung, HumanExpression)
    # compute expression divergence between gene pairs in sister species    
    SisterSpDiv = ComputeExpressionDivergenceGenePairs(SisterSpUnested, SisterSpExpression)
    # compute expression divergence between control pairs    
    ControlDiv = ComputeExpressionDivergenceGenePairs(ControlPairs, SisterSpExpression)    
    
    HumanDiv.sort()
    SisterSpDiv.sort()
    ControlDiv.sort()
         
    
  
    P = PermutationResampling(HumanDiv, SisterSpDiv, 10000, statistic = np.mean)
    print(len(HumanDiv), len(SisterSpDiv), np.mean(HumanDiv), np.mean(SisterSpDiv), P)
    P = PermutationResampling(HumanDiv, ControlDiv, 10000, statistic = np.mean)
    print(len(HumanDiv), len(ControlDiv), np.mean(HumanDiv), np.mean(ControlDiv), P)
    P = PermutationResampling(SisterSpDiv, ControlDiv, 10000, statistic = np.mean)
    print(len(SisterSpDiv), len(ControlDiv), np.mean(SisterSpDiv), np.mean(ControlDiv), P)
    
    


elif Analysis == 'orthos':
    # compare expression divergence between young nested genes and their un-nested orthologs
    # compare expression divergence between young host genes and their un-nested orthologs


## remove gene pairs with genes lacking expression
#HumanOld = FilterGenePairsWithoutExpression(HumanOld, HumanExpression)
#HumanYoung = FilterGenePairsWithoutExpression(HumanYoung, HumanExpression)
#ChimpOld = FilterGenePairsWithoutExpression(ChimpOld, ChimpExpression)
#ChimpYoung = FilterGenePairsWithoutExpression(ChimpYoung, ChimpExpression)



    # make sets of external and internal genes [hsaextold, hsaintold, hsaextyoung, hsaintyoung]
    HumanExtIntGenes = []
    for i in range(len(HumanInferredPairs)):
        external, internal = set(), set()
        for pair in HumanInferredPairs[i]:
            external.add(pair[0])
            internal.add(pair[1])
        HumanExtIntGenes.append(external)
        HumanExtIntGenes.append(internal)

    # make sets of external and internal genes [ptrextold, ptrintold, ptrextyoung, ptrintyoung]
    ChimpExtIntGenes = []
    for i in range(len(ChimpInferredPairs)):
        external, internal = set(), set()
        for pair in ChimpInferredPairs[i]:
            external.add(pair[0])
            internal.add(pair[1])
        ChimpExtIntGenes.append(external)
        ChimpExtIntGenes.append(internal)

    # for young external and internal, remove genes if ortholog is nested
    for i in range(2, len(HumanExtIntGenes)):
        to_remove = set()
        for gene in HumanExtIntGenes[i]:
            assert gene in OrthoPairs
            if OrthoPairs[gene] in NestedSets[1]:
                to_remove.add(gene)
            if gene in OrthoTrios:
                if OrthoTrios[gene][0] in NestedSets[1] or OrthoTrios[gene][1] in NestedSets[2]:
                    to_remove.add(gene)
        for gene in to_remove:
            HumanExtIntGenes[i].remove(gene)             

    # for young external and internal, remove genes if ortholog is nested
    for i in range(2, len(ChimpExtIntGenes)):
        to_remove = set()
        for gene in ChimpExtIntGenes[i]:
            assert gene in ChimpOrthos
            if ChimpOrthos[gene] in NestedSets[0]:
                to_remove.add(gene)
            if gene in ChimpOrthoTrios:
                if ChimpOrthoTrios[gene][0] in NestedSets[0] or ChimpOrthoTrios[gene][1] in NestedSets[2]:
                    to_remove.add(gene)
        for gene in to_remove:
            ChimpExtIntGenes[i].remove(gene)             

    to_remove = set()
    # remove non-overlapping genes if their ortholog are overlapping
    for gene in NonOverlappingSets[0]:
        # check if ortholog is overlapping in chimp or gorilla
        if gene in OrthoPairs and OrthoPairs[gene] in OverlapSets[1]:
            to_remove.add(gene)
        if gene in OrthoTrios:
            if OrthoTrios[gene][0] in OverlapSets[1] or OrthoTrios[gene][1] in OverlapSets[2]:
                to_remove.add(gene)
    for gene in to_remove:
        NonOverlappingSets[0].remove(gene)

    # remove non-overlapping genes genes without orthologs in chimp
    to_remove = [gene for gene in NonOverlappingSets[0] if gene not in OrthoPairs]
    for gene in to_remove:
        NonOverlappingSets[0].remove(gene)

    # remove non-overlapping genes without expression
    to_remove = [gene for gene in NonOverlappingSets[0] if gene not in HumanExpression]
    for gene in to_remove:
        NonOverlappingSets[0].remove(gene)

    #HsaGenes = [NonOverlappingSets[0]]
    #HsaGenes.extend(ExtIntGenes)

    # remove genes without expression
    for i in range(len(HumanExtIntGenes)):
        to_remove = [gene for gene in HumanExtIntGenes[i] if gene not in HumanExpression]
        for gene in to_remove:
            HumanExtIntGenes[i].remove(gene)
        HumanExtIntGenes[i] = set(HumanExtIntGenes[i])
    for i in range(len(ChimpExtIntGenes)):
        to_remove = [gene for gene in ChimpExtIntGenes[i] if gene not in ChimpExpression]
        for gene in to_remove:
            ChimpExtIntGenes[i].remove(gene)
        ChimpExtIntGenes[i] = set(ChimpExtIntGenes[i])
    
    # remove gene "duplicates" by removing chimp genes with ortologs already present in each group
    for i in range(len(ChimpExtIntGenes)):
        to_remove = []
        for gene in ChimpExtIntGenes[i]:
            if ChimpOrthos[gene] in HumanExtIntGenes[i]:
                to_remove.append(gene)
        for gene in to_remove:
            ChimpExtIntGenes[i].remove(gene)

    # add set of non-overlapping genes to list of external and internal genes
    HumanExtIntGenes.insert(0, NonOverlappingSets[0])
    # make lists with human and chimp orthologs
    HumanChimpPairs = []
    for i in range(len(HumanExtIntGenes)):
        pairs = [[gene, OrthoPairs[gene]] for gene in HumanExtIntGenes[i] if gene in HumanExpression and OrthoPairs[gene] in ChimpExpression]
        HumanChimpPairs.append(pairs)
    ChimpHumanPairs = []
    for i in range(len(ChimpExtIntGenes)):
        pairs = [[gene, ChimpOrthos[gene]] for gene in ChimpExtIntGenes[i] if gene in ChimpExpression and ChimpOrthos[gene] in HumanExpression]
        ChimpHumanPairs.append(pairs)

    ExpDivergence = []
    for i in range(len(HumanChimpPairs)):
        D = ComputeExpressionDivergenceOrthologs(HumanChimpPairs[i], HumanExpression, ChimpExpression)
        ExpDivergence.append(D)
    
    ## add divergence for chimp-specific nesting events
    #ExpDivergence[-2].extend(ComputeExpressionDivergenceOrthologs(ChimpHumanPairs[-2], ChimpExpression, HumanExpression))
    #ExpDivergence[-1].extend(ComputeExpressionDivergenceOrthologs(ChimpHumanPairs[-1], ChimpExpression, HumanExpression))
    
    
    for i in range(1, len(ExpDivergence)):
        P = PermutationResampling(ExpDivergence[0], ExpDivergence[i], 10000, statistic = np.mean)
        print(i, len(ExpDivergence[i]), np.mean(ExpDivergence[0]), np.mean(ExpDivergence[i]), P)

    print('merge young and old external and internal genes')

    # merge external and internal for the same age group
    Old, Young = [i for i in ExpDivergence[1]], [i for i in ExpDivergence[3]]
    Old.extend(ExpDivergence[2])
    Young.extend(ExpDivergence[4])

    ExpDiv = [ExpDivergence[0], Old, Young]
    for i in range(1, len(ExpDiv)):
        P = PermutationResampling(ExpDiv[0], ExpDiv[i], 10000, statistic = np.mean)
        print(i, len(ExpDiv[i]), np.mean(ExpDiv[0]), np.mean(ExpDiv[i]), P)
    
