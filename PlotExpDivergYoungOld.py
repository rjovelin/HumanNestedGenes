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
    print(len(pairs))
    AllPairs.append(pairs)
# make pairs of overlapping genes
HumanPairs = AllPairs[:2]
ChimpPairs = AllPairs[2:4]
GorillaPairs = AllPairs[4:]


# make lists of host and nested genes in each species
inthuman, exthuman, intchimp, extchimp, intgorilla, extgorilla = [], [] , [] ,[], [], []
for pair in HumanPairs[1]:
    exthuman.append(pair[0])
    inthuman.append(pair[1])
for pair in ChimpPairs[1]:
    extchimp.append(pair[0])
    intchimp.append(pair[1])
for pair in GorillaPairs[1]:
    extgorilla.append(pair[0])
    intgorilla.append(pair[1])


# make list with sets of non-overlapping genes
NonOverlappingSets = []
for i in range(3):
    j = i * 2
    print(i, j, jsonFiles[j])
    # make a set of non-overlapping gene
    nonoverlap = MakeNonOverlappingGeneSet(AllOverlap[j], AllCoordinates[i])
    NonOverlappingSets.append(nonoverlap)    

# make sets of host and nested nested genes
NestedSets = []
for i in range(1, len(AllOverlap), 2):
    nestedset = MakeFullPartialOverlapGeneSet(AllOverlap[i])
    NestedSets.append(nestedset)


# get 1:1 orthologs between human anc chimp
OrthoPairs = MatchOrthologPairs('HumanChimpOrthologs.txt')
# get 1:1 orthologs between human, chimp and gorilla {human:[chimp,gorilla]}
OrthoTrios = MatchOrthologTrios('HumanChimpGorillaOrthologs.txt')

## remove gene pairs lacking orthos
#for i in range(len(HumanPairs)):
#    to_remove = []
#    for pair in HumanPairs[i]:
#        if pair[0] not in Orthos or pair[1] not in Orthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        HumanPairs[i].remove(pair)
## remove chimp pairs lacking orthologs
#for i in range(len(ChimpPairs)):
#    to_remove = []
#    for pair in ChimpPairs[i]:
#        if pair[0] not in chimporthos or pair[1] not in chimporthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        ChimpPairs[i].remove(pair)
## remove gorilla pairs lacking orthologs
#for i in range(len(GorillaPairs)):
#    to_remove = []
#    for pair in GorillaPairs[i]:
#        if pair[0] not in gorillaorthos or pair[1] not in gorillaorthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        GorillaPairs[i].remove(pair)


# infer young and old nesting events in human and chimp
HumanOld, HumanYoung = InferYoungOldNestingEvents(OrthoPairs, OrthoTrios, ChimpPairs[1], GorillaPairs[1], HumanPairs[1])
print(len(HumanOld), len(HumanYoung))
    
# get expression profile in human and chimp
HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
# remove genes without expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
ChimpExpression = TransformRelativeExpression(ChimpExpression)

# remove gene pairs with genes lacking expression
HumanOld = FilterGenePairsWithoutExpression(HumanOld, HumanExpression)
HumanYoung = FilterGenePairsWithoutExpression(HumanYoung, HumanExpression)

InferredPairs = [HumanOld, HumanYoung]
for i in InferredPairs:
    print(len(i))

# compare expression divergence between young nested genes and their un-nested orthologs
# compare expression divergence between young host genes and their un-nested orthologs

# make sets of external and internal genes [hsaextold, hsaintold, hsaextyoung, hsaintyoung]
ExtIntGenes = []
for i in range(len(InferredPairs)):
    external, internal = set(), set()
    for pair in InferredPairs[i]:
        external.add(pair[0])
        internal.add(pair[1])
    ExtIntGenes.append(external)
    ExtIntGenes.append(internal)
for i in range(len(ExtIntGenes)):
    print(len(ExtIntGenes[i]))

# for young external and internal, remove genes if ortholog is nested
for i in range(1, len(ExtIntGenes), 2):
    to_remove = set()
    for gene in ExtIntGenes[i]:
        assert gene in OrthoPairs
        if OrthoPairs[gene] in NestedSets[1]:
            to_remove.add(gene)
        if gene in OrthoTrios:
            if OrthoTrios[gene][0] in NestedSets[1] or OrthoTrios[gene][1] in NestedSets[2]:
                to_remove.add(gene)
    for gene in to_remove:
        ExtIntGenes[i].remove(gene)             
for i in range(len(ExtIntGenes)):
    print(len(ExtIntGenes[i]))

to_remove = set()
# remove non-overlapping genes if their ortholog are overlapping
for gene in NonOverlappingSets[0]:
    # check if ortholog is overlapping in chimp or gorilla
    if gene in OrthoPairs and OrthoPairs[gene] in NonOverlappingSets[1]:
        to_remove.add(gene)
    if gene in OrthoTrios:
        if OrthoTrios[gene][0] in NonOverlappingSets[1] or OrthoTrios[gene][1] in NonOverlappingSets[2]:
            to_remove.add(gene)
for gene in to_remove:
    NonOverlappingSets[0].remove(gene)

# remove non-overlapping genes genes without orthologs in chimp
to_remove = [gene for gene in NonOverlappingSets[0] if gene not in OrthoPairs]
for gene in to_remove:
    NonOverlappingSets[0].remove(gene)

# make a list of human genes [non-overlapping, extold, intold, extyoung, intyoung]
HsaGenes = [NonOverlappingSets[0]]
HsaGenes.extend(ExtIntGenes)

# remove genes without expression
for i in range(len(HsaGenes)):
    to_remove = [gene for gene in HsaGenes[i] if gene not in HumanExpression]
    for gene in to_remove:
        HsaGenes.remove(gene)
for i in range(len(HsaGenes)):
    print(len(HsaGenes[i]))        
    
for i in range(len(HsaGenes)):
    HsaGenes[i] = list(HsaGenes[i])
    D = ComputeExpressionDivergenceOrthologs(HsaGenes[i], HumanExpression, ChimpExpression)
    print(np.mean(D))


##############################################################





## -*- coding: utf-8 -*-
#"""
#Created on Wed Mar  1 11:30:11 2017
#
#@author: RJovelin
#"""
#
#
## use this script to plot the distances between non-overlapping orthologs of human overlapping genes
#
## import modules
## use Agg backend on server without X server
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#from matplotlib import rc
#import matplotlib.gridspec as gridspec
#rc('mathtext', default='regular')
#import json
#import random
#import copy
#import sys
#import os
#import math
#import numpy as np
#from scipy import stats
#from HsaNestedGenes import *
#
#
## load dictionaries of overlapping genes
#jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
#             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json']
#
## make a list of dictionaries
#AllOverlap = []
## loop over files
#for i in range(len(jsonFiles)):
#    # load dictionary of overlapping gene pairs
#    json_data = open(jsonFiles[i])
#    overlapping = json.load(json_data)
#    json_data.close()
#    AllOverlap.append(overlapping)
#
## get GFF file
#GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3']
#
## make a list of gene coordinates       
#AllCoordinates, AllOrdered = [], []
## loop over GFF files
#for i in range(len(GFF)):
#    # get the coordinates of genes on each chromo
#    # {chromo: {gene:[chromosome, start, end, sense]}}
#    GeneChromoCoord = ChromoGenesCoord(GFF[i])
#    # map each gene to its mRNA transcripts
#    MapGeneTranscript = GeneToTranscripts(GFF[i])
#    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
#    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
#    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
#    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
#    AllCoordinates.append(GeneCoord)
#    AllOrdered.append(OrderedGenes)
#HumanOrdered, ChimpOrdered = AllOrdered[0], AllOrdered[1]
#HumanCoord, ChimpCoord = AllCoordinates[0], AllCoordinates[1]
#
## get 1:1 orthologs between human and chimp
#Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
#
## make pairs of overlapping genes
#AllPairs = []
#for i in range(len(AllOverlap)):
#    pairs = GetHostNestedPairs(AllOverlap[i])
#    AllPairs.append(pairs)
## get the gene pairs
#HumanPairs = AllPairs[:2]
#ChimpPairs = AllPairs[2:]
#
#
#
#extmissing, intmissing = set(), set()
#for pair in HumanPairs[1]:
#    if pair[0] not in Orthos:
#        extmissing.add(pair[0])
#    if pair[1] not in Orthos:
#        intmissing.add(pair[1])
#print(len(extmissing), len(intmissing))
#
#
#print(len(HumanPairs[0]), len(HumanPairs[1]))
## remove human genes lacking orthologs
#for i in range(len(HumanPairs)):
#    to_remove = []
#    for pair in HumanPairs[i]:
#        if pair[0] not in Orthos or pair[1] not in Orthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        HumanPairs[i].remove(pair)
## remove order for chimp gene pairs
#for i in range(len(ChimpPairs)):
#    for j in range(len(ChimpPairs[i])):
#        ChimpPairs[i][j] = set(ChimpPairs[i][j])
#print(len(HumanPairs[0]), len(HumanPairs[1]))
#
#
#SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG, NstCons = [], [], [], [], [], []
#Conserved, NotConserved = 0, 0
#
#
## loop over human nested gene pairs
#for pair in HumanPairs[1]:
#    # get the orthologs of the human genes
#    ortho1, ortho2 = Orthos[pair[0]], Orthos[pair[1]]
#    # check if human gene pairs is nested in chimp
#    if set([ortho1, ortho2]) not in ChimpPairs[1]:
#        NotConserved += 1
#        # get chromosomes of chimp orthologs
#        chromo1, chromo2 = ChimpCoord[ortho1][0], ChimpCoord[ortho2][0]        
#        if chromo1 == chromo2:
#            # get start positions of chimp orthologs
#            S1, S2 = ChimpCoord[ortho1][1], ChimpCoord[ortho2][1]
#            # get end positions of chimp orthologs
#            E1, E2 = ChimpCoord[ortho1][2], ChimpCoord[ortho2][2]            
#            # get indices of chimp orthologs in the ordred gene list
#            P1, P2 = ChimpOrdered[chromo1].index(ortho1), ChimpOrdered[chromo2].index(ortho2)
#            # check if orthologs are overlapping in chimp
#            if set([ortho1, ortho2]) in ChimpPairs[0]:
#                # chimp genes are overlapping
#                # check if they are adjcent
#                if S1 < S2:
#                    assert P1 < P2
#                    if P2 != P1 + 1:
#                        # overlapping non adjacent
#                        SepOvlp.append(pair)
#                    elif P2 == P1 + 1:
#                        # overlapping adjacent
#                        AdjOvlp.append(pair)
#                elif S1 > S2:
#                    assert P2 < P1
#                    if P1 != P2 + 1:
#                        # overlapping non adjacent
#                        SepOvlp.append(pair)
#                    elif P1 == P2 + 1:
#                        # overlapping adjacent
#                        AdjOvlp.append(pair)
#                elif S1 == S2:
#                    print('merde')
#            elif set([ortho1, ortho2]) not in ChimpPairs[0]:
#                # chimp genes are not overlapping
#                # check if they are adjacent
#                if S1 < S2:
#                    assert P1 < P2
#                    if P2 != P1 + 1:
#                        # non-overlapping non-adjacent
#                        SepNonOvlp.append(pair)
#                    elif P2 == P1 + 1:
#                        # non-overlapping adjacent
#                        AdjNonOvlp.append(pair)
#                elif S1 > S2:
#                    assert P2 < P1
#                    if P1 != P2 + 1:
#                        # non-overlapping non adjacent
#                        SepNonOvlp.append(pair)
#                    elif P1 == P2 + 1:
#                        # non-overlapping adjacent
#                        AdjNonOvlp.append(pair)
#                elif S1 == S2:
#                    print('shit')
#              
#    elif set([ortho1, ortho2]) in ChimpPairs[1]:
#        # orthologs of human nested genes are also nested
#        Conserved += 1            
#
## make a list of counts
#Pairs = [SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG, NstCons]
#
## make a list of labels
#Labels = ['Separated overlapping', 'Separated non-overlapping', 'Adjacent overlapping', 'Adjacent non-overlapping', 'Different chromosomes', 'Nested']
## remove counts and labels equal to 0
#to_remove = []
#for i in range(len(Pairs)):
#    if len(Pairs[i]) != 0:
#        print(Labels[i])
#        print(Pairs[i][0][0], Pairs[i][0][1], Orthos[Pairs[i][0][0]], Orthos[Pairs[i][0][1]])
#        print(HumanCoord[Pairs[i][0][0]], HumanCoord[Pairs[i][0][1]], ChimpCoord[Orthos[Pairs[i][0][0]]], ChimpCoord[Orthos[Pairs[i][0][1]]])
#        
