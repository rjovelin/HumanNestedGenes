# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""


# use this script to plot expression divergence between host-nested pairs and
# their un-nested orthologs and expression divergence between external genes and
# their un-nested orthologs and between internal and their un-nested orthologs

# usage python3 PlotExpDivergYoungOld.py [options]
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
Analysis = sys.argv[1]
assert Analysis in ['pairs', 'orthos']


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

# get 1:1 orthologs between human anc chimp
OrthoPairs = MatchOrthologPairs('HumanChimpOrthologs.txt')
# get 1:1 orthologs between human, chimp and gorilla {human:[chimp,gorilla]}
OrthoTrios = MatchOrthologTrios('HumanChimpGorillaOrthologs.txt')

# reverse dict with human and chimp orthologs
ChimpOrthos = {}
for gene in OrthoPairs:
    ChimpOrthos[OrthoPairs[gene]] = gene
# make a dict of ortho trios with chimp genes as key
ChimpOrthoTrios = {}
for gene in OrthoTrios:
    chimpgene, gorillagene = OrthoTrios[gene][0], OrthoTrios[gene][1]
    ChimpOrthoTrios[chimpgene] = [gene, gorillagene]

# infer young and old nesting events in human and chimp
HumanOld, HumanYoung = InferYoungOldNestingEvents(OrthoPairs, OrthoTrios, ChimpPairs[1], GorillaPairs[1], HumanPairs[1])
ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthos, ChimpOrthoTrios, HumanPairs[1], GorillaPairs[1], ChimpPairs[1])
# do some QC   
for pair in HumanOld:
    assert [OrthoPairs[pair[0]], OrthoPairs[pair[1]]] in ChimpOld    
    
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
ChimpOld = FilterGenePairsWithoutExpression(ChimpOld, ChimpExpression)
ChimpYoung = FilterGenePairsWithoutExpression(ChimpYoung, ChimpExpression)

HumanInferredPairs = [HumanOld, HumanYoung]
ChimpInferredPairs = [ChimpOld, ChimpYoung]

if Analysis == 'pairs':
    # compare expression divergence between human host and nested genes and their un-nested orthologs in chimp   
    # remove human pairs if orthologs are nested in chimp
    print(len(HumanYoung))    
    to_remove = []
    for pair in HumanYoung:
        if OrthoPairs[pair[0]] in NestedSets[1] or OrthoPairs[pair[1]] in NestedSets[1]:
            to_remove.append(pair)
    for pair in to_remove:
        HumanYoung.remove(pair)
    print(len(HumanYoung))    
    # remove human pairs if genes are not expressed
    to_remove = [pair for pair in HumanYoung if pair[0] not in HumanExpression or pair[1] not in HumanExpression]
    for pair in to_remove:
        HumanYoung.remove(pair)
    print(len(HumanYoung))
    # remove pairs if orthologs lack expression
    to_remove = []
    for pair in HumanYoung:
        if OrthoPairs[pair[0]] not in ChimpExpression or OrthoPairs[pair[1]] not in ChimpExpression:
            to_remove.append(pair)
    for pair in to_remove:
        HumanYoung.remove(pair)
    print(len(HumanYoung))
    # get the chimp un-nested orthologs
    ChimpUnested = []
    for pair in HumanYoung:
        ChimpUnested.append([OrthoPairs[pair[0]], OrthoPairs[pair[1]]])
    # compute expression divergence betwen human nested gene pairs
    HumanDiv = ComputeExpressionDivergenceGenePairs(HumanYoung, HumanExpression)
    ChimpDiv = ComputeExpressionDivergenceGenePairs(ChimpUnested, ChimpExpression)
    P = PermutationResampling(HumanDiv, ChimpDiv, 10000, statistic = np.mean)
    print(len(HumanDiv), len(ChimpDiv), np.mean(HumanDiv), np.mean(ChimpDiv), P)
    


elif Analysis == 'orthos':
    # compare expression divergence between young nested genes and their un-nested orthologs
    # compare expression divergence between young host genes and their un-nested orthologs

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
    
  




## use this function to create lists of orthologs with both genes expressed
#def ExpressedOrthologousPairs(Sp1Expression, Sp2Expression, Genes, Orthologs):
#    '''
#    (dict, dict, set, dict) -> list
#    Take the dictionaries of expression profiles for species 1 and 2, the set
#    of genes of interest in species 1, and the dictionary of orthologs and return
#    a list of expressed orthologous pairs
#    '''
#    # create a list of gene pairs
#    ExpressedOrthos = []
#    # loop over gene set of interest
#    for gene in Genes:
#        # check that gene has ortholog
#        if gene in Orthologs:
#            # check that gene and its orthologs are expressed
#            if gene in Sp1Expression and Orthologs[gene] in Sp2Expression:
#                ExpressedOrthos.append([gene, Orthologs[gene]])
#    return ExpressedOrthos













  
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
