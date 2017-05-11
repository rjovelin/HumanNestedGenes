# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""

#############

# modify the function to infer young/old nesting events
# include multiple outgroups
# for human-chimp: outgroups = baboon, gorilla, macaque, marmoset, orangoutan, mouse etc
# for human-mouse: outroups = dog, cat/cow, bat, shrews, hedgehog, marsupial, platypus
# get pairwise orthologs
# require presence of both genes in human and chimp or mouse and at least 1 outgroup
# if gene pair not present in chimp/mouse and not present in all outgroup then young nesting event
# how to deal with one_to_many and many_to_many orthologs?
# include all orthology type, if 1_to_many or many_to_many, consider nested if any ortholog is nested (at least 1 ortho of external and 1 ortho of internal must be nested)

# http://tolweb.org/Eutheria/15997
# http://www.pnas.org/content/suppl/2007/12/14/0705658104.DC1
# https://research.amnh.org/paleontology/perissodactyl/node/55
# http://useast.ensembl.org/info/genome/compara/homology_method.html

##############






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


if SisterSp == 'chimp':
    # make a list of GFF files
    GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
                 'Gorilla_gorilla.gorGor3.1.88.gff3', 'Macaca_mulatta.Mmul_8.0.1.88.gff3',
                 'Pongo_abelii.PPYG2.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
                 'Mus_musculus.GRCm38.88.gff3', 'Bos_taurus.UMD3.1.88.gff3',
                 'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
                 'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
                 'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
                 'Monodelphis_domestica.BROADO5.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3',
                 'Sorex_araneus.COMMON_SHREW1.88.gff3']

    # make a parallel list of Species names
    Species = ['Human', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
               'Mouse', 'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
               'Cat', 'Opossum', 'Platypus', 'Shrew']
elif SisterSp == 'mouse':
    
    
    # make a list of GFF files
    GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Mus_musculus.GRCm38.88.gff3',
                 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
                 'Gorilla_gorilla.gorGor3.1.88.gff3', 'Macaca_mulatta.Mmul_8.0.1.88.gff3',
                 'Pongo_abelii.PPYG2.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
                 'Bos_taurus.UMD3.1.88.gff3',
                 'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
                 'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
                 'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
                 'Monodelphis_domestica.BROADO5.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3',
                 'Sorex_araneus.COMMON_SHREW1.88.gff3']

    # make a parallel list of Species names
    Species = ['Human', 'Mouse', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
               'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
               'Cat', 'Opossum', 'Platypus', 'Shrew']    
    
    
    
    
    
#    # make a list of GFF files
#    GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Mus_musculus.GRCm38.88.gff3',
#                 'Bos_taurus.UMD3.1.88.gff3', 'Canis_familiaris.CanFam3.1.88.gff3',
#                 'Choloepus_hoffmanni.choHof1.88.gff3', 'Dasypus_novemcinctus.Dasnov3.0.88.gff3',
#                 'Equus_caballus.EquCab2.88.gff3', 'Erinaceus_europaeus.HEDGEHOG.88.gff3',
#                 'Felis_catus.Felis_catus_6.2.88.gff3', 'Monodelphis_domestica.BROADO5.88.gff3',
#                 'Ornithorhynchus_anatinus.OANA5.88.gff3', 'Sorex_araneus.COMMON_SHREW1.88.gff3']
#
#    # make a parallel list of Species names
#    Species = ['Human', 'Mouse', 'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse',
#               'Hedgehog', 'Cat', 'Opossum', 'Platypus', 'Shrew']

# make a parallel list of json files with nested genes
JsonFiles = [i + 'NestedGenes.json' for i in Species]
# make a parallel list of ortholog files
OrthoFiles = ['Human' + i + 'Orthologs.txt' for i in Species[1:]]

# make a list of dictionaries
AllNestedGenes = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    nested = json.load(json_data)
    json_data.close()
    AllNestedGenes.append(nested)
# make a list of dicts with orthologs
AllOrthologs = []
for i in range(len(OrthoFiles)):
    orthologs = MatchOrthologs(OrthoFiles[i])
    AllOrthologs.append(orthologs)
    

# get nested pairs in human and other species (sister species is 1st in list)
HumanNestedPairs = GetHostNestedPairs(AllNestedGenes[0])
SpeciesNestedPairs = [GetHostNestedPairs(AllNestedGenes[i]) for i in range(1, len(AllNestedGenes))]

# infer old and young nested events in human
OldNested, YoungNested = InferYoungOldNestingEvents(HumanNestedPairs, SpeciesNestedPairs, AllOrthologs) 
print(len(OldNested), len(YoungNested))


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





############################################################



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

# do some QC   
a = [set([OrthoPairs[pair[0]], OrthoPairs[pair[1]]]) for pair in HumanOld]
b = [set(pair) for pair in SisterSpOld]
assert len(a) == len(b)
for i in a:
    assert i in b    


HumanInferredPairs = [HumanOld, HumanYoung]



if Analysis == 'pairs':
    
    
    # count how many ortho pairs are nested in sister species


    
    
    
    
    # compare expression divergence between human host and nested genes and their un-nested orthologs in sister-species   
    # remove human pairs if orthologs are nested in sister-species
    to_remove = [pair for pair in HumanYoung if OrthoPairs[pair[0]] in NestedSets[1] or OrthoPairs[pair[1]] in NestedSets[1]]
    for pair in to_remove:
        HumanYoung.remove(pair)
    # remove pairs if orthologs are nested in outgroup
    to_remove = []
    for pair in HumanYoung:
        # check if ortholog is present in outgroup
        if pair[0] in OrthoTrios:
            # add pair to list of pairs to remove if pair not already added
            if OrthoTrios[pair[0]][1] in NestedSets[2] and pair not in to_remove:
                to_remove.append(pair)
        if pair[1] in OrthoTrios:
            if OrthoTrios[pair[1]][1] in NestedSets[2] and pair not in to_remove:
                to_remove.append(pair)
    print(len(to_remove))
    if len(to_remove) != 0:
        for pair in to_remove:
            HumanYoung.remove(pair)
    to_remove = []
    for pair in SisterSpYoung:
        # check if ortholog is present in outgroup
        if pair[0] in SisterOrthoTrios:
            # add pair to list of pairs to remove if pair not already added
            if SisterOrthoTrios[pair[0]][1] in NestedSets[2] and pair not in to_remove:
                to_remove.append(pair)
        if pair[1] in SisterOrthoTrios:
            if SisterOrthoTrios[pair[1]][1] in NestedSets[2] and pair not in to_remove:
                to_remove.append(pair)
    print(len(to_remove))
    if len(to_remove) != 0:
        for pair in to_remove:
            SisterSpYoung.remove(pair)
    
    # remove human pairs if genes are not expressed
    to_remove = [pair for pair in HumanYoung if pair[0] not in HumanExpression or pair[1] not in HumanExpression]
    for pair in to_remove:
        HumanYoung.remove(pair)
    # remove sister species pairs if genes are not expressed
    to_remove = [pair for pair in SisterSpYoung if pair[0] not in SisterSpExpression or pair[1] not in SisterSpExpression]
    for pair in to_remove:
        SisterSpYoung.remove(pair)
    
    # get the sister-species un-nested orthologs
    SisterSpUnested = []
    for pair in HumanYoung:
        # check if orthologs are expressed in sister species
        if OrthoPairs[pair[0]] in SisterSpExpression and OrthoPairs[pair[1]] in SisterSpExpression:
            SisterSpUnested.append([OrthoPairs[pair[0]], OrthoPairs[pair[1]]])
    # get the human un-nested orthologs
    HumanUnNested = []
    for pair in SisterSpYoung:
        # check that orthologs are expressed in human
        if SisterOrthos[pair[0]] in HumanExpression and SisterOrthos[pair[1]] in HumanExpression:
            HumanUnNested.append([SisterOrthos[pair[0]], SisterOrthos[pair[1]]])
    
    # generate a dict to draw random genes in sister-species
    SisterRandomGenes = GenerateAllUnNestedGenes(NestedSets[1], AllOrdered[1], SisterSpExpression)
    # generate a dict to draw genes in human    
    HumanRandomGenes = GenerateAllUnNestedGenes(NestedSets[0], AllOrdered[0], HumanExpression)
    
    # make a list of control un-nested pairs in sister species
    SisterSpControlPairs = []
    for pair in SisterSpUnested:
        # make a list of matching gene pairs (orientation, chromosome, distance)
        PairPool = GenerateMatchingPoolPairs(pair, SisterRandomGenes, AllCoordinates[1], 2000)
        # draw a matching gene pair at random
        i = random.randint(0, len(PairPool) -1)
        SisterSpControlPairs.append(PairPool[i])
    # make a list of control un-nested pairs in human
    HumanControlPairs = []    
    for pair in HumanUnNested:
        # make a list of matching gene pairs (orientation, chromosome, distance)
        PairPool = GenerateMatchingPoolPairs(pair, HumanRandomGenes, AllCoordinates[0], 2000)
        # draw a matching gene pair at random
        i = random.randint(0, len(PairPool) -1)
        HumanControlPairs.append(PairPool[i])

    # compute expression divergence betwen human nested gene pairs
    HumanNestedDiv = ComputeExpressionDivergenceGenePairs(HumanYoung, HumanExpression)
    # compute expression divergence betweennested gene pairs in sister species
    SisterSpNestedDiv = ComputeExpressionDivergenceGenePairs(SisterSpYoung, SisterSpExpression)
    # compute expression divergence between un-nested gene pairs in human
    HumanUnNestedDiv = ComputeExpressionDivergenceGenePairs(HumanUnNested, HumanExpression)
    # compute expression divergence between un-nested gene pairs in sister species    
    SisterSpUnNestedDiv = ComputeExpressionDivergenceGenePairs(SisterSpUnested, SisterSpExpression)
    # compute expression divergence between un-nested control pairs in human
    HumanControlDiv = ComputeExpressionDivergenceGenePairs(HumanControlPairs, HumanExpression)
    # compute expression divergence between un-nested control pairs in sister species    
    SisterSpControlDiv = ComputeExpressionDivergenceGenePairs(SisterSpControlPairs, SisterSpExpression)

    # merge gene expression divergence for the nested genes in human and sister species
    ExpDivNested = []
    ExpDivNested.extend(HumanNestedDiv)
    ExpDivNested.extend(SisterSpNestedDiv)
    # merge gene expression divergence for the un-nested genes in human and sister species
    ExpDivUnNested = []
    ExpDivUnNested.extend(HumanUnNestedDiv)
    ExpDivUnNested.extend(SisterSpUnNestedDiv)
    # merge gene expression divergence fot the control pairs
    ExpDivControl = []
    ExpDivControl.extend(HumanControlDiv)
    ExpDivControl.extend(SisterSpControlDiv)

    P = PermutationResampling(ExpDivNested, ExpDivUnNested, 1000, statistic = np.mean)
    print(len(ExpDivNested), len(ExpDivUnNested), np.mean(ExpDivNested), np.mean(ExpDivUnNested), P)
    P = PermutationResampling(ExpDivNested, ExpDivControl, 1000, statistic = np.mean)
    print(len(ExpDivNested), len(ExpDivControl), np.mean(ExpDivNested), np.mean(ExpDivControl), P)
    P = PermutationResampling(ExpDivUnNested, ExpDivControl, 1000, statistic = np.mean)
    print(len(ExpDivUnNested), len(ExpDivControl), np.mean(ExpDivUnNested), np.mean(ExpDivControl), P)


elif Analysis == 'orthos':
    # compare distances between expression profiles of internal/external-like genes and their orthologs
    # young internal genes and their un-nested orthologs
    # old internal genes and their nested orthologs
    # young external genes and their un-nested orthologs
    # old external genes and their un-nested orthologs

    # create a list of young external genes in human and sister species
    HumanYoungExt = list(set([pair[0] for pair in HumanYoung if pair[0] in HumanExpression]))
    SisterSpYoungExt = list(set([pair[0] for pair in SisterSpYoung if pair[0] in SisterSpExpression]))   
    # create a list of young internal genes in human and sister species    
    HumanYoungInt = list(set([pair[1] for pair in HumanYoung if pair[1] in HumanExpression]))
    SisterSpYoungInt = list(set([pair[1] for pair in SisterSpYoung if pair[1] in SisterSpExpression]))    
    # create a list of old external genes in human and sister species
    HumanOldExt = list(set([pair[0] for pair in HumanOld if pair[0] in HumanExpression]))
    SisterSpOldExt = list(set([pair[0] for pair in SisterSpOld if pair[0] in SisterSpExpression]))
    # create a list of old internal genes in human and sister species
    HumanOldInt = list(set([pair[1] for pair in HumanOld if pair[1] in HumanExpression]))
    SisterSpOldInt = list(set([pair[1] for pair in SisterSpOld if pair[1] in SisterSpExpression]))

    # remove young external and internal genes if their ortholog is nested
    to_remove = [gene for gene in HumanYoungExt if OrthoPairs[gene] in NestedSets[1]]
    for gene in to_remove:
        HumanYoungExt.remove(gene)
    to_remove = [gene for gene in HumanYoungInt if OrthoPairs[gene] in NestedSets[1]]
    for gene in to_remove:
        HumanYoungInt.remove(gene)
    to_remove = [gene for gene in SisterSpYoungExt if SisterOrthos[gene] in NestedSets[0]]    
    for gene in to_remove:
        SisterSpYoungExt.remove(gene)
    to_remove = [gene for gene in SisterSpYoungInt if SisterOrthos[gene] in NestedSets[0]]
    for gene in to_remove:
        SisterSpYoungInt.remove(gene)
    
    # remove young external and internal genes if their ortholog is nested in outgroup    
    to_remove = [gene for gene in HumanYoungExt if gene in OrthoTrios and OrthoTrios[gene][1] in NestedSets[2]]
    for gene in to_remove:
        HumanYoungExt.remove(gene)
    to_remove = [gene for gene in HumanYoungInt if gene in OrthoTrios and OrthoTrios[gene][1] in NestedSets[2]]
    for gene in to_remove:
        HumanYoungInt.remove(gene)
    to_remove = [gene for gene in SisterSpYoungExt if gene in SisterOrthoTrios and SisterOrthoTrios[gene][1] in NestedSets[2]]
    for gene in to_remove:
        SisterSpYoungExt.remove(gene)
    to_remove = [gene for gene in SisterSpYoungInt if gene in SisterOrthoTrios and SisterOrthoTrios[gene][1] in NestedSets[2]]
    for gene in to_remove:
        SisterSpYoungInt.remove(gene)
    
    # remove genes if their ortholog is not expressed
    to_remove = [gene for gene in HumanYoungExt if OrthoPairs[gene] not in SisterSpExpression]    
    for gene in to_remove:
        HumanYoungExt.remove(gene)
    to_remove = [gene for gene in HumanYoungInt if OrthoPairs[gene] not in SisterSpExpression]
    for gene in to_remove:
        HumanYoungInt.remove(gene)
    to_remove = [gene for gene in SisterSpYoungExt if SisterOrthos[gene] not in HumanExpression]    
    for gene in to_remove:
        SisterSpYoungExt.remove(gene)
    to_remove = [gene for gene in SisterSpYoungInt if SisterOrthos[gene] not in HumanExpression]
    for gene in to_remove:
        SisterSpYoungInt.remove(gene)
    
    # check that orthologs of old external and internal genes are nested    
    for gene in HumanOldExt:
        assert OrthoPairs[gene] in NestedSets[1]
    for gene in HumanOldInt:
        assert OrthoPairs[gene] in NestedSets[1]
    for gene in SisterSpOldExt:
        assert SisterOrthos[gene] in NestedSets[0]
    for gene in SisterSpOldInt:
        assert SisterOrthos[gene] in NestedSets[0]
    

    
    
    # generate a dict to draw random genes in sister-species
    SisterRandomGenes = GenerateAllUnNestedGenes(NestedSets[1], AllOrdered[1], SisterSpExpression)
    # generate a dict to draw genes in human    
    HumanRandomGenes = GenerateAllUnNestedGenes(NestedSets[0], AllOrdered[0], HumanExpression)

    # make list of control genes, match genes by chromosome and tissue specificity
    HumanExtLike = GenerateMatchingGenes(HumanYoungExt, AllCoordinates[0], HumanRandomGenes, HumanSpecificity, OrthoPairs, SisterSpExpression)
    HumanIntLike = GenerateMatchingGenes(HumanYoungInt, AllCoordinates[0], HumanRandomGenes, HumanSpecificity, OrthoPairs, SisterSpExpression)
    SisterSpExtLike = GenerateMatchingGenes(SisterSpYoungExt, AllCoordinates[1], SisterRandomGenes, SisterSpSpecificity, SisterOrthos, HumanExpression)    
    SisterSpYoungInt = GenerateMatchingGenes(SisterSpYoungInt, AllCoordinates[1], SisterRandomGenes, SisterSpSpecificity, SisterOrthos, HumanExpression)
    
    print(len(HumanExtLike))
    print(len(HumanIntLike))
    print(len(SisterSpExtLike))
    print(len(SisterSpYoungInt))
    
 
    # make list of gene pairs
    HumanControlExt = [[gene, OrthoPairs[gene]] for gene in HumanExtLike if gene in HumanExpression and OrthoPairs[gene] in SisterSpExpression] 
    HumanControlInt = [[gene, OrthoPairs[gene]] for gene in HumanIntLike if gene in HumanExpression and OrthoPairs[gene] in SisterSpExpression]

    


    YoungExtPairs = [[gene, OrthoPairs[gene]] for gene in HumanYoungExt if gene in HumanExpression and OrthoPairs[gene] in SisterSpExpression] 
    YoungIntPairs = [[gene, OrthoPairs[gene]] for gene in HumanYoungInt if gene in HumanExpression and OrthoPairs[gene] in SisterSpExpression]
    OldExtPairs = [[gene, OrthoPairs[gene]] for gene in HumanOldExt if gene in HumanExpression and OrthoPairs[gene] in SisterSpExpression] 
    OldIntPairs = [[gene, OrthoPairs[gene]] for gene in HumanOldInt if gene in HumanExpression and OrthoPairs[gene] in SisterSpExpression]

    HumanControlExtDiv = ComputeExpressionDivergenceOrthologs(HumanControlExt, HumanExpression, SisterSpExpression)
    HumanControlIntDiv = ComputeExpressionDivergenceOrthologs(HumanControlInt, HumanExpression, SisterSpExpression)    
    YoungExtDiv = ComputeExpressionDivergenceOrthologs(YoungExtPairs, HumanExpression, SisterSpExpression)
    YoungIntDiv = ComputeExpressionDivergenceOrthologs(YoungIntPairs, HumanExpression, SisterSpExpression)
    OldExtDiv = ComputeExpressionDivergenceOrthologs(OldExtPairs, HumanExpression, SisterSpExpression)
    OldIntDiv = ComputeExpressionDivergenceOrthologs(OldIntPairs, HumanExpression, SisterSpExpression)


##### check with median

    

    P = PermutationResampling(YoungExtDiv, HumanControlExtDiv, 1000, statistic = np.median)
    print(len(YoungExtDiv), len(HumanControlExtDiv), np.median(YoungExtDiv), np.median(HumanControlExtDiv), P)
    P = PermutationResampling(OldExtDiv, HumanControlExtDiv, 1000, statistic = np.median)    
    print(len(OldExtDiv), len(HumanControlExtDiv), np.median(OldExtDiv), np.median(HumanControlExtDiv), P)
    P = PermutationResampling(OldExtDiv, YoungExtDiv, 1000, statistic = np.median)    
    print(len(OldExtDiv), len(YoungExtDiv), np.median(OldExtDiv), np.median(YoungExtDiv), P)
    
    
    P = PermutationResampling(YoungIntDiv, HumanControlIntDiv, 1000, statistic = np.median)
    print(len(YoungIntDiv), len(HumanControlIntDiv), np.median(YoungIntDiv), np.median(HumanControlIntDiv), P)
    P = PermutationResampling(OldIntDiv, HumanControlIntDiv, 1000, statistic = np.median)    
    print(len(OldIntDiv), len(HumanControlIntDiv), np.median(OldIntDiv), np.median(HumanControlIntDiv), P)
    P = PermutationResampling(OldIntDiv, YoungIntDiv, 1000, statistic = np.median)    
    print(len(OldIntDiv), len(YoungIntDiv), np.median(OldIntDiv), np.median(YoungIntDiv), P)
    






    
#    # remove gene "duplicates" by removing chimp genes with ortologs already present in each group
#    for i in range(len(ChimpExtIntGenes)):
#        to_remove = []
#        for gene in ChimpExtIntGenes[i]:
#            if ChimpOrthos[gene] in HumanExtIntGenes[i]:
#                to_remove.append(gene)
#        for gene in to_remove:
#            ChimpExtIntGenes[i].remove(gene)


