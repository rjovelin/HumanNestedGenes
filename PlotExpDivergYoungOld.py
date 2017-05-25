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
    
    
#    # make a list of GFF files
#    GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Mus_musculus.GRCm38.88.gff3',
#                 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
#                 'Gorilla_gorilla.gorGor3.1.88.gff3', 'Macaca_mulatta.Mmul_8.0.1.88.gff3',
#                 'Pongo_abelii.PPYG2.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
#                 'Bos_taurus.UMD3.1.88.gff3',
#                 'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
#                 'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
#                 'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
#                 'Monodelphis_domestica.BROADO5.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3',
#                 'Sorex_araneus.COMMON_SHREW1.88.gff3']
#
#    # make a parallel list of Species names
#    Species = ['Human', 'Mouse', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
#               'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
#               'Cat', 'Opossum', 'Platypus', 'Shrew']    
    
    
    
    
    
    # make a list of GFF files
    GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Mus_musculus.GRCm38.88.gff3',
                 'Bos_taurus.UMD3.1.88.gff3', 'Canis_familiaris.CanFam3.1.88.gff3',
                 'Choloepus_hoffmanni.choHof1.88.gff3', 'Dasypus_novemcinctus.Dasnov3.0.88.gff3',
                 'Equus_caballus.EquCab2.88.gff3', 'Erinaceus_europaeus.HEDGEHOG.88.gff3',
                 'Felis_catus.Felis_catus_6.2.88.gff3', 'Monodelphis_domestica.BROADO5.88.gff3',
                 'Ornithorhynchus_anatinus.OANA5.88.gff3', 'Sorex_araneus.COMMON_SHREW1.88.gff3']

    # make a parallel list of Species names
    Species = ['Human', 'Mouse', 'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse',
               'Hedgehog', 'Cat', 'Opossum', 'Platypus', 'Shrew']

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
# get nested pairs in other species (human is focal species)
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


# make a list of gene coordinates in human, chimp and mouse      
AllCoordinates, AllOrdered = [], []
# loop over GFF files
for i in range(len(GFF_Files)):
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF_Files[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF_Files[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    AllCoordinates.append(GeneCoord)
    AllOrdered.append(OrderedGenes)

# generate a set of nested genes in human and in sister species
HumanNestedGenes = MakeFullPartialOverlapGeneSet(AllNestedGenes[0])
SisterSpNestedGenes = MakeFullPartialOverlapGeneSet(AllNestedGenes[1])
print(len(HumanNestedGenes), len(SisterSpNestedGenes))



if Analysis == 'pairs':
    
    # compare expression divergence between gene pairs for young nested pairs in
    # human with their ancestral non-nested pairs in siste species
    
    # remove pairs for which both genes are not expressed
    YoungNested = FilterGenePairsWithoutExpression(YoungNested, HumanExpression, 'strict')
    OldNested = FilterGenePairsWithoutExpression(OldNested, HumanExpression, 'strict')
    print(len(OldNested), len(YoungNested))
    
    # generate gene pairs with orthologs of human young nested in sister species
    AncestralPairs = []
    for pair in YoungNested:
        for ortho1 in AllOrthologs[0][pair[0]]:
            for ortho2 in AllOrthologs[0][pair[1]]:
                AncestralPairs.append([ortho1, ortho2])
                assert [ortho1, ortho2] not in SpeciesNestedPairs[0]
    # remove pairs for which both members are not expressed in sister species            
    SisterSpAncestralPairs = FilterGenePairsWithoutExpression(AncestralPairs, SisterSpExpression, 'strict')
    print(len(SisterSpAncestralPairs))    
    # remove pairs for which members are not valid genes (eg pseudogenes, etc)        
    to_remove = [pair for pair in SisterSpAncestralPairs if pair[0] not in AllCoordinates[1] or pair[1] not in AllCoordinates[1]]
    for pair in to_remove:
        SisterSpAncestralPairs.remove(pair)        
        
    # remove pairs if any gene is nested    
    to_remove = [pair for pair in SisterSpAncestralPairs if pair[0] in SisterSpNestedGenes or pair[1] in SisterSpNestedGenes]
    for pair in to_remove:
        SisterSpAncestralPairs.remove(pair)
    print(len(SisterSpAncestralPairs))
    
    # generate a list of control un-nested pairs
    SisterSpRandomGenes = GenerateAllUnNestedGenes(SisterSpNestedGenes, AllOrdered[1], SisterSpExpression)
    # make a list of control un-nested pairs in sister species
    SisterSpControlPairs = []
    to_remove = []
    for pair in SisterSpAncestralPairs:
        # make a list of matching gene pairs (orientation, chromosome, distance)
        PairPool = GenerateMatchingPoolPairs(pair, SisterSpRandomGenes, AllCoordinates[1], 2000)
        # draw a matching gene pair at random
        i = random.randint(0, len(PairPool) -1)
        SisterSpControlPairs.append(PairPool[i])
       
    # compute divergence in young human nested pairs
    HumanYoungDiv = ComputeExpressionDivergenceGenePairs(YoungNested, HumanExpression)    
    # compute divergence in ancestral un-nested pairs
    SisterSpAncestralDiv = ComputeExpressionDivergenceGenePairs(SisterSpAncestralPairs, SisterSpExpression)
    # compute expression divergence between un-nested control pairs in sister species    
    SisterSpControlDiv = ComputeExpressionDivergenceGenePairs(SisterSpControlPairs, SisterSpExpression)

    # compute P values using permutation tests
    P = PermutationResampling(HumanYoungDiv, SisterSpAncestralDiv, 1000, statistic = np.mean)
    print(len(HumanYoungDiv), len(SisterSpAncestralDiv), np.mean(HumanYoungDiv), np.mean(SisterSpAncestralDiv), P)
    P = PermutationResampling(SisterSpAncestralDiv, SisterSpControlDiv , 1000, statistic = np.mean)
    print(len(SisterSpAncestralDiv), len(SisterSpControlDiv), np.mean(SisterSpAncestralDiv), np.mean(SisterSpControlDiv), P)
    

elif Analysis == 'orthos':
    # compare distances between expression profiles of internal/external-like genes and their orthologs
    # young internal genes and their un-nested orthologs
    # old internal genes and their nested orthologs
    # young external genes and their un-nested orthologs
    # old external genes and their nested orthologs


    # create a list of young external genes in human and sister species
    HumanYoungExt = list(set([pair[0] for pair in YoungNested if pair[0] in HumanExpression]))
    # create a list of young internal genes in human and sister species    
    HumanYoungInt = list(set([pair[1] for pair in YoungNested if pair[1] in HumanExpression]))
    # create a list of old external genes in human and sister species
    HumanOldExt = list(set([pair[0] for pair in OldNested if pair[0] in HumanExpression]))
    # create a list of old internal genes in human and sister species
    HumanOldInt = list(set([pair[1] for pair in OldNested if pair[1] in HumanExpression]))
    
    # make pairs of orthologs
    YoungExtPairs, YoungIntPairs, OldExtPairs, OldIntPairs = [], [], [], []
    for gene in HumanYoungExt:
        for ortho in AllOrthologs[0][gene]:
            # check that ortholog is expressed, is valid and is not nested
            if ortho in SisterSpExpression and ortho in AllCoordinates[1] and ortho not in SisterSpNestedGenes:
                YoungExtPairs.append([gene, ortho])
    for gene in HumanYoungInt:
        for ortho in AllOrthologs[0][gene]:
            # check that ortholog is expressed, is valid and is not nested
            if ortho in SisterSpExpression and ortho in AllCoordinates[1] and ortho not in SisterSpNestedGenes:
                YoungIntPairs.append([gene, ortho])
    for gene in HumanOldExt:
        for ortho in AllOrthologs[0][gene]:
            # check that ortholog is expressed, is valid and is nested
            if ortho in SisterSpExpression and ortho in AllCoordinates[1] and ortho in SisterSpNestedGenes:
                OldExtPairs.append([gene, ortho])
    for gene in HumanOldInt:
        for ortho in AllOrthologs[0][gene]:
            # check that ortholog is expressed, is valid and is nested
            if ortho in SisterSpExpression and ortho in AllCoordinates[1] and ortho in SisterSpNestedGenes:
                OldIntPairs.append([gene, ortho])
        
    # generate a dict to draw genes in human    
    HumanRandomGenes = GenerateAllUnNestedGenes(HumanNestedGenes, AllOrdered[0], HumanExpression)

    # make lists of human genes
    HumanExternal = list(set([pair[0] for pair in YoungExtPairs]))
    HumanInternal = list(set([pair[0] for pair in YoungIntPairs]))
    
    # generate lists of control genes, match genes by chromosome and tissue specificity
    HumanExtLike = GenerateMatchingGenes(HumanExternal, AllCoordinates[0], HumanRandomGenes, HumanSpecificity, AllOrthologs[0], SisterSpExpression)
    HumanIntLike = GenerateMatchingGenes(HumanInternal, AllCoordinates[0], HumanRandomGenes, HumanSpecificity, AllOrthologs[0], SisterSpExpression)
    
    print(len(HumanExtLike))
    print(len(HumanIntLike))
    
    # make list of gene pairs
    HumanControlExt = []
    for gene in HumanExtLike:
        assert gene in HumanExpression
        assert gene in AllOrthologs[0]
        for ortho in AllOrthologs[0][gene]:
            if ortho in SisterSpExpression:
                HumanControlExt.append([gene, ortho])
    HumanControlInt = []
    for gene in HumanIntLike:
        assert gene in HumanExpression
        assert gene in AllOrthologs[0]
        for ortho in AllOrthologs[0][gene]:
            if ortho in SisterSpExpression:
                HumanControlInt.append([gene, ortho])
       
    # compute expression divergence between orthologs
    HumanControlExtDiv = ComputeExpressionDivergenceOrthologs(HumanControlExt, HumanExpression, SisterSpExpression)
    HumanControlIntDiv = ComputeExpressionDivergenceOrthologs(HumanControlInt, HumanExpression, SisterSpExpression)    
    YoungExtDiv = ComputeExpressionDivergenceOrthologs(YoungExtPairs, HumanExpression, SisterSpExpression)
    YoungIntDiv = ComputeExpressionDivergenceOrthologs(YoungIntPairs, HumanExpression, SisterSpExpression)
    OldExtDiv = ComputeExpressionDivergenceOrthologs(OldExtPairs, HumanExpression, SisterSpExpression)
    OldIntDiv = ComputeExpressionDivergenceOrthologs(OldIntPairs, HumanExpression, SisterSpExpression)
    
    # compare expression divergence among genes
    P = PermutationResampling(YoungExtDiv, HumanControlExtDiv, 1000, statistic = np.mean)
    print(len(YoungExtDiv), len(HumanControlExtDiv), np.mean(YoungExtDiv), np.mean(HumanControlExtDiv), P)
    P = PermutationResampling(YoungIntDiv, HumanControlIntDiv, 1000, statistic = np.mean)
    print(len(YoungIntDiv), len(HumanControlIntDiv), np.mean(YoungIntDiv), np.mean(HumanControlIntDiv), P)
    P = PermutationResampling(YoungExtDiv, YoungIntDiv, 1000, statistic = np.mean)
    print(len(YoungExtDiv), len(YoungIntDiv), np.mean(YoungExtDiv), np.mean(YoungIntDiv), P)

#    P = PermutationResampling(OldExtDiv, HumanControlExtDiv, 1000, statistic = np.median)    
#    print(len(OldExtDiv), len(HumanControlExtDiv), np.median(OldExtDiv), np.median(HumanControlExtDiv), P)
#    P = PermutationResampling(OldExtDiv, YoungExtDiv, 1000, statistic = np.median)    
#    print(len(OldExtDiv), len(YoungExtDiv), np.median(OldExtDiv), np.median(YoungExtDiv), P)
#    
#    
#    P = PermutationResampling(OldIntDiv, HumanControlIntDiv, 1000, statistic = np.median)    
#    print(len(OldIntDiv), len(HumanControlIntDiv), np.median(OldIntDiv), np.median(HumanControlIntDiv), P)
#    P = PermutationResampling(OldIntDiv, YoungIntDiv, 1000, statistic = np.median)    
#    print(len(OldIntDiv), len(YoungIntDiv), np.median(OldIntDiv), np.median(YoungIntDiv), P)
    






    