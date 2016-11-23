# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 17:24:13 2016

@author: RJovelin
"""


# use this script to plot the number of introns, or the intron length between host,
# nested and un-nested genes

# usage python3 PlotIntronComparison.py [option]:


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
with open('OrangOutanHostNestedGenes.json') as orangoutan_json_data:
    OrangOutanHostGenes = json.load(orangoutan_json_data)
with open('MacaqueHostNestedGenes.json') as macaque_json_data:
    MacaqueHostGenes = json.load(macaque_json_data)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
PabGFF = 'Pongo_abelii.PPYG2.86.gff3'
MmlGFF = 'Macaca_mulatta.Mmul_8.0.1.86.gff3' 

# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF, PabGFF, MmlGFF]
# make a list of species names
SpeciesNames = ['human', 'chimp', 'gorilla', 'orangoutan', 'macaque']
# make a list of host:nested genes dictionaries
HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes, OrangOutanHostGenes, MacaqueHostGenes]

# make parallel lists to store expression specificity of host and nested genes [[human], [chimp], [gorilla], [orangutan], [macaque]]
HostSpecificity, NestedSpecificity = [], []

# loop over GFF files, find nested and intronic=nested genes in each species 
for i in range(len(GFFs)):
    print(GFFs[i][:GFFs[i].index('.')], SpeciesNames[i])
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    SpGeneChromoCoord = ChromoGenesCoord(GFFs[i])
    # map each gene to its mRNA transcripts
    SpMapGeneTranscript = GeneToTranscripts(GFFs[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    SpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(SpGeneChromoCoord, SpMapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    SpGeneCoord = FromChromoCoordToGeneCoord(SpGeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    SpOrderedGenes = OrderGenesAlongChromo(SpGeneChromoCoord)
    # get expression profile of the species genes
    SpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[i])
    # remove genes wuthout expression
    SpExpression = RemoveGenesLackingExpression(SpExpression)
    # get relative expression
    SpExpression = TransformRelativeExpression(SpExpression)
    # compute expression specificity for all genes
    SpTau = ExpressionSpecificity(SpExpression)    
    # generate un-nested gene pairs on ech chromosome with length < 500 bp
    SpControlPairs = GenerateUnNestedGenePairs(HostGenes[i], SpGeneCoord, SpOrderedGenes, SpExpression)
    # make a list of host-nested gene pairs
    SpHostNestedPairs = GetHostNestedPairs(HostGenes[i])
    # remove gene pairs with genes lacking expression
    SpHostNestedPairs = FilterGenePairsWithoutExpression(SpHostNestedPairs, SpExpression)
    # create lists to store specificity for host and nested genes for that species
    tauhost, taunested = [], []
    # loop over gene pairs:
    for pair in SpHostNestedPairs:
        # populate lists
        tauhost.append(SpTau[pair[0]])
        taunested.append(SpTau[pair[1]])
    HostSpecificity.append(tauhost)
    NestedSpecificity.append(taunested)


# make lists of gene specificity for each species     
HumanTau = [HostSpecificity[0], NestedSpecificity[0]]
ChimpTau = [HostSpecificity[1], NestedSpecificity[1]]
GorillaTau = [HostSpecificity[2], NestedSpecificity[2]]
OrangutanTau = [HostSpecificity[3], HostSpecificity[3]]
Macaque = [HostSpecificity[4], HostSpecificity[4]]

# make a list with all the list data
AllData = [HumanTau, ChimpTau, GorillaTau, OrangutanTau, MacaqueTau]

# create a function to get the mean and SEM of items in a list
def GetMeanSEM(L):
    '''
    (list) -> (list, list)
    Take a list of inner lists of numbers and return a list with mean values
    and a parallel list with SEM values for each item of the outter list
    '''
    # create lists of mean and SEM
    MeanVal, SEMVal = [], []
    # loop over the outter ist
    for i in range(len(L)):
        MeanVal.append(np.mean(L[i]))
        SEMVal.append(np.std(L[i]) / math.sqrt(len(L[i])))
    return MeanVal, SEMVal
    
# create lists with means and with SEM
HumanMeans, HumanSEM = GetMeanSEM(HumanTau)
ChimpMeans, ChimpSEM = GetMeanSEM(ChimpTau)
GorillaMeans, GorillaSEM = GetMeanSEM(GorillaTau)
OrangutanMeans, OrangutanSEM = GetMeanSEM(OrangutanTau)
MacaqueMeans, MacaqueSEM = GetMeanSEM(MacaqueTau)

# perform statistical tests between gene categories in all species
# create dict to store results {species: [P_host-nested, P_host-unnested, P_nested-unnested]}
PValues = {}
# loop over inner lists in data list
for i in range(len(AllData)):
    # initialize dict with empty list
    PValues[SpeciesNames[i]] = []
    # loop over inner list, compare gene categories
    for j in range(0, len(AllData[i]) -1):
        for k in range(j+1, len(AllData[i])):
            P = stats.ranksums(AllData[i][j], AllData[i][k])[1]
            PValues[SpeciesNames[i]].append(P)
# print p values
for sp in PValues:
    print(sp, PValues[sp])




#    # print number of genes on each chromo
#    chromohost = {}
#    expressedpairs = 0
#    for pair in SpHostNestedPairs:
#        if pair[0] in SpExpression and pair[1] in SpExpression:
#            chromo = SpGeneCoord[pair[0]][0]
#            if chromo in chromohost:
#                chromohost[chromo] += 1
#            else:
#                chromohost[chromo] = 1
#    for chromo in chromohost:
#        print(chromo, chromohost[chromo], len(SpControlPairs[chromo]))
        



