# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""


# use this script to plot expression divergence before and after nesting

# usage python3 PlotExpDivergBeforeAfterNesting.py 


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


# get expression profile of the species genes
HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', SpeciesNames[0])
# remove genes wuthout expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
# make a set of host and nested genes (include all host and nested even if not expressed)    
HumanNestedConformation = MakeHostNestedGeneSet(HostGenes[0])    
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
ChimpNestedConformation = MakeHostNestedGeneSet(HostGenes[1])    
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
GorillaNestedConformation = MakeHostNestedGeneSet(HostGenes[2])    
# make a list of host-nested gene pairs
GorillaHostNestedPairs = GetHostNestedPairs(HostGenes[2])
# remove gene pairs with genes lacking expression
GorillaHostNestedPairs = FilterGenePairsWithoutExpression(GorillaHostNestedPairs, GorillaExpression)


HumanOrthologs = ParseOrthologFile('HumanChimpGorillaOrthologs.txt')
# create a dicti with chimp genes as key {chimp_gene: [human_orto, gorilla_ortho]}
ChimpOrthologs = {}
for gene in HumanOrthologs:
    ChimpOrthologs[HumanOrthologs[gene][0]] = [gene, HumanOrthologs[gene][1]]

# make lists of old and yound host:nested pairs
HumanOld, HumanYoung = SortYoungOldNestingEvents(HumanOrthologs, HostGenes[0], HostGenes[1], HostGenes[2], HumanHostNestedPairs)
ChimpOld, ChimpYoung = SortYoungOldNestingEvents(ChimpOrthologs, HostGenes[1], HostGenes[0], HostGenes[2], ChimpHostNestedPairs)
 

# make list of ancestral un-nested gene pairs
humanancestral, chimpancestral = [], []

for pair in HumanYoung:
    ortho1, ortho2 = HumanOrthologs[pair[0]][0], HumanOrthologs[pair[1]][0]
    assert ortho1 in ChimpExpression
    assert ortho2 in ChimpExpression
    chimpancestral.append([ortho1, ortho2])

for pair in ChimpYoung:
    ortho1, ortho2 = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
    assert ortho1 in HumanExpression
    assert ortho2 in HumanExpression
    humanancestral.append([ortho1, ortho2])


print(len(HumanOld), len(HumanYoung), len(humanancestral))
print(len(ChimpOld), len(ChimpYoung), len(chimpancestral))


## compute divergence
#crmderivedDiv = ComputeExpressionDivergenceGenePairs(crmyoung, CrmRelativeExpression)
#cbrancestral = ComputeExpressionDivergenceGenePairs(cbrancestral, CbrRelativeExpression)
#print(len(crmderivedDiv), np.mean(crmderivedDiv), len(cbrancestral), np.mean(cbrancestral), stats.ranksums(crmderivedDiv, cbrancestral)[1])
#    
#for pair in cbryoung:
#    hostorthosp2 = CbrOrthos[pair[0]][1]
#    nestedorthosp2 = CbrOrthos[pair[1]][1]
#    crmancestral.append([hostorthosp2, nestedorthosp2])
#    
#    
#cbrderivedDiv = ComputeExpressionDivergenceGenePairs(cbryoung, CbrRelativeExpression)
#crmancestralDiv = ComputeExpressionDivergenceGenePairs(crmancestral, CrmRelativeExpression)
#print(len(cbrderivedDiv), np.mean(cbrderivedDiv), len(crmancestralDiv), np.mean(crmancestralDiv), stats.ranksums(cbrderivedDiv, crmancestralDiv))