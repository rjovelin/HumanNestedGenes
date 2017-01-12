# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 20:59:27 2017

@author: Richard
"""

# use this script to plot expression divergence between orthologs 

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
import matplotlib.gridspec as gridspec
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


# load dictionary of overlapping gene pairs
json_data = open('HumanOverlappingGenes.json')
Overlapping = json.load(json_data)
json_data.close()
# load dictionary of nested gene pairs
json_data = open('HumanNestedGenes.json')
Nested = json.load(json_data)
json_data.close()
# load dictionary of pibbyback gene pairs
json_data = open('HumanPiggyBackGenes.json')
Piggyback = json.load(json_data)
json_data.close()
# load dictionary of convergent gene pairs
json_data = open('HumanConvergentGenes.json')
Convergent = json.load(json_data)
json_data.close()
# load dictionary of divergent gene pairs
json_data = open('HumanDivergentGenes.json')
Divergent = json.load(json_data)
json_data.close()

# get GFF file
GFF = 'Homo_sapiens.GRCh38.86.gff3'
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)

# generate gene sets
NestedGenes  = MakeFullPartialOverlapGeneSet(Nested)
OverlappingGenes = MakeFullPartialOverlapGeneSet(Overlapping)
ConvergentGenes = MakeFullPartialOverlapGeneSet(Convergent)
DivergentGenes = MakeFullPartialOverlapGeneSet(Divergent)
PiggyBackGenes = MakeFullPartialOverlapGeneSet(Piggyback)

# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)

# create lists of nested gene pairs with same and opposite directions
NestedPairs = GetHostNestedPairs(Nested)
same, opposite = [], []
for pair in NestedPairs:
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        opposite.append(pair)
    elif len(set(orientation)) == 1:
        same.append(pair)
# create sets of internal and external nested genes depending on orientation 
InternalSameGenes, InternalOppositeGenes, ExternalSameGenes, ExternalOppositeGenes = set(), set(), set(), set()
for pair in same:
    ExternalSameGenes.add(pair[0])
    InternalSameGenes.add(pair[1])
for pair in opposite:
    ExternalOppositeGenes.add(pair[0])
    InternalOppositeGenes.add(pair[1])


# get expression profile of human genes
HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
# remove genes wuthout expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)

# get expression profile of chimp genes
ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
# remove genes wuthout expression
ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
# get relative expression
ChimpExpression = TransformRelativeExpression(ChimpExpression)

# get 1:1 orthologs between human and chimp
Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
 
# generate lists of ortholog pairs for each gene category

# use this function to create lists of orthologs with both genes expressed
def ExpressedOrthologousPairs(Sp1Expression, Sp2Expression, Genes, Orthologs):
    '''
    (dict, dict, set, dict) -> list
    Take the dictionaries of expression profiles for species 1 and 2, the set
    of genes of interest in species 1, and the dictionary of orthologs and return
    a list of expressed orthologous pairs
    '''
    # create a list of gene pairs
    ExpressedOrthos = []
    # loop over gene set of interest
    for gene in Genes:
        # check that gene has ortholog
        if gene in Orthologs:
            # check that gene and its orthologs are expressed
            if gene in Sp1Expression and Orthologs[gene] in Sp2Expression:
                ExpressedOrthos.append([gene, Orthologs[gene]])
    return ExpressedOrthos
    

NstOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, NestedGenes, Orthos)
OvlOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, OverlappingGenes, Orthos)
ConOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, ConvergentGenes, Orthos)
DivOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, DivergentGenes, Orthos) 
PbkOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, PiggyBackGenes, Orthos)  
CisIntOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, InternalSameGenes, Orthos) 
TransIntOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, InternalOppositeGenes, Orthos) 
CisExtOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, ExternalSameGenes, Orthos) 
TransExtOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, ExternalOppositeGenes, Orthos) 
NoOvlOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, NonOverlappingGenes, Orthos) 
 

InternalOrthos = [pair for pair in CisIntOrthos]
InternalOrthos.extend(TransIntOrthos)
ExternalOrthos = [pair for pair in CisExtOrthos] 
ExternalOrthos.extend(TransExtOrthos)
 
 
 

a = [NstOrthos, InternalOrthos, ExternalOrthos, OvlOrthos, ConOrthos, DivOrthos, PbkOrthos, CisIntOrthos, TransIntOrthos, CisExtOrthos, TransExtOrthos, NoOvlOrthos]


b = ['NstOrthos', 'InternalOrthos', 'ExternalOrthos', 'OvlOrthos', 'ConOrthos', 'DivOrthos', 'PbkOrthos', 'CisIntOrthos', 'TransIntOrthos', 'CisExtOrthos', 'TransExtOrthos', 'NoOvlOrthos']


# compute expression divergence between orthologs for each gene category
ExpressionDivergence = []
for i in range(len(a)):
    Div = ComputeExpressionDivergenceOrthologs(a[i], HumanExpression, ChimpExpression)
    ExpressionDivergence.append(Div)
    

for i in range(len(ExpressionDivergence)):
    print(b[i], len(a[i]), np.mean(ExpressionDivergence[i]))    


for i in range(0, len(ExpressionDivergence)-1):
    P = stats.ranksums(ExpressionDivergence[i], ExpressionDivergence[-1])[1]
    print(b[i], b[-1], P)







    
#NstOrthos, 463, 0.222439523431
#OvlOrthos, 3464, 0.195316554697
#ConOrthos, 1591, 0.193772155802
#DivOrthos, 1387, 0.185638705044
#PbkOrthos, 554, 0.203905539341
#CisIntOrthos, 26, 0.265285371342
#TransIntOrthos, 181, 0.248025865567
#CisExtOrthos, 52, 0.178863021383
#TransExtOrthos, 216, 0.205683446829
#NoOvlOrthos, 8724, 0.21121773543
































## generate lists of gene pairs separated by distance 
#Proximal, Moderate, Intermediate, Distant = GenerateSetsGenePairsDistance(GeneCoord, OrderedGenes, ExpressionProfile)
#
## filter gene pairs lacking expression
#AllPairs = [NestedPairs, PiggybackPairs, ConvergentPairs, DivergentPairs]
#for i in range(len(AllPairs)):
#    AllPairs[i] = FilterGenePairsWithoutExpression(AllPairs[i], ExpressionProfile)
#
## add gene pairs to Allpairs list
#AllPairs.append(Proximal)
#AllPairs.append(Moderate)
#AllPairs.append(Intermediate)
#AllPairs.append(Distant)
#
## compute expression divergence between pairs of genes
#Divergence = []
#for i in range(len(AllPairs)):
#    Div = ComputeExpressionDivergenceGenePairs(AllPairs[i], ExpressionProfile)
#    Divergence.append(Div)
#
## make a list of gene category names parallel to the list of gene pairs
#GeneCats = ['Nst', 'Pbk', 'Conv', 'Div', 'Prox', 'Mod', 'Int', 'Dist']
#
## create lists with means and SEM for each gene category
#MeanExpDiv, SEMExpDiv = [], []
## loop over lists in Divergence list
#for i in range(len(Divergence)):
#    MeanExpDiv.append(np.mean(Divergence[i]))
#    SEMExpDiv.append(np.std(Divergence[i]) / math.sqrt(len(Divergence[i])))
#
## create figure
#fig = plt.figure(1, figsize = (3, 2))
#
## add a plot to figure (N row, N column, plot N)
#ax = fig.add_subplot(1, 1, 1)
## set colors
#colorscheme = ['grey','grey','grey','grey', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
## plot nucleotide divergence
#ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
#       edgecolor = 'black', linewidth = 0.5,
#       error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
## set font for all text in figure
#FigFont = {'fontname':'Arial'}   
## write y axis label
#ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
## add ticks and lebels
#plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
## add a range for the Y and X axes
#plt.ylim([0, 0.61])
#plt.xlim([0, 2.45])
## do not show lines around figure  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(True)  
## edit tick parameters    
#plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                right = 'off', left = 'on', labelbottom='on',
#                colors = 'black', labelsize = 7, direction = 'out')  
## Set the tick labels font name
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')   
#      
#
## perform statistical tests between gene categories
## create list to store the p-values
#PValues = []
## loop over inner list, compare gene categories
#for i in range(0, len(Divergence) -1):
#    for j in range(i+1, len(Divergence)):
#        P = stats.ranksums(Divergence[i], Divergence[j])[1]
#        PValues.append(P)
## print p values
#for p in PValues:
#    print(p)
#
## convert p-values to star significance level
#Significance = []
#for pvalue in PValues:
#    if pvalue >= 0.05:
#        Significance.append('')
#    elif pvalue < 0.05 and pvalue >= 0.01:
#        Significance.append('*')
#    elif pvalue < 0.01 and pvalue >= 0.001:
#        Significance.append('**')
#    elif pvalue < 0.001:
#        Significance.append('***')
#
#
## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
#ypos = [0.55] * 8
#xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
#for i in range(len(Diff)):
#    ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#    
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')
