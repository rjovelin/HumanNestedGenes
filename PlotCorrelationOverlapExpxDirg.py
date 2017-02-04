# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 11:13:11 2017

@author: RJovelin
"""


# use this script to plot the correlation between overlap length and expression divergence

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

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)
NestedPairs = GetHostNestedPairs(Nested)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)

# create a list of lists of gene pairs
AllPairs = [OverlappingPairs, NestedPairs, PiggybackPairs, ConvergentPairs, DivergentPairs]
# create a list of dicts with the overlap length for each pair {'gene1:gene2' : overlaplength}
AllLength = []
# loop over lists of gene pairs for each overlapping group
for i in range(len(AllPairs)):
    # initialize dict
    Length = {}
    for pair in AllPairs[i]:
        assert GeneCoord[pair[0]][0] == GeneCoord[pair[1]][0]
        # get the gene pair
        genepair = pair[0] + ':' + pair[1]        
        # get the coordinates of the gene pairs
        coord1 = set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2]))
        coord2 = set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2]))
        L = len(coord1.intersection(coord2))
        assert genepair not in Length        
        Length[genepair] = L
    # store dict of overlap length
    AllLength.append(Length)
# get the dicts of overlap length for each group
OverlapLength, NestedLength = AllLength[0], AllLength[1]
PiggybackLength, ConvergentLength, DivergentLength = AllLength[2], AllLength[3], AllLength[4]



# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseGTEXExpression('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# filter gene pairs lacking expression
for i in range(len(AllPairs)):
    AllPairs[i] = FilterGenePairsWithoutExpression(AllPairs[i], ExpressionProfile)

# compute expression divergence between pairs of genes
Divergence = []
for i in range(len(AllPairs)):
    Div = {} 
    # loop over the list of gene pairs
    for j in range(len(AllPairs[i])):
        # compute divergence between the 2 genes in the given pair
        gene1, gene2 = AllPairs[i][j][0], AllPairs[i][j][1]        
        D = EuclidianDistance(ExpressionProfile[gene1], ExpressionProfile[gene2])
        genepair = gene1 + ':' + gene2
        assert genepair not in Div
        Div[genepair] = D
    Divergence.append(Div)

OverlapDivg, NestedDivg = Divergence[0], Divergence[1]
PiggybackDivg, ConvergentDivg, DivergentDivg = Divergence[2], Divergence[3], Divergence[4]

   
# use this function to match the overlap length and the expression divergence for the same pair
def GetExpxDigOverlap(Overlap, ExpDiv):
    '''
    (dict, dict) -> list
    Take the dictionaries of overlap length and expression divergence
    and return a list of lists with matched overlap and length and expression divergence
    '''
    # check that pairs are recorded once
    a = list(Overlap.keys())
    b = list(ExpDiv.keys())
    for i in range(len(a)):
        a[i] = set(a[i].split(':'))
    for i in range(len(b)):
        b[i] = set(b[i].split(':'))
    c, d = [], []
    for i in a:
        c.append(b.count(i))
    for i in b:
        d.append(a.count(i))
    c, d = set(c), set(d)
    assert sum(c) <= 1 and sum(d) <= 1
    
    # create a list with overlap and expression divergence
    OverlpLengthDiv = []
    for i in OverlapDivg:
        k = set(i.split(':'))
        # find the overlap length for the corresponding pair
        for j in OverlapLength:
            m = set(j.split(':'))
            if k == m:
                OverlpLengthDiv.append([OverlapDivg[i], OverlapLength[j]])
    return OverlpLengthDiv


# create lists with matched overlap length and expression divergence
OverlapLengthDiv = GetExpxDigOverlap(OverlapLength, OverlapDivg)
NestedLengthDiv = GetExpxDigOverlap(NestedLength, NestedDivg)
PiggybackLengthDiv = GetExpxDigOverlap(PiggybackLength, PiggybackDivg)
ConvergentLengthDiv = GetExpxDigOverlap(ConvergentLength, ConvergentDivg)
DivergentLengthDiv = GetExpxDigOverlap(DivergentLength, DivergentDivg)



# make a list of gene category names parallel to the list of gene pairs
#GeneCats = ['Nst', 'Pbk', 'Conv', 'Div', 'Prox', 'Mod', 'Int', 'Dist']
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
#
#
#
#
#
#
#
#
#
#
#
###################################
#
#
## convert bp to Kbp
#ToKb = lambda x: x / 1000
#OverlapLength = list(map(ToKb, OverlapLength))
#NestedLength = list(map(ToKb, NestedLength))
#PiggybackLength = list(map(ToKb, PiggybackLength))
#ConvergentLength = list(map(ToKb, ConvergentLength))
#DivergentLength = list(map(ToKb, DivergentLength))
