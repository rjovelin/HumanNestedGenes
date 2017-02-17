# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 11:56:27 2017

@author: RJovelin
"""

# use this script to compare the proportions of human overlapping genes conserved in chimp and mouse
# and the proportions of overlapping genes in chimp and mouse conserved in human


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

# load dictionaries of overlapping genes
jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json', 'HumanPiggyBackGenes.json',
             'HumanConvergentGenes.json', 'HumanDivergentGenes.json',
             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json', 'ChimpPiggyBackGenes.json',
             'ChimpConvergentGenes.json', 'ChimpDivergentGenes.json',
             'MouseOverlappingGenes.json', 'MouseNestedGenes.json', 'MousePiggyBackGenes.json',
             'MouseConvergentGenes.json', 'MouseDivergentGenes.json']

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
GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3', 'Mus_musculus.GRCm38.86.gff3']

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
HumanOrdered, ChimpOrdered, MouseOrdered = AllOrdered[0], AllOrdered[1], AllOrdered[2]
HumanCoord, ChimpCoord, MouseCoord = AllCoordinates[0], AllCoordinates[1], AllCoordinates[2]

# get 1:1 orthologs between human and other species
HsaPtrOrthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
HsaMmuOrthos = MatchOrthologPairs('HumanMouseOrthologs.txt')

# make a reverse dictionary of orthologs 
ChimpOrthos = {}
for gene in HsaPtrOrthos:
    assert HsaPtrOrthos[gene] not in ChimpOrthos
    ChimpOrthos[HsaPtrOrthos[gene]] = gene
MouseOrthos = {}
for gene in HsaMmuOrthos:
    assert HsaMmuOrthos[gene] not in MouseOrthos
    MouseOrthos[HsaMmuOrthos[gene]] = gene

# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# make lissts of gene pairs
HumanPairs = AllPairs[:5]
HsaPairs = copy.deepcopy(HumanPairs)
ChimpPairs = AllPairs[5:10]
MousePairs = AllPairs[10:]

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in HsaPtrOrthos or pair[1] not in HsaPtrOrthos:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
for i in range(len(HsaPairs)):
    to_remove = []
    for pair in HsaPairs[i]:
        if pair[0] not in HsaMmuOrthos or pair[1] not in HsaMmuOrthos:
            to_remove.append(pair)
    for pair in to_remove:
        HsaPairs[i].remove(pair)
# remove chimp genes lacking orthologs
for i in range(len(ChimpPairs)):
    to_remove = []
    for pair in ChimpPairs[i]:
        if pair[0] not in ChimpOrthos or pair[1] not in ChimpOrthos:
            to_remove.append(pair)
    for pair in to_remove:
        ChimpPairs[i].remove(pair)
# remove mouse genes lacking orthologs
for i in range(len(MousePairs)):
    to_remove = []
    for pair in MousePairs[i]:
        if pair[0] not in MouseOrthos or pair[1] not in MouseOrthos:
            to_remove.append(pair)
    for pair in to_remove:
        MousePairs[i].remove(pair)

# make a list of sets of gene pairs
HumanSets, HsaSets, ChimpSets, MouseSets = [], [], [], []
for i in range(len(HumanPairs)):
    HumanSets.append([set(j) for j in HumanPairs[i]])
for i in range(len(HsaPairs)):
    HsaSets.append([set(j) for j in HsaPairs[i]])    
for i in range(len(ChimpPairs)):
    ChimpSets.append([set(j) for j in ChimpPairs[i]])
for i in range(len(MousePairs)):
    MouseSets.append([set(j) for j in MousePairs[i]])

# do qc
for i in range(1, len(ChimpSets)):
    for pair in ChimpSets[i]:
        assert pair in ChimpSets[0]
for i in range(1, len(MouseSets)):
    for pair in MouseSets[i]:
        assert pair in MouseSets[0]

# make a list of counts of conserved and non-conserved human overlapping genes in chimp
HumanConserved = []
for i in range(len(HumanPairs)):
    conserved, divergent = 0, 0
    for pair in HumanPairs[i]:
        if set([HsaPtrOrthos[pair[0]], HsaPtrOrthos[pair[1]]]) in ChimpSets[i]:
            conserved += 1
        else:
            divergent += 1
    HumanConserved.append([conserved, divergent])
# make a list of counts of conserved and non-conserved human overlapping genes in mouse
HsaConserved = []
for i in range(len(HsaPairs)):
    conserved, divergent = 0, 0
    for pair in HsaPairs[i]:
        if set([HsaMmuOrthos[pair[0]], HsaMmuOrthos[pair[1]]]) in MouseSets[i]:
            conserved += 1
        else:
            divergent += 1
    HsaConserved.append([conserved, divergent])
# make a list of counts of conserved and non-conserved chimp overlapping genes in human
ChimpConserved = []
for i in range(len(ChimpPairs)):
    conserved, divergent = 0, 0
    for pair in ChimpPairs[i]:
        if set([ChimpOrthos[pair[0]], ChimpOrthos[pair[1]]]) in HumanSets[i]:
            conserved += 1
        else:
            divergent += 1
    ChimpConserved.append([conserved, divergent])
# make a list of counts of conserved and non-conserved mouse overlapping genes in human
MouseConserved = []
for i in range(len(MousePairs)):
    conserved, divergent = 0, 0
    for pair in MousePairs[i]:
        if set([MouseOrthos[pair[0]], MouseOrthos[pair[1]]]) in HsaSets[i]:
            conserved += 1
        else:
            divergent += 1
    MouseConserved.append([conserved, divergent])

# create lists with proportions of overlapping genes with conserved topologies
# [prop human overlapping genes conserved, prop chimp overlapping genes conserved in human, etc]
HumanProp, HsaProp = [], []
for i in range(len(HumanConserved)):
    man = HumanConserved[i][0] / sum(HumanConserved[i])    
    sp2 = ChimpConserved[i][0] / sum(ChimpConserved[i])
    diff = (sp2 - man) * 100 
    HumanProp.append(diff) 
for i in range(len(HsaConserved)):
    man = HsaConserved[i][0] / sum(HsaConserved[i])
    sp2 = MouseConserved[i][0] / sum(MouseConserved[i])
    diff = (sp2 - man) * 100
    HsaProp.append(diff)
  

Differences = []
for i in range(len(HumanProp)):
    Differences.append(HumanProp[i])
    Differences.append(HsaProp[i])

# create figure
fig = plt.figure(1, figsize = (3, 2))

# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# plot data
Colors = ['black', 'lightgrey'] * 5    
BarPos = [0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4]
ax.bar(BarPos, Differences, width = 0.2, color = Colors, edgecolor = 'black', linewidth = 0.7)
ax.plot([0, 2.8], [0, 0], lw = 0.7, color = 'black')


# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('% excess of orthologous gene pairs\nwith conserved topology', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
XTickpos = [0.4, 0.9, 1.4, 1.9, 2.4]    
XTicklabels = ['all', 'nst', 'pbk', 'conv', 'div']
plt.xticks(XTickpos, XTicklabels, rotation = 0, size = 7, color = 'black', ha = 'center', **FigFont)
# edit y axis ticks
plt.yticks(np.arange(-20, 120, 20))   
# add a range for the Y and X axes
plt.ylim([-20, 100])    
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)  

# make sure the y axis crosses the x axis at 0
ax.spines['left'].set_position('zero')

# edit tick parameters    
plt.tick_params(axis='both', which='both', bottom='off', top='off',
                right = 'off', left = 'on', labelbottom='on',
                colors = 'black', labelsize = 7, direction = 'out')  
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   

# add legend
mouse = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'mouse')
chimp = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'chimp')
ax.legend(handles = [chimp, mouse], loc = (0.2, 0.8), fontsize = 7, frameon = False, ncol = 2)

# save figure to file
fig.savefig('DifferencesTopologyConservation.pdf', bbox_inches = 'tight')
fig.savefig('DifferencesTopologyConservation.eps', bbox_inches = 'tight')