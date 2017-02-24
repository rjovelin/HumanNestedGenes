# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:46:17 2017

@author: RJovelin
"""

# use this script to plot the % of orthologous gene pairs with same topology between human and chimp

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

print('ordered genes')
print('got gene coordinates')

# get 1:1 orthologs between human and chimp
HsaPtrOrthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
HsaMmuOrthos = MatchOrthologPairs('HumanMouseOrthologs.txt')

# make a list of dictionaries with orthologs
Orthos = [HsaPtrOrthos, HsaMmuOrthos]

print('mapped orthologs')


# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# get the gene pairs
HumanPairs = AllPairs[:5]
# copy human gene pairs
HsaPairs = copy.deepcopy(HumanPairs)
ChimpPairs = AllPairs[5:10]
MousePairs = AllPairs[10:]

print('made lists of gene pairs')

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos[0] or pair[1] not in Orthos[0]:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
for i in range(len(HsaPairs)):
    to_remove = []
    for pair in HsaPairs[i]:
        if pair[0] not in Orthos[1] or pair[1] not in Orthos[1]:
            to_remove.append(pair)
    for pair in to_remove:
        HsaPairs[i].remove(pair)        
print('removed human gene pairs lacking orthologs')

# create a dict with lists of gene pairs
AllGenePairs = {}
AllGenePairs['hsa_ptr'] = copy.deepcopy(HumanPairs)
AllGenePairs['hsa_mmu'] = copy.deepcopy(HsaPairs)
AllGenePairs['ptr'] = copy.deepcopy(ChimpPairs)
AllGenepairs['mmu'] = copy.deepcopy(MousePairs)

# replace human genes by their orthologs
species = ['hsa_ptr', 'hsa_mmu']
for i in range(len(species)):
    for j in range(len(AllGenePairs[species[i]])):
        for k in range(len(AllGenePairs[species[i]][j])):
            AllGenePairs[species[i]][j][k][0] = Orthos[i][AllGenePairs[species[i]][j][k][0]]
            AllGenePairs[species[i]][j][k][1] = Orthos[i][AllGenePairs[species[i]][j][k][1]]

# sort gene pairs
for sp in AllGenePairs:
    for i in range(len(AllGenePairs[sp])):
        for j in range(len(AllGenePairs[sp][i])):
            AllGenePairs[sp][i][j].sort()

# replace gene lists with strings
for sp in AllGenePairs:
    for i in range(len(AllGenePairs[sp])):
        for j in range(len(AllGenePairs[sp][i])):
            AllGenePairs[sp][i][j] = ':'.join(AllGenePairs[sp][i][j])

# do qc
for sp in AllGenePairs:
    for i in range(1, len(AllGenePairs[sp])):
        for pair in AllGenePairs[sp][i]:
            assert pair in AllGenePairs[sp][0]

print('done with QC')


# make lists of adjacent gene pairs with varying distance in human conserved in chimp and mouse
# [[gene1, gene2], ....[gene n, gene n+1]]
for i in range(len(species)):
    HsaPairsDist = []    
    for j in range(7):
        HsaPairsDist.append([])
    # loop over chromosomes
    for chromo in HumanOrdered:
        # loop over the list of ordered genes
        for k in range(len(HumanOrdered[chromo]) - 1):
            # get the end position of gene 1
            EndGene1 = HumanCoord[HumanOrdered[chromo][k]][2]                
            # get the start position of adjacent gene 2
            StartGene2 = HumanCoord[HumanOrdered[chromo][k+1]][1]
            # computance distance between genes
            D = StartGene2 - EndGene1
            # assign infinity value to k
            m = float('inf')
            if D < 0:
                m = 0
            elif D >= 0 and D < 1000:
                m = 1
            elif D >= 1000 and D < 10000:
                m = 2                
            elif D >= 10000 and D < 50000:
                m = 3
            elif D >= 50000 and D < 100000:
                m = 4
            elif D >= 100000 and D < 150000:
                m = 5
            elif D >= 150000:
                m = 6
            if HumanOrdered[chromo][k] in Orthos[i] and HumanOrdered[chromo][k+1] in Orthos[i] and m in range(7):
                # add the human orthologs
                pair = [Orthos[i][HumanOrdered[chromo][k]], Orthos[i][HumanOrdered[chromo][k+1]]]
                # sort gene pair
                pair.sort()
                # concert list to string
                pair = ':'.join(pair)
                HsaPairsDist[m].append(pair)
    # add the pairs of non-overlapping genes to the lists of gene pairs
    AllGenePaies[species[i]].extend(HsaPairsDist)

print('generated human gene pairs by distance')


# make list of adjacent gene pairs regardless of distance in chimp and mouse [gene1:gene2, ....gene_n:gene_n+1]
sp2 = ['ptr', 'mmu']
Sp2Coord = [ChimpCoord, MouseCoord]
Sp2Ordered = [ChimpOrdered, MouseOrdered]
for i in range(len(sp2)):
    Sp2PairsDist = []
    # loop over chromosomes
    for chromo in Sp2Ordered[i]:
        # loop over the list of ordered genes
        for j in range(len(Sp2Ordered[i][chromo]) - 1):
            # get the end position of gene 1
            EndGene1 = Sp2Coord[i][Sp2Ordered[i][chromo][j]][2]                
            # get the start position of adjacent gene 2
            StartGene2 = Sp2Coord[i][Sp2Ordered[i][chromo][j+1]][1]
            # check that genes are not overlapping
            D = StartGene2 - EndGene1
            if D >= 0:
                # get gene pair
                pair = [Sp2Ordered[i][chromo][j], Sp2Ordered[i][chromo][j+1]]
                # sort pair
                pair.sort()
                # convert list to string
                pair = ':'.join(pair)
                Sp2PairsDist.append(pair)
    AllGenePairs[sp2[i]].append(Sp2PairsDist)

print('generated species 2 gene pairs by distance')


# convert lists to numpy arrays
for sp in AllGenePairs:
    for i in range(len(AllGenePairs[sp])):
        AllGenePairs[sp][i] = np.array(AllGenePairs[sp][i])

for sp in AllGenePairs:
    for i in range(len(AllGenePairs[sp])):
        print(sp, i, len(AllGenePairs[sp][i]))

################## continue here
##### a, distribution of conservation by distance in chimp and mouse
##### b, distribution of overlapping conservation in each category and overlapping in chimp
##### c, distribution of overlapping conservation in each category and overlapping in mouse
##### d differences between conservation opf human overalpping in other species and overlapping genes in opther species in human


# count the number of pairs with conserved topology
CountPairs = []
for i in range(len(HumanPairs)):
    if i < 5:
        total = sum(np.in1d(HumanPairs[i], Sp2Pairs[i], invert = False))
    else:
        total = sum(np.in1d(HumanPairs[i], Sp2Pairs[5]))
    CountPairs.append([total, len(HumanPairs[i])])



#Sp2Genes = []
#for i in range(len(Sp2Pairs)):
#    Sp2Genes.extend(Sp2Pairs[i])
#CountPairs = []
#for i in range(len(HumanPairs)):
#    total = sum(np.in1d(HumanPairs[i], Sp2Genes, invert = False))
#    CountPairs.append([total, len(HumanPairs[i])])



## count the number of pairs with conserved topology
#CountPairs = []
#for i in range(len(HumanPairs)):
#    if i < 5:
#        total = sum(np.in1d(HumanPairs[i], Sp2Pairs[0], invert = False))
#        CountPairs.append([total, len(HumanPairs[i])])


# create a list of overlapping gene categories parallel to the list of overlapping pairs
GeneCats = ['overlapping', 'nested', 'piggyback', 'convergent', 'divergent',
            'proximal', 'moderate', 'intermediate', 'distant']

for i in range(len(CountPairs)):
    print(str(i) + '\t' + str(CountPairs[i][0] / CountPairs[i][1]))


#newfile = open('PairsCounts.txt', 'a')
#header = '\t'.join(['species', 'genes', 'conserved', 'total', 'ratio'])
#newfile.write(header + '\n')
#
#for i in range(len(GeneCats)):
#    line = '\t'.join([species, GeneCats[i], str(PairCounts[GeneCats[i]][0]), str(PairCounts[GeneCats[i]][1]), str(PairCounts[GeneCats[i]][0] / PairCounts[GeneCats[i]][1])])
#    newfile.write(line + '\n')
#newfile.close()


#newfile = open('test.txt', 'w')
#newfile.write('a\n')
#newfile.write(str(sum(np.in1d(a, a))) + '\n')
#newfile.write(str(sum(np.in1d(a, a, invert = False)) / len(a)) + '\n')
#newfile.write('\n\n')
#
#
#truc = '\t'.join([str(sum(np.in1d(b, b))), str(sum(np.in1d(b, b[:15]))), str(sum(np.in1d(b, b[:20]))), str(sum(np.in1d(b, b[:30]))), str(sum(np.in1d(b, b[:50])))])
#newfile.write(truc + '\n')
#newfile.close()




##############################


# make a reverse dictionary of orthologs 
ChimpOrthos = {}
for gene in HsaPtrOrthos:
    assert HsaPtrOrthos[gene] not in ChimpOrthos
    ChimpOrthos[HsaPtrOrthos[gene]] = gene
MouseOrthos = {}
for gene in HsaMmuOrthos:
    assert HsaMmuOrthos[gene] not in MouseOrthos
    MouseOrthos[HsaMmuOrthos[gene]] = gene


# remove gene pairs lacking orthologs
for i in range(len(ChimpPairs)):
    to_remove = []
    for pair in ChimpPairs[i]:
        if pair[0] not in ChimpOrthos or pair[1] not in ChimpOrthos:
            to_remove.append(pair)
    for pair in to_remove:
        ChimpPairs[i].remove(pair)
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
fig.savefig('truc.pdf', bbox_inches = 'tight')
