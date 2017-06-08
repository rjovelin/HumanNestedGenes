# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 18:25:55 2017

@author: RJovelin
"""

# use this script to plot conservation of overlapping genes

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


# 1) plot the proportion of gene in each overlapping gene category that has orthologs in other species

# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
           'Mouse', 'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
           'Cat', 'Opossum', 'Platypus', 'Shrew']

# make a list of json files
JsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json', 'HumanPiggyBackGenes.json', 'HumanConvergentGenes.json', 'HumanDivergentGenes.json']

# make a parallel list of ortholog files
OrthoFiles = ['Human' + i + 'Orthologs.txt' for i in Species[1:]]

# make lists of dictionaries for each type of overlapping gene
Overlap  = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary from json file
    json_data = open(JsonFiles[i])
    Overlap.append(json.load(json_data))
    json_data.close()

# make a list of gene sets
AllGenes = []
for i in range(len(Overlap)):
    AllGenes.append(MakeFullPartialOverlapGeneSet(Overlap[i]))

# make a set of non-overlapping genes
# get GFF file
GFF = 'Homo_sapiens.GRCh38.88.gff3'
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

# make sets of external and internal genes
External, Internal = set(), set()
# make a list of nested gene pairs
NestedGenePairs = GetHostNestedPairs(Overlap[1])
for pair in NestedGenePairs:
     External.add(pair[0])
     Internal.add(pair[1])
    
# make a list of dicts with orthologs
AllOrthologs = []
for i in range(len(OrthoFiles)):
    orthologs = MatchOrthologs(OrthoFiles[i])
    AllOrthologs.append(orthologs)
    
# compute proportions of overlapping genes with orthologs in each species
# make a list of sets of genes of interest [nonovl, nst, ext, int, pbk, con, div]
GenesInterest = [NonOverlappingGenes, AllGenes[1], External, Internal, AllGenes[2], AllGenes[3], AllGenes[4]]
Conserved = []
for i in range(len(GenesInterest)):
    # make a list to store proportions of given overlapping gene conserved across each species
    ConservedAccrossSpecies = [0] * len(AllOrthologs)
    for gene in GenesInterest[i]:
        # loop over dicts of orthologous genes
        for j in range(len(AllOrthologs)):
            # check if gene as orthologs
            if gene in AllOrthologs[j]:
                # update counter for given species
                ConservedAccrossSpecies[j] += 1
    # divide by number of genes to get proportion of human genes
    for j in range(len(ConservedAccrossSpecies)):
        ConservedAccrossSpecies[j] = ConservedAccrossSpecies[j] / len(GenesInterest[i])
    # add list for given gene category
    Conserved.append(ConservedAccrossSpecies)

# convert list to numpy array
Conserved = np.array(Conserved)
# transpose array to get gene categories as columns and species as rows
Conserved = np.transpose(Conserved)


# 2) plot conservation of gene pairs, full, partial or none between human and mouse

# make parallel lists of json files with each overlapping classes [human, mouse]
NestedFiles = [i + 'NestedGenes.json' for i in ['Human', 'Mouse']]
PiggybackFiles = [i + 'PiggyBackGenes.json' for i in ['Human', 'Mouse']]
ConvergentFiles = [i + 'ConvergentGenes.json' for i in ['Human', 'Mouse']]
DivergentFiles = [i + 'DivergentGenes.json' for i in ['Human', 'Mouse']]

# make lists of dictionaries for each type of overlapping gene
AllOverlapGenes  = []
# loop over files
AllFiles = [NestedFiles, PiggybackFiles, ConvergentFiles, DivergentFiles]
for i in range(len(AllFiles)):
    OvDicts = []
    for j in range(len(AllFiles[i])):
        # load dictionary from json file
        json_data = open(AllFiles[i][j])
        OvDicts.append(json.load(json_data))
        json_data.close()
    AllOverlapGenes.append(OvDicts)
for i in range(len(AllOverlapGenes)):
    assert len(AllOverlapGenes[i]) == 2
    
# for each overlapping gene class, count the proportion of gene pairs with both genes
# in the same configuration, 1 gene only or none
HumanMouseConservation = []

# get nested pairs 
NstPairs = [GetHostNestedPairs(AllOverlapGenes[0][i]) for i in range(len(AllOverlapGenes[0]))]
PbkPairs = [GetHostNestedPairs(AllOverlapGenes[1][i]) for i in range(len(AllOverlapGenes[1]))]
ConPairs = [GetHostNestedPairs(AllOverlapGenes[2][i]) for i in range(len(AllOverlapGenes[2]))]
DivPairs = [GetHostNestedPairs(AllOverlapGenes[3][i]) for i in range(len(AllOverlapGenes[3]))]

# get orthologs between human and mouse
Orthologs = MatchOrthologs('HumanMouseOrthologs.txt')
# reverse dictionary 
MouseOrthologs = {}
for gene in Orthologs:
    for ortho in Orthologs[gene]:
        if ortho not in MouseOrthologs:
            MouseOrthologs[ortho] = [gene]
        else:
            MouseOrthologs[ortho].append(gene)

# make pairs of human orthologs for each mouse gene pair
PairsOrthos = []
for L in [NstPairs, PbkPairs, ConPairs, DivPairs]:
    # create a list to store the gene pairs
    GenePairs = []
    # loop over mouse gene pairs
    for pair in L[1]:
        # count only pairs in which both genes have orthos
        if pair[0] in MouseOrthologs and pair[1] in MouseOrthologs:
            for ortho1 in MouseOrthologs[pair[0]]:
                for ortho2 in MouseOrthologs[pair[1]]:
                    # remove order
                    GenePairs.append(set([ortho1, ortho2]))
    PairsOrthos.append(GenePairs)

# make a list with pair counts for each overlapping gene type
PairCounts = []
AllPairs = [NstPairs, PbkPairs, ConPairs, DivPairs]

# loop over human overlapping gene classes
for i in range(len(AllPairs)):
    # create a list with pair counts [both gene conserved, 1 conserved, non conserved]
    counts = [0, 0, 0]
    # make a set of orthologous overlapping gene
    orthosovlp = set()
    for pair in PairsOrthos[i]:
        pair = list(pair)
        for item in pair:
            orthosovlp.add(item)
    # loop over human gene pairs
    for pair in AllPairs[i][0]:
        if set(pair) in PairsOrthos[i]:
            counts[0] += 1
        else:
            if pair[0] in orthosovlp or pair[1] in orthosovlp:
                counts[1] += 1
            elif pair[0] not in orthosovlp and pair[1] not in orthosovlp:
                counts[2] += 1
    # divide by number of gene pairs to get proportions    
    assert sum(counts) == len(AllPairs[i][0])
    for j in range(len(counts)):
        counts[j] = counts[j] / len(AllPairs[i][0])
    PairCounts.append(counts)
    assert sum(counts) == 1


# 3) plot the % of orthologous gene pairs with same topology between human and mouse

# load dictionaries of overlapping genes
JsonFiles = ['Overlapping', 'Nested', 'PiggyBack', 'Convergent', 'Divergent']
# make a list of dictionaries
HsaAllOverlap, MmuAllOverlap = [], []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open('Human' + JsonFiles[i] + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    HsaAllOverlap.append(overlapping)
    json_data = open('Mouse' + JsonFiles[i] + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    MmuAllOverlap.append(overlapping)

# get GFF file
GFF = ['Homo_sapiens.GRCh38.88.gff3', 'Mus_musculus.GRCm38.88.gff3']

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
HumanOrdered, MouseOrdered = AllOrdered[0], AllOrdered[1]
HumanCoord, MouseCoord = AllCoordinates[0], AllCoordinates[1]

# get orthologs between human and mouse
Orthos = MatchOrthologs('HumanMouseOrthologs.txt')
# reverse dictionary 
MouseOrthologs = {}
for gene in Orthos:
    for ortho in Orthos[gene]:
        if ortho not in MouseOrthologs:
            MouseOrthologs[ortho] = [gene]
        else:
            MouseOrthologs[ortho].append(gene)

# record adjacent gene pairs in human (remove order)
HumanAdjacentgenePairs = {}
for i in ['Proximal', 'Moderate', 'Intermediate', 'Distant']:
    HumanAdjacentgenePairs[i] = []
# loop over chromo in human
for chromo in HumanOrdered:
    # loop over the list of ordered genes
    for i in range(len(HumanOrdered[chromo]) -1):
        gene1, gene2 = HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]
        # only consider gene pairs with orthologs
        if gene1 in Orthos and gene2 in Orthos:
            # get the end position of gene 1
            EndGene1 = HumanCoord[gene1][2]
            # get the start position of adjacent gene 2
            StartGene2 = HumanCoord[gene2][1]
            # computance distance between genes
            D = StartGene2 - EndGene1
            if D >= 0 and D < 1000:
                HumanAdjacentgenePairs['Proximal'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))
            elif D >= 1000 and D < 10000:
                HumanAdjacentgenePairs['Moderate'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))
            elif D >= 10000 and D < 50000:
                HumanAdjacentgenePairs['Intermediate'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))            
            elif D >= 50000:
                HumanAdjacentgenePairs['Distant'].append(set([HumanOrdered[chromo][i], HumanOrdered[chromo][i+1]]))

# make pairs of overlapping genes
for i in range(1, len(HsaAllOverlap)):
    pairs = GetHostNestedPairs(HsaAllOverlap[i])
    # remove pairs if any gene is lacking an ortholog
    to_remove = [L for L in pairs if L[0] not in Orthos or L[1] not in Orthos]    
    for L in to_remove:
        pairs.remove(L)
    HumanAdjacentgenePairs[JsonFiles[i]] = pairs    
    
# count gene pairs conserved in mouse    
ConservedPairs = {}
for GeneType in HumanAdjacentgenePairs:
    # initialize counters
    ConservedPairs[GeneType] = 0
    if GeneType in ['Proximal', 'Moderate', 'Intermediate', 'Distant']:
        # check if mouse orthologs are adjacent
        for pair in HumanAdjacentgenePairs[GeneType]:
            pair = list(pair)
            # make a list of orthologs in mouse
            L = []
            for ortho1 in Orthos[pair[0]]:
                for ortho2 in Orthos[pair[1]]:
                    L.append([ortho1, ortho2])
            for genes in L:
                gene1, gene2 = genes[0], genes[1]
                # check if genes are valid mouse genes
                if gene1 in MouseCoord and gene2 in MouseCoord:
                    # check if genes are on the same chromo and that genes are different
                    if (MouseCoord[gene1][0] == MouseCoord[gene2][0]) and (gene1 != gene2):
                        # get indices of gene1 and gene2
                        I1, I2 = MouseOrdered[MouseCoord[gene1][0]].index(gene1), MouseOrdered[MouseCoord[gene2][0]].index(gene2)
                        assert I1 != I2                        
                        # check if genes are adjacent
                        if I1 == I2 + 1 or I2 == I1 + 1:
                            ConservedPairs[GeneType] += 1
                        # record only 1 pair of orthologs if multiple co-ortholiogs exit, exit loop 
                        break
                
# make pairs of human orthologs for each mouse gene pair
PairsOrthos = []
for i in range(1, len(MmuAllOverlap)):
    # create a list to store the gene pairs
    GenePairs = []
    Pairs = GetHostNestedPairs(MmuAllOverlap[i])
    # loop over mouse gene pairs
    for pair in Pairs:
        # check that both genes have coordinates
        assert pair[0] in MouseCoord and pair[1] in MouseCoord
        # count only pairs in which both genes have orthologs in human
        if pair[0] in MouseOrthologs and pair[1] in MouseOrthologs:
            for ortho1 in MouseOrthologs[pair[0]]:
                for ortho2 in MouseOrthologs[pair[1]]:
                    # remove order
                    GenePairs.append(set([ortho1, ortho2]))
    PairsOrthos.append(GenePairs)
# loop over human overlapping gene classes    
for i in range(1, len(JsonFiles)):
    for pair in HumanAdjacentgenePairs[JsonFiles[i]]:
        if set(pair) in PairsOrthos[i-1]:
            ConservedPairs[JsonFiles[i]] += 1
# compute proportions
for i in ConservedPairs:
    ConservedPairs[i] = ConservedPairs[i] / len(HumanAdjacentgenePairs[i])    
    
# test differences among gene categories
# write P values to file
newfile = open('TopologyConservationPvalues.txt', 'w')
GeneTypes = JsonFiles[1:] + ['Proximal', 'Moderate', 'Intermediate', 'Distant']
Header = '\t'.join(['\t']+GeneTypes)
newfile.write(Header + '\n')
for i in range(len(GeneTypes)):
    # create a list with pvalues for each line
    Line = [GeneTypes[i]]
    for j in range(len(GeneTypes)):
        # add p values for each pairwise comparison
        Line.append(str(stats.fisher_exact([[ConservedPairs[GeneTypes[i]], len(HumanAdjacentgenePairs[GeneTypes[i]])-ConservedPairs[GeneTypes[i]]], [ConservedPairs[GeneTypes[j]], len(HumanAdjacentgenePairs[GeneTypes[j]])-ConservedPairs[GeneTypes[j]]]])[1]))
    newfile.write('\t'.join(Line) + '\n')
newfile.close()
    


# create a function to format the subplots
def CreateAx(NColumns, NRows, Grid1, Grid2, RowPos,ColPos, figure, gs, Data, GraphType):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot using gridspec
    ax = plt.subplot2grid((Grid1,Grid2), (RowPos,ColPos), colspan=NColumns, rowspan=NRows)

    if GraphType == 'heatmap':
        # plot heatmap (use vmin and vmax to get the full range of values)
        heatmap = ax.imshow(Data, interpolation = 'nearest', cmap = 'YlGn')
        # add heatmap scale 
        cbar = plt.colorbar(heatmap)
        # edit tcik parameters of the heatmap scale
        cbar.ax.tick_params(labelsize=7)
        cbar.ax.tick_params(direction = 'out')
        # edit xticks
        plt.xticks([0,1,2,3,4,5,6], ['Not', 'Nst', 'Ext', 'Int', 'Pbk', 'Con', 'Div'])
        plt.yticks([i for i in range(16)], ['Chimp', 'Gorilla', 'Orangutan', 'Macaque',
                    'Marmoset', 'Hedgehog', 'Shrew', 'Cat', 'Dog', 'Mouse', 'Cow', 'Horse',
                    'Sloth', 'Armadillo','Opossum', 'Platypus'])
        # do not show lines around figure  
        ax.spines["top"].set_visible(False)    
        ax.spines["bottom"].set_visible(False)    
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)  
        # edit tick parameters    
        plt.tick_params(axis='both', which='both', bottom='on', top='off',
                        right = 'off', left = 'on', labelbottom='on', labelleft = 'on',
                        colors = 'black', labelsize = 7, direction = 'out')  
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')   
    
    elif GraphType == 'pairs':
        a, b, c = [i[0] for i in Data], [i[1] for i in Data], [i[2] for i in Data]
        # make a list for added values for a and b
        d = [a[i] + b[i] for i in range(len(a))]
        ## Create a bar plot for proportions of conserved gene pairs
        ax.bar([0, 0.3, 0.6, 0.9], a, width = 0.2, label = '2 conserved', color= '#9e9ac8', linewidth = 0.7)
        ax.bar([0, 0.3, 0.6, 0.9], b, width = 0.2, bottom = a, label = '1 conserved', color= '#fd8d3c', linewidth = 0.7)
        ax.bar([0, 0.3, 0.6, 0.9], c, width = 0.2, bottom = d, label = '0 conserved', color= '#78c679', linewidth = 0.7)
        LabelSize = 7
        # set font for all text in figure
        FigFont = {'fontname':'Arial'}   
        # write label for y and x axis
        ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = LabelSize, ha = 'center', **FigFont)
        # write label for x axis
        plt.xticks([0.1, 0.4, 0.7, 1], ['Nst', 'Pbk', 'Con', 'Div'], ha = 'center', fontsize = LabelSize, **FigFont)
        # limit the y axis value range
        plt.ylim([0, 1])   
        # do not show lines around figure  
        ax.spines["top"].set_visible(False)    
        ax.spines["bottom"].set_visible(True)    
        ax.spines["right"].set_visible(False)    
        ax.spines["left"].set_visible(True)  
        # do not show ticks
        plt.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
                        left = 'on', labelbottom='on', colors = 'black', labelsize = LabelSize, direction = 'out')  
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')   
        # create a margin around the x axis
        plt.margins(0.1)
        # add legend
        Two = mpatches.Patch(facecolor = '#9e9ac8' , edgecolor = 'black', linewidth = 0.7, label= '2')
        One = mpatches.Patch(facecolor = '#fd8d3c' , edgecolor = 'black', linewidth = 0.7, label= '1')
        Zero = mpatches.Patch(facecolor = '#78c679' , edgecolor = 'black', linewidth = 0.7, label= '0')
        ax.legend(handles = [Two, One, Zero], loc = (0.1, 1.05), fontsize = LabelSize, frameon = False, ncol = 3)
    elif GraphType == 'distance':
        # make a list of gene category names parallel to the list of gene pairs
        GeneCats = ['Nst', 'Pbk', 'Conv', 'Div', 'Prox', 'Mod', 'Int', 'Dist']
        GeneTypes = ['Nested', 'PiggyBack', 'Convergent', 'Divergent', 'Proximal', 'Moderate', 'Intermediate', 'Distant']
        # set colors
        colorscheme = ['#f03b20', '#43a2ca', '#fee391', '#74c476', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
        # plot proportions of gene pairs
        ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], [Data[i] for i in GeneTypes], 0.2, color = colorscheme,
               edgecolor = 'black', linewidth = 0.7)
        # set font for all text in figure
        FigFont = {'fontname':'Arial'}   
        # write y axis label
        ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = 7, ha = 'center', **FigFont)
        # add ticks and lebels
        plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)

        ## add a range for the Y and X axes
        plt.ylim([0, 1])
        # edit y axis ticks
        plt.yticks(np.arange(0, 1.2, 0.2)) 
        plt.xlim([0, 2.45])
        # do not show lines around figure  
        ax.spines["top"].set_visible(False)    
        ax.spines["bottom"].set_visible(True)    
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(True)  
        # edit tick parameters    
        plt.tick_params(axis='both', which='both', bottom='on', top='off',
                        right = 'off', left = 'on', labelbottom='on',
                        colors = 'black', labelsize = 7, direction = 'out')  
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')   
      
    return ax

# create figure
figure = plt.figure(1, figsize = (6.5, 4))
# set up grid
gs = gridspec.GridSpec(2, 3, width_ratios=[3, 1]) 
# plot data
ax1 = CreateAx(2, 2, 2, 3, 0, 0, figure, gs, Conserved, 'heatmap')
ax2 = CreateAx(1, 1, 2, 3, 0, 2, figure, gs, PairCounts, 'pairs')
ax3 = CreateAx(1, 1, 2, 3, 1, 2, figure, gs, ConservedPairs, 'distance')

# add subplot labels
ax2.text(-2.5, 1.1, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 8)
ax2.text(-0.3, 1.1, 'B', horizontalalignment='center', verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
ax3.text(-0.3, 1.1, 'C', horizontalalignment='center', verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)

# make sure subplots do not overlap
plt.tight_layout()
# save figure
figure.savefig('truc.pdf', bbox_inches = 'tight')
