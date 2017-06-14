# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 18:58:58 2017

@author: RJovelin
"""


#  use this script to plot conservation of overlapping pairs in mammals

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


# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan',
           'Macaque', 'Marmoset', 'Hedgehog', 'Shrew',
           'Cat', 'Dog', 'Mouse', 'Cow', 'Horse', 'Sloth',
           'Armadillo', 'Opossum', 'Platypus']

# make a list of json files
JsonFiles = ['Overlapping', 'Nested', 'PiggyBack', 'Convergent', 'Divergent']

# make a list of lists of dictionaries for each type of overlapping gene in each species
AllOverlap = []
# loop over files
for i in range(len(Species)):
    ovlp = []
    for j in range(len(JsonFiles)):
        # load dictionary from json file
        json_data = open(Species[i] + JsonFiles[j] + 'Genes.json')
        ovlp.append(json.load(json_data))
        json_data.close()
    AllOverlap.append(ovlp)

# make a parallel list of ortholog files
OrthoFiles = ['Human' + i + 'Orthologs.txt' for i in Species[1:]]

# get coordinates of human genes
GFF = 'Homo_sapiens.GRCh38.88.gff3'
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
HumanCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    
# make a list of lists of pairs of genes
AllGenePairs = []
for i in range(len(AllOverlap)):
    sppairs = []
    for j in range(len(AllOverlap[i])):
        sppairs.append(GetHostNestedPairs(AllOverlap[i][j]))
    AllGenePairs.append(sppairs)

# extract human overlapping gene pairs
HumanGenePairs, HumanAllOverlap = AllGenePairs.pop(0), AllOverlap.pop(0)
 
assert len(AllOverlap) == len(AllGenePairs) == len(OrthoFiles) == len(Species[1:])


# 1) plot the proportion of gene pairs that are overlapping in each species 

# make a list with counts of conserved pairs in each species for overlapping genes defined by strand orientation
ConservedOrientation = []

# make a list of gene pairs by pooling same and opposite strand overlapping gene pairs [[same], [opposite]]
HumanOrientation = [[], []]
for i in range(1, len(HumanAllOverlap)):
    # make a list of human pairs
    pairs = GetHostNestedPairs(HumanAllOverlap[i])
    if i == 1:
        # check if nested pairs are same or opposite strand
        for L in pairs:
            if len(set(GenePairOrientation(L, HumanCoord))) == 2:
                HumanOrientation[1].append(L)
            elif len(set(GenePairOrientation(L, HumanCoord))) == 1:
                HumanOrientation[0].append(L)
    elif i == 2:
        # add all pbk pairs to same strand list
        for L in pairs:
            HumanOrientation[0].append(L)
    elif i == 3 or i == 4:
        # add convergent or divergent pairs to opposite strand list
        for L in pairs:
            HumanOrientation[1].append(L)

# make a list of all overlapping gene pairs for each species
SpeciesOverlap = []
for i in range(len(AllOverlap)):
    SpeciesOverlap.append(GetHostNestedPairs(AllOverlap[i][0]))
    
# loop over species
for i in range(len(SpeciesOverlap)):
    # make a list with coutns of conserved pairs for each orientation [same, opposite]
    ConservedPairs = [0, 0]
    # get the orthologs for that species
    Orthos = MatchOrthologs(OrthoFiles[i])
    # reverse dictionary
    SpeciesOrthos = {}
    for gene in Orthos:
        for ortho in Orthos[gene]:
            if ortho not in SpeciesOrthos:
                SpeciesOrthos[ortho] = [gene]
            else:
                SpeciesOrthos[ortho].append(gene)
    # make pairs of human orthologs for each species gene pairs
    PairsOrthos = []            
    for pair in SpeciesOverlap[i]:
        # count only pairs in which both genes have orthos in human
        if pair[0] in SpeciesOrthos and pair[1] in SpeciesOrthos:
            for ortho1 in SpeciesOrthos[pair[0]]:
                for ortho2 in SpeciesOrthos[pair[1]]:
                    # remove order
                    PairsOrthos.append(set([ortho1, ortho2]))
    # loop over lists of human pairs
    for j in range(len(HumanOrientation)):
        # make a copy of the list iof human pairs
        HumanPairs = copy.deepcopy(HumanOrientation[j])
        # remove pairs if any gene is lacking an ortholog
        to_remove = [L for L in HumanPairs if L[0] not in Orthos or L[1] not in Orthos]
        for L in to_remove:
            HumanPairs.remove(L)
        # check if human pairs are conserved
        for pair in HumanPairs:
            if set(pair) in PairsOrthos:
                # update counter for the given gene class
                ConservedPairs[j] += 1
        # compute proportion
        ConservedPairs[j] = ConservedPairs[j] / len(HumanPairs)
    # populate l;ist for the given species
    ConservedOrientation.append(ConservedPairs)
        
# convert list to numpy array
ConservedOrientation = np.array(ConservedOrientation) 
    


# 2) plot the proportion of gene pairs overlapping in each otehr species

ConservedOverlap = []
# keep gene pairs for the 4 overlapping classes
HumanGenePairs = HumanGenePairs[1:]

# loop over species
for i in range(len(SpeciesOverlap)):
    # make a list with coutns of conserved pairs for each type of overlapping genes
    ConservedPairs = [0] * len(HumanGenePairs)
    # get the orthologs for that species
    Orthos = MatchOrthologs(OrthoFiles[i])
    # reverse dictionary
    SpeciesOrthos = {}
    for gene in Orthos:
        for ortho in Orthos[gene]:
            if ortho not in SpeciesOrthos:
                SpeciesOrthos[ortho] = [gene]
            else:
                SpeciesOrthos[ortho].append(gene)
    # make pairs of human orthologs for each species gene pairs
    PairsOrthos = []            
    for pair in SpeciesOverlap[i]:
        # count only pairs in which both genes have orthos in human
        if pair[0] in SpeciesOrthos and pair[1] in SpeciesOrthos:
            for ortho1 in SpeciesOrthos[pair[0]]:
                for ortho2 in SpeciesOrthos[pair[1]]:
                    # remove order
                    PairsOrthos.append(set([ortho1, ortho2]))
    # loop over lists of human pairs
    for j in range(len(HumanGenePairs)):
        # make a copy of the list iof human pairs
        HumanPairs = copy.deepcopy(HumanGenePairs[j])
        # remove pairs if any gene is lacking an ortholog
        to_remove = [L for L in HumanPairs if L[0] not in Orthos or L[1] not in Orthos]
        for L in to_remove:
            HumanPairs.remove(L)
        # check if human pairs are conserved
        for pair in HumanPairs:
            if set(pair) in PairsOrthos:
                # update counter for the given gene class
                ConservedPairs[j] += 1
        # compute proportion
        ConservedPairs[j] = ConservedPairs[j] / len(HumanPairs)
    # populate l;ist for the given species
    ConservedOverlap.append(ConservedPairs)
        
# convert list to numpy array
ConservedOverlap = np.array(ConservedOverlap) 


  
# create a function to format the subplots
def CreateAx(NColumns, NRows, Grid1, Grid2, RowPos,ColPos, figure, gs, Data, GraphType):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot using gridspec
    ax = plt.subplot2grid((Grid1,Grid2), (RowPos,ColPos), colspan=NColumns, rowspan=NRows)

    if GraphType == 'pairs':
        # plot heatmap (use vmin and vmax to get the full range of values)
        heatmap = ax.imshow(Data, interpolation = 'nearest', cmap = 'Blues', vmin = 0, vmax = 0.5)
        # add heatmap scale 
        cbar = plt.colorbar(heatmap)
        # edit tcik parameters of the heatmap scale
        cbar.ax.tick_params(labelsize=7)
        cbar.ax.tick_params(direction = 'out')
        # set font for all text in figure
        FigFont = {'fontname':'Arial'}
        # edit xticks
        plt.xticks([0,1,2,3], ['Nst', 'Pbk', 'Con', 'Div'], size = 7, color = 'black', ha = 'center', **FigFont)
    elif GraphType == 'orientation':
        # plot heatmap (use vmin and vmax to get the full range of values)
        heatmap = ax.imshow(Data, interpolation = 'nearest', cmap = 'Reds', vmin = 0, vmax = 0.5)
        # add heatmap scale 
        cbar = plt.colorbar(heatmap)
        # edit tcik parameters of the heatmap scale
        cbar.ax.tick_params(labelsize=7)
        cbar.ax.tick_params(direction = 'out')
        # set font for all text in figure
        FigFont = {'fontname':'Arial'}
        # edit xticks
        plt.xticks([0,1], ['S', 'O'], size = 7, color = 'black', ha = 'center', **FigFont)
        # add x axis label
        ax.set_xlabel('Orientation', color = 'black', size = 7, ha = 'center', **FigFont)  
   
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
    return ax

# create figure
figure = plt.figure(1, figsize = (5, 4))
# set up grid
gs = gridspec.GridSpec(1, 2) 
# plot data
ax1 = CreateAx(1, 1, 1, 2, 0, 0, figure, gs, ConservedOverlap, 'pairs')
ax2 = CreateAx(1, 1, 1, 2, 0, 1, figure, gs, ConservedOrientation, 'orientation')

## add subplot labels
ax2.text(-16, -1, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 8)
ax2.text(-3.5, -1, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 8)

# make sure subplots do not overlap
plt.tight_layout()
# save figure
figure.savefig('truc.pdf', bbox_inches = 'tight')
