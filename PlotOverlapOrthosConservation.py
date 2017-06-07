# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 17:39:43 2017

@author: RJovelin
"""

# use this script to plot the proportion of gene in each overlapping gene category
# that has orthologs in each other species


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

# create figure
figure = plt.figure(1, figsize = (5, 5))
# add a plot to figure (N row, N column, plot N)
ax = figure.add_subplot(1, 1, 1)

# plot heatmap (use vmin and vmax to get the full range of values)
heatmap = ax.imshow(Conserved, interpolation = 'nearest', cmap = 'YlGn')
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
   
# save figure
figure.savefig('truc.pdf', bbox_inches = 'tight')


#https://pythonhosted.org/DendroPy/primer/trees.html
#http://huboqiang.cn/2016/02/13/PyHeatMapHcl
#https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib
#http://seaborn.pydata.org/generated/seaborn.heatmap.html
#https://matplotlib.org/examples/pylab_examples/pcolor_small.html