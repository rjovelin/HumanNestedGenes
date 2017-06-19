# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 17:05:50 2017

@author: RJovelin
"""


# use this script to plot the proportion of nested gene pairs with and without disease genes for same and opposite strand orientation


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


# load dictionary of nested genes
json_data = open('HumanNestedGenes.json')
Nested = json.load(json_data)
json_data.close()

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

# map ensembl gene IDs to gene names {gene_ID: Name}
GeneIDToNames = MapNametoID('Homo_sapiens.GRCh38.88.gff3')
# reverse dictionary {name: ID}
GeneNamesToID = {}
for ID in GeneIDToNames:
    GeneNamesToID[GeneIDToNames[ID]] = ID

# make a set of complex disease genes
GAD = ParseComplexDisease('GADCDC_data.tsv', GeneNamesToID)
# make a set of GWAS genes 
GWAS = ParseGWASDisease('gwas_catalog_v1.0-associations_e87_r2017-01-09.tsv', 'TraitsToRemove.txt', GeneNamesToID)
# make a set of cancer driver genes
Drivers = ParseCosmicFile('Census_allMon_Mar27_2017.tsv', GeneNamesToID)
# mnake a set of mendelean disease genes
OMIM = ParseOMIMDisease('mimTitles.txt', 'morbidmap.txt', GeneNamesToID)

# create a set with all disease genes
DiseaseGenes = GAD.union(GWAS).union(Drivers).union(OMIM)

# make a list of nested gene pairs
NestedPairs = GetHostNestedPairs(Nested)

# make a list with counts of disease, non disease for same and opposite strand for each disease
PairCounts = []
# loop over disease types
for D in [GAD, GWAS, Drivers, OMIM, DiseaseGenes]:
    SameDisease, SameNoDisease, OppDisease, OppNoDisease = 0, 0, 0, 0
    for pair in NestedPairs:
        # check gene orientation
        Orientation = GenePairOrientation(pair, GeneCoord)
        if len(set(Orientation)) == 1:
            # same strand pair
            # check if at least 1 gene is disease gene
            if pair[0] in D or pair[1] in D:
                SameDisease += 1
            else:
                SameNoDisease += 1
        elif len(set(Orientation)) ==2:
            # opposite strand pair
            # check if at least 1 gene is disease gene
            if pair[0] in D or pair[1] in D:
                OppDisease += 1
            else:
                OppNoDisease += 1
    assert SameDisease + SameNoDisease + OppDisease + OppNoDisease == len(NestedPairs)
    PairCounts.append([[SameDisease, SameNoDisease], [OppDisease, OppNoDisease]])
    
# make a list with Proportions of pairs with disease genes for same and oppoiste strands
ProportionDisease = []
for L in PairCounts:
    print(L[0][0] / sum(L[0]), L[1][0] / sum(L[1]), stats.fisher_exact(L)[1])    
    ProportionDisease.append([L[0][0] / sum(L[0]), L[1][0] / sum(L[1])])
    
    

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel):
    '''
    Returns a ax instance in figure
    '''    
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # Create a horizontal bar plot for proportions of disease genes
    ax.bar([0, 0.2], Data, width = 0.2, label = 'disease', color= ['#2b8cbe', '#88419d'], edgecolor = 'black', linewidth = 0.7)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    YLabel = 'Proportions of pairs\nwith disease genes'
    ax.set_ylabel(YLabel, color = 'black',  size = 6, ha = 'center', **FigFont)
    # edit y axis ticks
#    plt.yticks(np.arange(0, 1, 0.25))    
    # add ticks and lebels
    plt.xticks([0.1, 0.3], ['S', 'O'], rotation = 0, size = 6, color = 'black', ha = 'center', **FigFont)
    # add x label
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # add a range for the Y axis
    # plt.ylim([0, 1])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 3))           
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='on', top='off',
                    right = 'off', left = 'on', labelbottom='on',
                    colors = 'black', labelsize = 6, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    # add margins
    plt.margins(0.1)
    return ax


# make a figure with proportion of disease and non-disease genes

# create figure
fig = plt.figure(1, figsize = (6.2, 1.8))
# plot data
ax1 = CreateAx(5, 1, 1, fig, ProportionDisease[0], 'complex')
ax2 = CreateAx(5, 1, 2, fig, ProportionDisease[1], 'GWAS')
ax3 = CreateAx(5, 1, 3, fig, ProportionDisease[2], 'tumors')
ax4 = CreateAx(5, 1, 4, fig, ProportionDisease[3], 'mendelian')
ax5 = CreateAx(5, 1, 5, fig, ProportionDisease[4], 'all')


# add legend
S = mpatches.Patch(facecolor = '#2b8cbe', edgecolor = 'black', linewidth = 0.7, label= 'same strand')
O = mpatches.Patch(facecolor = '#88419d', edgecolor = 'black', linewidth = 0.7, label= 'opposite strand')
ax1.legend(handles = [S, O], loc = (0.05, 1.2), fontsize = 6, frameon = False, ncol = 2)


# make sure subplots do not overlap
plt.tight_layout()

## save figure to file
#outputfile = 'truc'
#for extension in ['.pdf', '.eps', '.png']:
#    fig.savefig(outputfile + extension, bbox_inches = 'tight')


# save figure to file
fig.savefig('truc.pdf', bbox_inches = 'tight')

