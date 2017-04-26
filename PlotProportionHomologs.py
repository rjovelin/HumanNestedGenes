# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 22:00:06 2017

@author: Richard
"""


# use this script to plot the proportion of genes with homologs 

# usage python3 PlotProportionHomologs.py [options]
# -[chimp/mouse]: use human-chimp or human-mouse comparisons

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

Species = sys.argv[1]
assert Species in ['chimp', 'mouse']


# load dictionaries of overlapping genes
JsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',
             'HumanPiggyBackGenes.json', 'HumanConvergentGenes.json',
             'HumanDivergentGenes.json']
# make a list of dictionaries
Overlap = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    Overlap.append(overlapping)

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

# generate gene sets
OverlappingGeneSets = []
for i in range(len(Overlap)):
    OverlappingGeneSets.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

# create sets of internal and external nested gene pairs
NestedPairs = GetHostNestedPairs(Overlap[1])
InternalGenes, ExternalGenes = set(), set()
for pair in NestedPairs:
    ExternalGenes.add(pair[0])
    InternalGenes.add(pair[1])

# get orthologs
if Species == 'chimp':
    Orthos = MatchOrthologs('HumanChimpOrthologs.txt')
elif Species == 'mouse':
    Orthos = MatchOrthologs('HumanMouseOrthologs.txt')

# create a list of gene categories
GeneCats = ['Not', 'Nst', 'Int', 'Ext', 'Pbk', 'Con', 'Div']

# create lists of orthologous pairs for each gene category 
AllPairs = []
AllGenes = [NonOverlappingGenes, OverlappingGeneSets[1], InternalGenes, ExternalGenes,
            OverlappingGeneSets[2], OverlappingGeneSets[3], OverlappingGeneSets[4]] 
# loop over gene sets
for i in range(len(AllGenes)):
    # create a list of gene pairs
    orthologs = []
    # loop over genes in given gene set
    for gene in AllGenes[i]:
        # check that gene has ortholog
        if gene in Orthos:
            for homolog in Orthos[gene]:
                orthologs.append([gene, homolog])
    AllPairs.append(orthologs)


# compare proportions of genes with homologs for each category

# create a set of human genes that have any homologs
Homologs = set()
if Species == 'chimp':
    infile = open('HumanChimpOrthologs.txt')
elif Species == 'mouse':
    infile = open('HumanMouseOrthologs.txt')
infile.readline()
for line in infile:
    if 'ortholog' in line:
        line = line.rstrip().split('\t')
        assert 'ENS' in line[0], 'gene id is not valid'
        assert 'ortholog' in line[4], 'ortholog should be in homology type'
        Homologs.add(line[0])
infile.close()            

# count genes with and without homologs
GeneCounts = []
# loop over gene sets
for i in range(len(AllGenes)):
    # initialize counters
    homo, nohomo = 0, 0
    # loop over genes in given set
    for gene in AllGenes[i]:
        if gene in Homologs:
            homo += 1
        else:
            nohomo += 1
    GeneCounts.append([homo, nohomo])    

# compare the proportions of gene with and without homologs
# create a list to store the P-values
PProp = []
for i in range(1, len(GeneCounts)):
    p = stats.fisher_exact([GeneCounts[0], GeneCounts[i]])[1]
    PProp.append(p)
# replace P values by significance
PProp = ConvertPToStars(PProp)

# get the proportions of genes with and without homologs
WithHomolog, NoHomolog = [], []
for i in range(len(GeneCounts)):
    WithHomolog.append(GeneCounts[i][0] / sum(GeneCounts[i]))
    NoHomolog.append(GeneCounts[i][1] / sum(GeneCounts[i]))
    assert sum(GeneCounts[i]) == len(AllGenes[i])

# plot sequence divergence, proportion of homologs, expression divergence in a single figure
fig = plt.figure(1, figsize = (2.2, 1.8))

# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)

# Create a horizontal bar plot for proportions of genes with homologs
ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], WithHomolog, width = 0.2, label = 'homolog', color= '#dadaeb', linewidth = 0.7)
# Create a horizontal bar plot for proportions of same strand pairs
ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], NoHomolog, width = 0.2, bottom = WithHomolog, label = 'no homolog', color= '#d9f0a3', linewidth = 0.7)

# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Proportion', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9], GeneCats, rotation = 0, size = 7, color = 'black', ha = 'right', **FigFont)
# add a range for the Y and X axes
plt.ylim([0, 1])    
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
# add margins
plt.margins(0.1)

# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
for i in range(len(PProp)):
    ax.text(xpos[i], 1.02, PProp[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)

# add legend
NoH = mpatches.Patch(facecolor = '#d9f0a3' , edgecolor = 'black', linewidth = 0.7, label= 'no homolog')
WiH = mpatches.Patch(facecolor = '#dadaeb' , edgecolor = 'black', linewidth = 0.7, label= 'homolog')
ax.legend(handles = [WiH, NoH], loc = (0, 1.1), fontsize = 6, frameon = False, ncol = 2)

# build figure name with option
FigName = 'ProportionsHomologs' + Species.title()
# save figure
fig.savefig(FigName + '.pdf', bbox_inches = 'tight')
fig.savefig(FigName + '.eps', bbox_inches = 'tight')

