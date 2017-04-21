# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 11:58:03 2017

@author: Richard
"""

# use this script to plot dN/dS and between orthologs and proportion of genes with homologs
# usage python3 PlotSelectiveConstraints.py [options]
# -[chimp/mouse]: compute divergence between human and mouse or chimp orthologs



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
GeneCats = ['NoOv', 'Nst', 'Int', 'Ext', 'Pbk', 'Conv', 'Div']
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

    
# create a dict with divergence values {human_gene: {ortholog: divergence}}
SeqDiv = {}
if Species == 'chimp':
    infile = open('HumanChimpSeqDiverg.txt')
elif Species == 'mouse':
    infile = open('HumanMouseSeqDiverg.txt')
infile.readline()
for line in infile:
    if line.startswith('ENSG'):
        line = line.rstrip().split('\t')
        if line[-1] != 'NA':
            if line[0] not in SeqDiv:
                SeqDiv[line[0]] = {}
            SeqDiv[line[0]][line[1]] = float(line[-1])
infile.close()            
    
# make list of divergence for each gene category
Divergence = []
for i in range(len(AllPairs)):
    nucldiv = []
    for pair in AllPairs[i]:
        if pair[0] in SeqDiv:
            if pair[1] in SeqDiv[pair[0]]:
                nucldiv.append(SeqDiv[pair[0]][pair[1]])
    Divergence.append(nucldiv)
    
 # create lists with means and SEM for divergence for each gene category
MeanDiverg, SEMDiverg = [], []
for i in range(len(Divergence)):
    MeanDiverg.append(np.mean(Divergence[i]))
    SEMDiverg.append(np.std(Divergence[i]) / math.sqrt(len(Divergence[i])))

# perform statistical tests between gene categories using permutation tests
# create list to store the p-values
PValDiverg = []
for i in range(1, len(Divergence)):
    # compare each gene category to non-overlapping genes
    P = PermutationResampling(Divergence[0], Divergence[i], 1000, statistic = np.mean)
    PValDiverg.append(P)
# replace P values with significance
for i in range(len(PValDiverg)):
    if PValDiverg[i] >= 0.05:
        PValDiverg[i] = ''
    elif PValDiverg[i] < 0.05 and PValDiverg[i] >= 0.01:
        PValDiverg[i] = '*'
    elif PValDiverg[i] < 0.01 and PValDiverg[i] >= 0.001:
        PValDiverg[i] = '**'
    elif PValDiverg[i] < 0.001:
        PValDiverg[i] = '***'

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
for i in range(len(PProp)):
    if PProp[i] >= 0.05:
        PProp[i] = ''
    elif PProp[i] < 0.05 and PProp[i] >= 0.01:
        PProp[i] = '*'
    elif PProp[i] < 0.01 and PProp[i] >= 0.001:
        PProp[i] = '**'
    elif PProp[i] < 0.001:
        PProp[i] = '***'

# get the proportions of genes with and without homologs
WithHomolog, NoHomolog = [], []
for i in range(len(GeneCounts)):
    WithHomolog.append(GeneCounts[i][0] / sum(GeneCounts[i]))
    NoHomolog.append(GeneCounts[i][1] / sum(GeneCounts[i]))
    assert sum(GeneCounts[i]) == len(AllGenes[i])


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YLabel, DataType, YMax):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # check type of graphic    
    if DataType == 'divergence':
        # set colors
        colorscheme = ['black','lightgrey','lightgrey','lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
        # plot nucleotide divergence
        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[0], 0.2, yerr = Data[1], color = colorscheme,
               edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    elif DataType == 'proportion':
        ## Create a horizontal bar plot for proportions of genes with homologs
        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[0], width = 0.2, label = 'homolog', color= 'black', linewidth = 0.7)
        # Create a horizontal bar plot for proportions of same strand pairs
        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[1], width = 0.2, bottom = Data[0], label = 'no homolog', color= 'lightgrey', linewidth = 0.7)

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # add ticks and lebels
    plt.xticks([0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9], XLabel, rotation = 30, size = 7, color = 'black', ha = 'right', **FigFont)
    # add a range for the Y and X axes
    plt.ylim([0, YMax])    
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
    return ax


# make a figure with mean divergence and with proportion of gene with homologs
# create figure
fig = plt.figure(1, figsize = (4.5, 2))
# plot data
ax1 = CreateAx(2, 1, 1, fig, [MeanDiverg, SEMDiverg], GeneCats, 'Nucleotide divergence (dN/dS)', 'divergence', 0.50)
ax2 = CreateAx(2, 1, 2, fig, [WithHomolog, NoHomolog], GeneCats, 'Proportion', 'proportion', 1)

# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
ypos = [0.47, 0.50, 0.47, 0.47, 0.45, 0.45]
xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
for i in range(len(PValOmega)):
    ax1.text(xpos[i], ypos[i], PValDiverg[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
for i in range(len(PProp)):
    ax2.text(xpos[i], 1.02, PProp[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)

# add legend
NoH = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'no homolog')
WiH = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'homolog')
ax2.legend(handles = [WiH, NoH], loc = (0, 1.1), fontsize = 6, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('SelectiveConstraints.pdf', bbox_inches = 'tight')
fig.savefig('SelectiveConstraints.eps', bbox_inches = 'tight')