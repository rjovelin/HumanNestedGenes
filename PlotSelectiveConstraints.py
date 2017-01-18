# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 11:58:03 2017

@author: Richard
"""

# use this script to plot dN/dS and between orthologs and proportion of genes with homologs


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

# create sets of internal and external nested gene pairs
NestedPairs = GetHostNestedPairs(Nested)
InternalGenes, ExternalGenes = set(), set()
for pair in NestedPairs:
    ExternalGenes.add(pair[0])
    InternalGenes.add(pair[1])
    
# get 1:1 orthologs between human and chimp
Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')

# create lists of orthologous pairs for each gene category 
GeneCats = ['NoOv', 'Nst', 'Int', 'Ext', 'Pbk', 'Conv', 'Div']
AllPairs = []
AllGenes = [NonOverlappingGenes, NestedGenes, InternalGenes, ExternalGenes,
            PiggyBackGenes, ConvergentGenes, DivergentGenes] 
# loop over gene sets
for i in range(len(AllGenes)):
    # create a list of gene pairs
    orthologs = []
    # loop over genes in given gene set
    for gene in AllGenes[i]:
        # check that gene has ortholog
        if gene in Orthos:
            orthologs.append([gene, Orthos[gene]])
    AllPairs.append(orthologs)

    
# create a dict with divergence values
SeqDiv = {}
infile = open('HumanChimpSeqDiverg.txt')
infile.readline()
for line in infile:
    if line.startswith('ENSG'):
        line = line.rstrip().split('\t')
        if line[-1] != 'NA':
            SeqDiv[line[0]] = float(line[-1])
infile.close()            
    
# make list of dN/dS for each gene category
Omega = []
for i in range(len(AllPairs)):
    ratio = [] 
    for pair in AllPairs[i]:
        if pair[0] in SeqDiv:
            ratio.append(SeqDiv[pair[0]])
    Omega.append(ratio)
    
 # create lists with means and SEM for dN/dS for each gene category
MeanOmega, SEMOmega = [], []
for i in range(len(Omega)):
    MeanOmega.append(np.mean(Omega[i]))
    SEMOmega.append(np.std(Omega[i]) / math.sqrt(len(Omega[i])))

# perform statistical tests between gene categories using Kolmogorov-Smirnof test
# create list to store the p-values
PValOmega = []
for i in range(1, len(Omega)):
    # compare each gene category to non-overlapping genes
    val, P = stats.ks_2samp(Omega[0], Omega[i])
    PValOmega.append(P)
# replace P values with significance
for i in range(len(PValOmega)):
    if PValOmega[i] >= 0.05:
        PValOmega[i] = ''
    elif PValOmega[i] < 0.05 and PValOmega[i] >= 0.01:
        PValOmega[i] = '*'
    elif PValOmega[i] < 0.01 and PValOmega[i] >= 0.001:
        PValOmega[i] = '**'
    elif PValOmega[i] < 0.001:
        PValOmega[i] = '***'

# create a set of human genes that have homologs
# include 1:1 and 1 to many or many to many orthologs
Homologs = set()
infile = open('HumanChimpOrthologs.txt')
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



# make a figure with mean dN/dS and with proportion of gene with homologs

# create figure
fig = plt.figure(1, figsize = (4.5, 2))
# plot data
ax1 = CreateAx(2, 1, 1, fig, [MeanOmega, SEMOmega], GeneCats, 'Nucleotide divergence (dN/dS)', 'divergence', 0.50)
ax2 = CreateAx(2, 1, 2, fig, [WithHomolog, NoHomolog], GeneCats, 'Proportion', 'proportion', 1)

# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
ypos = [0.47, 0.50, 0.47, 0.47, 0.45, 0.45]
xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
for i in range(len(PValOmega)):
    ax1.text(xpos[i], ypos[i], PValOmega[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
for i in range(len(PProp)):
    ax2.text(xpos[i], 1.02, PValOmega[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)

# add legend
NoH = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'no homolog')
WiH = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'homolog')
ax2.legend(handles = [WiH, NoH], loc = (0, 1.1), fontsize = 6, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('SelectiveConstraints.pdf', bbox_inches = 'tight')
fig.savefig('SelectiveConstraints.eps', bbox_inches = 'tight')