# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 21:38:52 2017

@author: Richard
"""

# use this script to plot dN and dS between orthologs 

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
        SeqDiv[line[0]] = [float(line[2]), float(line[3])]
infile.close()            
    
# make lists of dN, dS and omega values for each gene category
dN, dS = [], []
for i in range(len(AllPairs)):
    nonsyn, syn = [], []
    for pair in AllPairs[i]:
        if pair[0] in SeqDiv:
            nonsyn.append(SeqDiv[pair[0]][0])
            syn.append(SeqDiv[pair[0]][1])
    dN.append(nonsyn)
    dS.append(syn)
    
 # create lists with means and SEM for dN for each gene category
MeandN, SEMdN = [], []
for i in range(len(dN)):
    MeandN.append(np.mean(dN[i]))
    SEMdN.append(np.std(dN[i]) / math.sqrt(len(dN[i])))
MeandS, SEMdS = [], []
for i in range(len(dS)):
    MeandS.append(np.mean(dS[i]))
    SEMdS.append(np.std(dS[i]) / math.sqrt(len(dS[i])))

# perform statistical tests between gene categories using Kolmogorov-Smirnof test
# create list to store the p-values
PValdN, PValdS = [], []
for i in range(1, len(dN)):
    # compare each gene category to non-overlapping genes
    val, P = stats.ks_2samp(dN[0], dN[i])
    PValdN.append(P)
for i in range(1, len(dS)):
    # compare each gene category to non-overlapping genes
    val, P = stats.ks_2samp(dS[0], dS[i])
    PValdS.append(P)

# replace P values with significance
for i in range(len(PValdN)):
    if PValdN[i] >= 0.05:
        PValdN[i] = ''
    elif PValdN[i] < 0.05 and PValdN[i] >= 0.01:
        PValdN[i] = '*'
    elif PValdN[i] < 0.01 and PValdN[i] >= 0.001:
        PValdN[i] = '**'
    elif PValdN[i] < 0.001:
        PValdN[i] = '***'
for i in range(len(PValdS)):
    if PValdS[i] >= 0.05:
        PValdS[i] = ''
    elif PValdS[i] < 0.05 and PValdS[i] >= 0.01:
        PValdS[i] = '*'
    elif PValdS[i] < 0.01 and PValdS[i] >= 0.001:
        PValdS[i] = '**'
    elif PValdS[i] < 0.001:
        PValdS[i] = '***'



# create a function to format the subplots
def CreateAx(Columns, Rows, Position, fig, Data, XLabel, YLabel, DataType, YMax):
    '''
    Returns a ax instance in figure
    '''    

    # add a plot to figure (N row, N column, plot N)
    ax = fig.add_subplot(Rows, Columns, Position)
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


# make a figure with mean dN and dS 

# create figure
figure = plt.figure(1, figsize = (4.5, 2))
# plot data
plot1 = CreateAx(2, 1, 1, figure, [MeandN, SEMdN], GeneCats, 'Nucleotide divergence (dN)', 'divergence', 0.03)
plot2 = CreateAx(2, 1, 2, figure, [MeandS, SEMdS], GeneCats, 'Nucleotide divergence (dS)', 'divergence', 0.06)

# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
ypos = [0.03, 0.03, 0.03, 0.03, 0.03, 0.03]
xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
for i in range(len(PValdN)):
    plot1.text(xpos[i], ypos[i], PValdN[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
ypos = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
for i in range(len(PValdS)):
    plot2.text(xpos[i], ypos[i], PValdS[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
figure.savefig('truc.pdf', bbox_inches = 'tight')

