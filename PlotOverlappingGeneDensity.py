# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 13:07:16 2016

@author: Richard
"""


# use this script to plot the density of overlapping genes on each chromo

# usage PlotOverlappingGeneDensity.py

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
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
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)

# make sets of genes
OverlapGenes = MakeFullPartialOverlapGeneSet(Overlapping)

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)

# get the gene counts on each chromo
GeneCounts = {}
for chromo in GeneChromoCoord:
    GeneCounts[chromo] = len(GeneChromoCoord[chromo])

# count the number of overlapping genes per chromo
OverlapCounts = {}
for gene in OverlapGenes:
    # get chromo
    chromo = GeneCoord[gene][0]
    if chromo in OverlapCounts:
        OverlapCounts[chromo] += 1
    else:
        OverlapCounts[chromo] = 1

# get frequencies of overlapping genes on each chromosome
Freq = {}
for chromo in OverlapCounts:
    Freq[chromo] = round((OverlapCounts[chromo] / len(OverlapGenes)) * 100, 4)

# create a dictionary with the chromosome sequences {chromo: dna}
Genome = {}    
infile = open('Homo_sapiens.GRCh38.dna.primary_assembly.fa')
content = infile.read().rstrip()
content = content.split('>')
content.remove('')
for i in range(len(content)):
    content[i] = content[i].split('\n')
for i in range(len(content)):
    chromo = content[i][0].split()[0]
    dna = ''.join(content[i][1:])
    Genome[chromo] = dna

  
# set up minimum overlapping gene frequency (in %)
MinimumFreq = 1
# consider chromosomes with frequency of overlapping genes > minimum frequency
Chromosomes = [chromo for chromo in Freq if Freq[chromo] > MinimumFreq]
Chromosomes.sort()

# look for clusters of overlapping genes on each chromo
# create a dict {chromo: [positions]}
OverlapStart = {}
for gene in OverlapGenes:
    chromo = GeneCoord[gene][0]
    if GeneCoord[gene][-1] == '+':
        start = GeneCoord[gene][1]
    elif GeneCoord[gene][-1] == '-':
        start = GeneCoord[gene][2] - 1
    if chromo in OverlapStart:
        OverlapStart[chromo].append(start)
    else:
        OverlapStart[chromo] = [start]
# sort positions
for chromo in OverlapStart:
    OverlapStart[chromo].sort()

# get the start position of each gene on chromo
GeneStart = {}
for gene in GeneCoord:
    chromo = GeneCoord[gene][0]
    # only consider chromosomes with enough overlapping genes
    if chromo in Chromosomes:
        # check orientation
        if GeneCoord[gene][-1] == '+':
            start = GeneCoord[gene][1]
        elif GeneCoord[gene][-1] == '-':
            start = GeneCoord[gene][2] - 1
        if chromo in GeneStart:
            GeneStart[chromo].append(start)
        else:
            GeneStart[chromo] = [start]
# sort positions
for chromo in GeneStart:
    GeneStart[chromo].sort()

    
# set up interval length in bp
Interval = int(10000000 / 2)


   
# get the count of genes per window
GeneWindowCount = {}
for chromo in Chromosomes:
    counts, binedges = np.histogram(GeneStart[chromo], bins = range(0, len(Genome[chromo]), Interval))
    GeneWindowCount[chromo] = [list(counts), list(binedges)]
# get the count of overlaps per window
OverlapWindowCount = {}
for chromo in Chromosomes:
    counts, binedges = np.histogram(OverlapStart[chromo], bins = range(0, len(Genome[chromo]), Interval))
    OverlapWindowCount[chromo] = [list(counts), list(binedges)]

# add 0 to match dimensions between bin edges and counts
for chromo in Chromosomes:
    GeneWindowCount[chromo][0].insert(0,0)
    OverlapWindowCount[chromo][0].insert(0,0) 

# convert counts to frequencies by dividing each count in window by number of genes on chromosome
for chromo in Chromosomes:
    for i in range(len(GeneWindowCount[chromo][0])):
        GeneWindowCount[chromo][0][i] = GeneWindowCount[chromo][0][i] / GeneCounts[chromo]
    for i in range(len(OverlapWindowCount[chromo][0])):
        OverlapWindowCount[chromo][0][i] = OverlapWindowCount[chromo][0][i] / OverlapCounts[chromo]


# get maximum frequency
Maximum = 0
for chromo in GeneWindowCount:
    for i in GeneWindowCount[chromo][0]:
        if i >= Maximum:
            Maximum = i
    for i in OverlapWindowCount[chromo][0]:
        if i >= Maximum:
            Maximum = i


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, GeneWindowCount, OverlapWindowCount, chromo, YMax, YLabel):
    '''
    (int, int, int, figure_object, dict, dict, str, float, bool)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, the dictionaries with gene densities on each chromosomes
    for all genes and for overlapping genes only, the chromosome of interest,
    the maximum frequency value, and a boolean indicating whether the Y axis should
    be represented and return a subplot in figure
    '''    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot gene density per window
    ax.plot(GeneWindowCount[chromo][1], GeneWindowCount[chromo][0], linewidth = 1, linestyle = '-', color = 'black', alpha = 1)
    ax.plot(OverlapWindowCount[chromo][1], OverlapWindowCount[chromo][0], linewidth = 1, linestyle = '-', color = 'grey', alpha = 1)
       
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # set axis labels
    if YLabel == True:
        ax.set_ylabel('Gene density', size = 10, ha = 'center', **FigFont)
    ax.set_xlabel(chromo, size = 10, color = 'black', ha = 'center', **FigFont )        
        
    # add a range for the Y axis
    plt.ylim([0, YMax])
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False) 
    if YLabel == True:
        ax.spines["left"].set_visible(True)  
    elif YLabel == False:
        ax.spines["left"].set_visible(False)
    
    if YLabel == True:
        # edit tick paramters
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'on',          
            labelbottom='off', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 10,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
           
    elif YLabel == False:
        # edit tick paramters
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='off', # labels along the bottom edge are on
            labelleft = 'off',            
            colors = 'black',
            labelsize = 8,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
           
    # set x axis ticks
    # plt.xticks(GeneWindowCount[chromo][1], list(map(lambda x: str(x), GeneWindowCount[chromo][1])), rotation = 0, fontsize = 8, fontname = 'Arial')
    plt.xticks([], [])
    
    if YLabel == True:
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')   
    # create a margin around the x axis
    #plt.margins(0.1)
    return ax      


# create figure
fig = plt.figure(1, figsize = (10, 5))

j = 1
LG = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '19', '20', '21', '22', 'X']
for i in range(len(LG)):
    if i == 0 or i == 11:
        YLabel = True
    else:
        YLabel = False
    ax = CreateAx(11, 2, j, fig, GeneWindowCount, OverlapWindowCount, LG[i], Maximum, YLabel)
    j += 1

# make sure subplots do not overlap
plt.tight_layout()

fig.savefig('truc.pdf', bbox_inches = 'tight')