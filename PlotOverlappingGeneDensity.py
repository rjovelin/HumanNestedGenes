# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 13:07:16 2016

@author: Richard
"""


# use this script to plot the density of overlapping genes on each chromo

# usage PlotOverlappingGeneDensity.py [options]
# -[All/NoOverlap]: consider all genes or non-overlapping genes as background
# -[Freq/Numbs]: plot frequencies or counts

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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


# get background gene density from command
# ("All": all genes or "NoOverlap": non-overlapping genes)
BackGround = sys.argv[1]
# get type of values to consider (frequencies or counts)
Vals = sys.argv[2]
# check that options have approporiate values
assert BackGround in ['All', 'NoOverlap']
assert Vals in ['Freq', 'Numbs']

# load dictionary of overlapping gene pairs
json_data = open('HumanOverlappingGenes.json')
Overlapping = json.load(json_data)
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
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)

# make sets of genes
OverlapGenes = MakeFullPartialOverlapGeneSet(Overlapping)

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)

# make a set of non-overlapping genes
NonOverlapGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)


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

# count the number of non-overlapping genes per chromo
NonOverlapCounts = {}
for gene in NonOverlapGenes:
    # get chromo
    chromo = GeneCoord[gene][0]
    if chromo in NonOverlapCounts:
        NonOverlapCounts[chromo] += 1
    else:
        NonOverlapCounts[chromo] = 1

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

# get the start position of each non-overlapping gene on each chromo
NonOverlapStart = {}
for gene in NonOverlapGenes:
    chromo = GeneCoord[gene][0]
    if GeneCoord[gene][-1] == '+':
        start = GeneCoord[gene][1]
    elif GeneCoord[gene][-1] == '-':
        start = GeneCoord[gene][2] - 1
    if chromo in NonOverlapStart:
        NonOverlapStart[chromo].append(start)
    else:
        NonOverlapStart[chromo] = [start]
# sort positions
for chromo in NonOverlapStart:
    NonOverlapStart[chromo].sort()

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
Interval = 1000000 


   
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
# get the count of non-overlapping genes per window
NonOverlapWindowCount = {}
for chromo in Chromosomes:
    counts, binedges = np.histogram(NonOverlapStart[chromo], bins = range(0, len(Genome[chromo]), Interval))
    NonOverlapWindowCount[chromo] = [list(counts), list(binedges)]

# add 0 to match dimensions between bin edges and counts
for chromo in Chromosomes:
    GeneWindowCount[chromo][0].insert(0,0)
    OverlapWindowCount[chromo][0].insert(0,0) 
    NonOverlapWindowCount[chromo][0].insert(0,0)


if Vals == 'Freq':
    # convert counts to frequencies by dividing each count in window by number of genes on chromosome
    for chromo in Chromosomes:
        for i in range(len(GeneWindowCount[chromo][0])):
            GeneWindowCount[chromo][0][i] = GeneWindowCount[chromo][0][i] / GeneCounts[chromo]
        for i in range(len(OverlapWindowCount[chromo][0])):
            OverlapWindowCount[chromo][0][i] = OverlapWindowCount[chromo][0][i] / OverlapCounts[chromo]
        for i in range(len(NonOverlapWindowCount[chromo][0])):
            NonOverlapWindowCount[chromo][0][i] = NonOverlapWindowCount[chromo][0][i] / NonOverlapCounts[chromo]


# get maximum frequency
Maximum = 0
for chromo in GeneWindowCount:
    for i in GeneWindowCount[chromo][0]:
        if i >= Maximum:
            Maximum = i
    for i in OverlapWindowCount[chromo][0]:
        if i >= Maximum:
            Maximum = i
    for i in NonOverlapWindowCount[chromo][0]:
        if i >= Maximum:
            Maximum = i


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, chromo, YMax, YLabel):
    '''
    (int, int, int, figure_object, list, str, float, bool)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of dictionaries with gene densities on each chromosomes
    or the difference between gene frequencies, the chromosome of interest,
    the maximum frequency value, and a boolean indicating whether the Y axis should
    be represented and return a subplot in figure
    '''    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    
    # plot all gene or non-overlapping gene density first
    ax.plot(Data[0][chromo][1], Data[0][chromo][0], linewidth = 1, linestyle = '-', color = 'black', alpha = 0.7)
    # plot overlapping genes second        
    ax.plot(Data[1][chromo][1], Data[1][chromo][0], linewidth = 1, linestyle = '-', color = 'red', alpha = 0.7)
        
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
    plt.xticks([], [])
    
    if YLabel == True:
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')   
    
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
    if BackGround == 'All':
        Data = [GeneWindowCount, OverlapWindowCount]
    elif BackGround == 'NoOverlap':
        Data = [NonOverlapWindowCount, OverlapWindowCount]
    ax = CreateAx(11, 2, j, fig, Data, LG[i], Maximum, YLabel)
    j += 1
    if i == 10:
        # create legend
        if BackGround == 'All':
            Label1 = 'All'
        elif BackGround == 'NoOverlap':
            Label1 = 'Non-overlapping'
        black_line = mlines.Line2D([], [], color='black', marker='', linewidth = 1.2, label = Label1)
        grey_line = mlines.Line2D([], [], color='red', marker='', linewidth = 1.2, label = 'Overlapping')
        plt.legend(handles=[black_line, grey_line], bbox_to_anchor=(-12, 0.8), loc = 3, ncol = 2, fontsize = 10, frameon = False, borderaxespad = 0.)


# make sure subplots do not overlap
plt.tight_layout()

# save figure
if BackGround == 'All' and Vals == 'Freq':
    outputfile = 'AllGenesDensityFreq'
elif BackGround == 'NoOverlap' and Vals == 'Freq':
    outputfile = 'NonOverlapGeneDensityFreq'
elif BackGround == 'All' and Vals == 'Numbs':
    outputfile = 'AllGenesDensityCounts'
elif BackGround == 'NoOverlap' and Vals == 'Numbs':
    outputfile = 'NonOverlapGeneDensityCounts'
     
fig.savefig(outputfile + '.pdf', bbox_inches = 'tight')
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')