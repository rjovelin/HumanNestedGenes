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
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)

# make sets of genes
OverlapGenes = MakeFullPartialOverlapGeneSet(Overlapping)
NestedGenes = MakeFullPartialOverlapGeneSet(Nested)
PiggybackGenes = MakeFullPartialOverlapGeneSet(Piggyback)
ConvergentGenes = MakeFullPartialOverlapGeneSet(Convergent)
DivergentGenes = MakeFullPartialOverlapGeneSet(Divergent)

# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)
NestedPairs = GetHostNestedPairs(Nested)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)


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
    
    
   
for chromo in Freq:
    if Freq[chromo] > 1:
        print(chromo, Freq[chromo], str(len(Genome[chromo])) + ' bp', str(len(Genome[chromo]) // 1000000) + ' Mb', str(len(Genome[chromo]) // 100000) + ' per 100Kb window')



# make a list of chromosomes with overlapping genes
Chromosomes = set()
for pair in OverlappingPairs:
    Chromosomes.add(GeneCoord[pair[0]][0])
Chromosomes = list(Chromosomes)    
Chromosomes.sort()

# look for clusters o overlapping genes on each chromo

# get the first position overlapping between pairs of genes
# create a dict {chromo: [positions]}
OverlapStart = {}
for pair in OverlappingPairs:
    chromo = GeneCoord[pair[0]][0]
    coord1 = set(range(GeneCoord[pair[0]][1], GeneCoord[pair[0]][2]))
    coord2 = set(range(GeneCoord[pair[1]][1], GeneCoord[pair[1]][2]))
    common = list(coord1.intersection(coord2))
    common.sort()
    start = common[0]
    if chromo in OverlapStart:
        OverlapStart[chromo].append(start)
    else:
        OverlapStart[chromo] = [start]

# sort positions
for chromo in OverlapStart:
    OverlapStart[chromo].sort()
















######## continue here


















   
## create a function to count the number of chemo per window interval
#def CountChemoWindow(chemo_start, genome, chromo, window):
#    '''
#    (dict, str, int) -> list
#    Take the dictionary with chemo gene start positions per chromo, 
#    the dict with genome sequence, a chromosome and a window interval in bp
#    and return a list with the number of chemo genes on chromo per window interval'
#    '''
#    # create list with count of gene o repeat in 100000 bp windows
#    range_counts = [0] * (len(genome[chromo]) // window)
#    for start in chemo_start[chromo]:
#        which_range = start // window
#        if which_range == len(range_counts):
#            which_range -= 1
#        # count repeats
#        range_counts[which_range] += 1
#    return range_counts
#    
#    
## set up interval length in bp
#Interval = 100000
#print('Interval window: {0} bp'.format(Interval))    
#
## get the count of chemo gene per window
#WindowCount = {}
#for chromo in chemo_start:
#    range_counts = CountChemoWindow(chemo_start, genome, chromo, Interval)     
#    WindowCount[chromo] = range_counts
#print('got chemo count per window')    
#
## create a list with the position of each window interval
#Positions = {}
#for chromo in WindowCount:
#    pos = [i for i in range(len(WindowCount[chromo]))] 
#    Positions[chromo] = pos     
#    print('position', chromo, len(pos))
#    print('interval length', len(genome[chromo]) // Interval)
#print('got positions of window intervals')
#
#
## create figure
#fig = plt.figure(1, figsize = (4, 1))
## add a plot to figure (1 row, 1 column, 1 plot)
#ax = fig.add_subplot(1, 1, 1)  
#
## find the longtest chromo
#chromoLength = [[len(genome[chromo]), chromo] for chromo in HighFreq]
#chromoLength.sort()
#size = [chromoLength[i][0] for i in range(len(chromoLength))]
#LG = [chromoLength[i][1] for i in range(len(chromoLength))]
#longest, maxlength = chromoLength[-1][-1], chromoLength[-1][0] 
#
## create a list of colors
##colorscheme = ['#a6cee3','#1f78b4','#b2df8a','#33a02c']
#colorscheme = ['#1b9e77','#d95f02','#7570b3']
#
#Graph = {}
## loop over chromo, from chromo with lowest to highest count
#for i in range(len(LG)):
#    print(LG[i])
#    # plot the repeat of gene density per window
#    graph = ax.plot(Positions[LG[i]], WindowCount[LG[i]], linewidth = 1.2, color = colorscheme[i], alpha = 0.7)
#    Graph[LG[i]] = graph
#    
#ax.set_ylabel('GPCRs /100 Kb', size = 10, ha = 'center', fontname = 'Arial')
# 
## set x axis label
#ax.set_xlabel('Position along linkage group (Mb)', size = 10, ha = 'center', fontname = 'Arial')
#
## do not show lines around figure, keep bottow line  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)    
#ax.spines["left"].set_visible(False)      
## offset the spines
#for spine in ax.spines.values():
#  spine.set_position(('outward', 5))
#  
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
## hide these grids behind plot objects
#ax.set_axisbelow(True)
#
## do not show ticks on 1st graph
#ax.tick_params(
#    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='on',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='on', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
## do not show ticks
#ax.tick_params(
#    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='off',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='off', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
#print('longest chromo', longest, maxlength)
#
## determine tick position on x axis
#xpos =  [j for j in range(0, len(Positions[longest]), 10)]
## convert interval windows numbers to genomic positions
#xtext = list(map(lambda x : (x * Interval) / 1000000, xpos))
#Xtext = []
#for i in xtext:
#    if i % 2 == 0:
#        Xtext.append(str(int(i)))
#    else:
#        Xtext.append('')
#
## set up tick positions and labels
#plt.xticks(xpos, Xtext, rotation = 0, fontsize = 10, fontname = 'Arial')
#
## add lines
#lns = Graph[LG[0]]
#for chromo in LG[1:]:
#    lns += Graph[chromo]
## get labels
#labs = []
#for chromo in LG:
#    assert chromo.count('_') == 2
#    lg = chromo[chromo.index('_', chromo.index('_')+1, -1)+1:]
#    labs.append('LG' + lg)
## plot legend
#ax.legend(lns, labs, loc=1, fontsize = 8, frameon = False)
#
#fig.savefig('truc.pdf', bbox_inches = 'tight')