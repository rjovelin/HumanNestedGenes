# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:14:49 2016

@author: RJovelin
"""

# use this script to plot the proportion of overlapping and non-overlapping genes
# with highest expression in each tissue

# usage python3 PlotHighestExpression.py 


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


# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)
NestedPairs = GetHostNestedPairs(Nested)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)

# generate gene sets
NestedGenes  = MakeFullPartialOverlapGeneSet(Nested)
OverlappingGenes = MakeFullPartialOverlapGeneSet(Overlapping)
ConvergentGenes = MakeFullPartialOverlapGeneSet(Convergent)
DivergentGenes = MakeFullPartialOverlapGeneSet(Divergent)
PiggyBackGenes = MakeFullPartialOverlapGeneSet(Piggyback)

# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)

# create lists of nested gene pairs with same and opposite directions
same, opposite = [], []
for pair in NestedPairs:
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        opposite.append(pair)
    elif len(set(orientation)) == 1:
        same.append(pair)
# create sets of internal and external nested genes depending on orientation 
InternalSameGenes, InternalOppositeGenes, ExternalSameGenes, ExternalOppositeGenes = set(), set(), set(), set()
for pair in same:
    ExternalSameGenes.add(pair[0])
    InternalSameGenes.add(pair[1])
for pair in opposite:
    ExternalOppositeGenes.add(pair[0])
    InternalOppositeGenes.add(pair[1])


# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseGTEXExpression('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)


# make a list of all gene sets
AllGeneSets = [NonOverlappingGenes, InternalSameGenes, InternalOppositeGenes, ExternalSameGenes,
               ExternalOppositeGenes, PiggyBackGenes, ConvergentGenes, DivergentGenes]

# make a list of tissues
infile = open('GTEX_Median_Normalized_FPKM.txt')
header = infile.readline().rstrip().split('\t')
Tissues = header[1:]
# replace spaces in tissue names
for i in range(len(Tissues)):
    if ' ' in Tissues[i]:
        Tissues[i] = Tissues[i].replace(' ', '_')


# make a list of genes with expression
Expressed = list(ExpressionProfile.keys())
# check that all genes have the same number of tissues
for gene in Expressed:
    assert len(ExpressionProfile[gene]) == len(Tissues)


# make a list of gene category names parallel to the list of gene sets
GeneCats = ['NoOv', 'CisInt', 'TransInt', 'CisExt', 'TransExt', 'Pbk', 'Conv', 'Div']


# create a dict to count the number of genes in each category with highest expression in each tissue
HighestExpression = {}
# inititialize dict with list of 0s
for i in GeneCats:
    HighestExpression[i] = [0] * len(Tissues)










######

#0       Gene    ENSG00000139133
#1       Adipose Tissue  3490.06030824
#2       Adrenal Gland   4311.17097459
#3       Bladder 3067.69288609
#4       Blood   1240.69549277
#5       Blood Vessel    2331.66199801
#6       Brain   1122.53401727
#7       Breast  3635.65341453
#8       Cervix Uteri    2848.49316676
#9       Colon   2148.5069133
#10      Esophagus       2407.01134852
#11      Fallopian Tube  2721.57143351
#12      Heart   1680.920066
#13      Kidney  2205.3354964
#14      Liver   1522.58459886
#15      Lung    3151.35898168
#16      Muscle  1422.70689811
#17      Nerve   3465.28093996
#18      Ovary   4142.44386425
#19      Pancreas        2988.7953959
#20      Pituitary       2642.30711158
#21      Prostate        2307.78470022
#22      Salivary Gland  4183.08771799
#23      Skin    3336.5844231
#24      Small Intestine 3058.60868661
#25      Spleen  3017.8843581
#26      Stomach 2306.89539816
#27      Testis  6008.86105018
#28      Thyroid 3603.82870762
#29      Uterus  3827.06925676
#30      Vagina  2828.22684982


######



    # create lists to count the number of genes with highest expression in each tissues
    # [brain, cerebellum, heart, kidney, liver, testis]   
    hosthigh, nestedhigh, controlhigh = [0]*6, [0]*6, [0]*6
    # loop over the control genes, determine the tissue with highest expression, update counter
    for control in SpControl:
        # get the index of the maximum expression value
        pos = SpExpression[control].index(max(SpExpression[control]))
        # check that there is a single maximum expression value
        assert SpExpression[control].count(max(SpExpression[control])) == 1, 'control has > 1 max value'
        # update counter list at pos index
        controlhigh[pos] += 1
    # loop over gene pairs:
    for pair in SpHostNestedPairs:
        # get the index of the maximum expression value for the host and nested genes
        poshost = SpExpression[pair[0]].index(max(SpExpression[pair[0]]))         
        posnested = SpExpression[pair[1]].index(max(SpExpression[pair[1]]))
        # check that there is a single maximum expression value
        assert SpExpression[pair[0]].count(max(SpExpression[pair[0]])) == 1, 'host has > 1 max value'        
        assert SpExpression[pair[1]].count(max(SpExpression[pair[1]])) == 1, 'nested has > 1 max value'
        # update counter list as position index
        hosthigh[poshost] += 1
        nestedhigh[posnested] += 1
    # divide counts by total of genes in each category to get the proportions
    for m in range(len(hosthigh)):
        hosthigh[m] = hosthigh[m] / len(SpHostNestedPairs)
    for m in range(len(nestedhigh)):
        nestedhigh[m] = nestedhigh[m] / len(SpHostNestedPairs)
    for m in range(len(controlhigh)):
        controlhigh[m] = controlhigh[m] / len(SpControl)
    # populate lists
    HostHighest.append(hosthigh)
    NestedHighest.append(nestedhigh)    
    ControlHighest.append(controlhigh)
    
# make lists with the proportions of host, nested and control; genes in each tissues
# [brain, cerebellum, heart, kidney, liver, testis]   
HumanExp = []
for i in range(len(HostHighest[0])):
    HumanExp.append(HostHighest[0][i])
    HumanExp.append(NestedHighest[0][i])
    HumanExp.append(ControlHighest[0][i])
ChimpExp = []
for i in range(len(HostHighest[1])):
    ChimpExp.append(HostHighest[1][i])
    ChimpExp.append(NestedHighest[1][i])
    ChimpExp.append(ControlHighest[1][i])
GorillaExp = []
for i in range(len(HostHighest[2])):
    GorillaExp.append(HostHighest[2][i])
    GorillaExp.append(NestedHighest[2][i])
    GorillaExp.append(ControlHighest[2][i])
OrangutanExp = []
for i in range(len(HostHighest[3])):
    OrangutanExp.append(HostHighest[3][i])
    OrangutanExp.append(NestedHighest[3][i])
    OrangutanExp.append(ControlHighest[3][i])
MacaqueExp = []
for i in range(len(HostHighest[4])):
    MacaqueExp.append(HostHighest[4][i])
    MacaqueExp.append(NestedHighest[4][i])
    MacaqueExp.append(ControlHighest[4][i])

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Expression, XLabel, YLabel, YMax):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['#a6cee3','#1f78b4','#b2df8a']
    # plot nucleotide divergence
    ax.bar([0, 0.2, 0.4,
            0.7, 0.9, 1.1,
            1.4, 1.6, 1.8,
            2.1, 2.3, 2.5,
            2.8, 3, 3.2,
            3.5, 3.7, 3.9], Expression, 0.2, color = colorscheme,
           edgecolor = 'black', linewidth = 1,
           error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    
    # write x ticks
    plt.xticks([0.3, 1, 1.7, 2.4, 3.1, 3.8], ['br', 'cb', 'ht', 'kd', 'lv', 'ts'], ha = 'center', fontsize = 8, **FigFont)
    
    # add a range for the Y axis
    plt.ylim([0, YMax])
       
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(True)  
        
    # do not show ticks
    plt.tick_params(
        axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        right = 'off',
        left = 'on',          
        labelbottom='on', # labels along the bottom edge are on
        colors = 'black',
        labelsize = 8,
        direction = 'out') # ticks are outside the frame when bottom = 'on'  
         
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
     
    # create a margin around the x axis
    #plt.margins(0.1)
    return ax      


# create figure
fig = plt.figure(1, figsize = (3.5, 7.5))

# plot data
ax1 = CreateAx(1, 5, 1, fig, HumanExp, 'Human', 'Proportion of genes', 0.41)
ax2 = CreateAx(1, 5, 2, fig, ChimpExp, 'Chimp', 'Proportion of genes', 0.41)
ax3 = CreateAx(1, 5, 3, fig, GorillaExp, 'Gorilla', 'Proportion of genes', 0.41)
ax4 = CreateAx(1, 5, 4, fig, OrangutanExp, 'Orangutan', 'Proportion of genes', 0.41)
ax5 = CreateAx(1, 5, 5, fig, MacaqueExp, 'Macaque', 'Proportion of genes', 0.41)

# add legend relative to ax1 using ax1 coordinates
H = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'Hosts')
N = mpatches.Patch(facecolor = '#1f78b4', edgecolor = 'black', linewidth = 1, label= 'Nested')
U = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'Control')
ax1.legend(handles = [H, N, U], loc = (0, 1), fontsize = 8, frameon = False, ncol = 3)

# make sure subplots do not overlap
plt.tight_layout()

fig.savefig('truc.pdf', bbox_inches = 'tight')
#fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
