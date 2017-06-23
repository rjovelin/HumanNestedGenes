# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 11:00:15 2017

@author: RJovelin
"""

# use this script to plot expression divergence between overlapping genes for disease and non-disease genes


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
import string
import numpy as np
from scipy import stats
from HsaNestedGenes import *


# load dictionaries
Overlap = []
for filename in ['PiggyBack', 'Convergent', 'Divergent']:
    json_data = open('Human' + filename + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    Overlap.append(overlapping)

# make a list of lists of overlapping gene pairs
OverlapPairs = []
for overlapping in Overlap:
    OverlapPairs.append(GetHostNestedPairs(overlapping))

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
print('parsed diseases')

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)
print('extracted expression')

# make a list of lists with gene pairs for each overlapping group and disease
AllPairs = []
for i in range(len(OverlapPairs)):
    DiseasePairs = []    
    # loop over disease
    for disease in [GAD, GWAS, Drivers, OMIM, DiseaseGenes]:
        # create list with at least 1 disease gene and none being disease genes
        atleast, none = [], []
        # loop over gene pairs
        for pair in OverlapPairs[i]:
            # check that both genes are expressed
            if pair[0] in ExpressionProfile and pair[1] in ExpressionProfile:
                # check which genes are disease genes
                if pair[0] in disease or pair[1] in disease:
                    atleast.append(pair)
                elif pair[0] not in disease and pair[1] not in disease:
                    none.append(pair)
        # add list to list for the given disease
        DiseasePairs.append([atleast, none])
    # add list with all gene pairs for each disease for the given overlapping group
    AllPairs.append(DiseasePairs)

# make a list with disease types
DiseaseGroups = ['complex', 'GWAS', 'cancer', 'mendelian', 'all']
# make a list of disease groups for each overlapping group
DiseaseClasses = DiseaseClasses = [copy.deepcopy(DiseaseGroups) for i in range(3)]

# filter diseases if any pair count < 20
# loop over overlapping groups
for i in range(len(AllPairs)):
    DiseaseToRemove, PairsToRemove = [], []
    # loop over disease in current overlapping group
    for j in range(len(AllPairs[i])):
        # check each pair count for current disease
        for k in range(len(AllPairs[i][j])):
            if len(AllPairs[i][j][k]) < 20:
                if AllPairs[i][j] not in PairsToRemove:
                    PairsToRemove.append(AllPairs[i][j])
                if DiseaseClasses[i][j] not in DiseaseToRemove:
                    DiseaseToRemove.append(DiseaseClasses[i][j])
    # remove lists and diseases      
    for item in DiseaseToRemove:
        DiseaseClasses[i].remove(item)
    for L in PairsToRemove:
        AllPairs[i].remove(L)
print(DiseaseClasses)

# do QC
for i in range(len(AllPairs)):
    for j in range(len(AllPairs[i])):
        for k in range(len(AllPairs[i][j])):
            assert len(AllPairs[i][j][k]) >= 20
   
# compute expression divergence between pairs of genes
Divergence = []
# loop over overlapping groups
for i in range(len(AllPairs)):
    Div = []
    # loop over disease groups
    for j in range(len(AllPairs[i])):
        Ddisease = ComputeExpressionDivergenceGenePairs(AllPairs[i][j][0], ExpressionProfile)
        Dnondisease = ComputeExpressionDivergenceGenePairs(AllPairs[i][j][1], ExpressionProfile)        
        Div.append([Ddisease, Dnondisease])
    Divergence.append(Div)

# create lists with means and SEM for each gene category
Means, SEMs = [], []
# loop over overlapping groups
for i in range(len(Divergence)):
    # create lists with means and sem
    average, error = [], []
    # loop over diseases
    for j in range(len(Divergence[i])):
        average.append(np.mean(Divergence[i][j][0]))
        error.append(np.std(Divergence[i][j][0]) / math.sqrt(len(Divergence[i][j][0])))
        average.append(np.mean(Divergence[i][j][1]))
        error.append(np.std(Divergence[i][j][1]) / math.sqrt(len(Divergence[i][j][1])))
        assert len(Divergence[i][j]) == 2
    Means.append(average)
    SEMs.append(error)        

# perform statistical tests between gene categories
PValues = []
# loop over overlapping groups
for i in range(len(Divergence)):
    pvals = []
    # loop over disease in current overlapping group
    for j in range(len(Divergence[i])):
        P = PermutationResampling(Divergence[i][j][0], Divergence[i][j][1], 1000, statistic = np.mean)
        pvals.append(P)
    PValues.append(pvals)
# convert P to stars
for i in range(len(PValues)):
    PValues[i] = ConvertPToStars(PValues[i])
   
# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XTickLabel, YRange, Title):
    '''
    return an ax object part of figure
    '''

    # add a plot to figure (N row, N column, plot N)
    ax = fig.add_subplot(Rows, Columns, Position)
    # set colors
    colorscheme = ['#225ea8', '#e31a1c']
    # plot nucleotide divergence
    if len(XTickLabel) == 4:
        ax.bar([0.05, 0.25, 0.55, 0.75, 1.05, 1.25, 1.55, 1.75], Data[0], 0.2, yerr = Data[1], color = colorscheme,
               edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    elif len(XTickLabel) == 5:
        ax.bar([0.05, 0.25, 0.55, 0.75, 1.05, 1.25, 1.55, 1.75, 2.05, 2.25], Data[0], 0.2, yerr = Data[1], color = colorscheme,
               edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
    # add ticks and lebels
    if len(XTickLabel) == 4:
        plt.xticks([0.25, 0.75, 1.25, 1.75], XTickLabel, size = 7, color = 'black', ha = 'center', **FigFont)
    elif len(XTickLabel) == 5:
        plt.xticks([0.25, 0.75, 1.25, 1.75, 2.25], XTickLabel, size = 7, color = 'black', ha = 'center', **FigFont)
    # add title
    ax.set_title(Title, color = 'black', size = 7, ha = 'center', **FigFont)    
    # add a range for the Y and X axes
    plt.yticks(YRange)    
    if len(XTickLabel) == 4:
        plt.xlim([0, 2])
    elif len(XTickLabel) == 5:
        plt.xlim([0, 2.5])
    plt.ylim([0, 0.4])
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
    return ax  


Titles = ['Piggyback', 'Convergent', 'Divergent']
# create figure
fig = plt.figure(1, figsize = (3, 4.5))
# plot data
for i in range(len(Means)):
    # get positions to annotate p values
    XPos = np.array([0.15, 0.35, 0.25])
    ax = CreateAx(1, 3, i+1, fig, [Means[i], SEMs[i]], DiseaseClasses[i], np.arange(0, 0.6, 0.2), Titles[i])
    # add significance
    for j in range(len(PValues[i])):
        # update positions of p values (add 0.50 to all positions)
        if j != 0:
            XPos += 0.50
        if PValues[i][j] != '':
            # add stars for significance
            ax = AddSignificanceToBars(ax, PValues[i][j], XPos[0], XPos[1], 0.36, XPos[2], 0.39)
    # add legend
    if i == 0:
        D = mpatches.Patch(facecolor = '#225ea8', edgecolor = 'black', linewidth = 0.7, label= 'disease')
        N = mpatches.Patch(facecolor = '#e31a1c', edgecolor = 'black', linewidth = 0.7, label= 'non-disease')
        ax.legend(handles = [D, N], loc = (0.05, 1.15), fontsize = 6, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()    
    
# save figure
outputfile = 'ExpDivOvlpDisease'
for extension in ['.pdf', '.eps', '.png']:
    fig.savefig(outputfile + extension, bbox_inches = 'tight')
