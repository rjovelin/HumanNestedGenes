# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:54:24 2017

@author: RJovelin
"""

# use this script to test for enrichement of overlapping genes among disease genes
# save data as table and figure

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
GeneSets = []
for i in range(len(Overlap)):
    GeneSets.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

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

# create a list of gene sets
AllGenes = [NonOverlappingGenes]
for i in range(1, len(GeneSets)):
    AllGenes.append(GeneSets[i])
# create a parallel list of gene categories
GeneCats = ['NoOvl', 'Nst', 'Pgk', 'Con', 'Div'] 

# make a list of counts for each disease and gene category[[[disease, non_disease],... ], ....]
GeneCounts = CountDiseaseGenes(AllGenes, [GAD, GWAS, Drivers, OMIM, DiseaseGenes])

# test significance by comparing each overlapping gene class to non-overlapping genes for each disease
Significance = []
for i in range(len(GeneCounts)):
    # access disease type
    PVals = []
    for j in range(1, len(GeneCounts[i])):
        # compare each gene class to non-overlapping genes for that disease class
        P = stats.fisher_exact([GeneCounts[i][0], GeneCounts[i][j]])[1]
        PVals.append(P)
    # convert p values to star significance
    Significance.append(ConvertPToStars(PVals))


# compute proportions of disease and non-disease for each class in each disease class
DisProp, NonDisProp = [], []
for i in range(len(GeneCounts)):
    # access disease type
    disease, nondisease = [], []
    for j in range(len(GeneCounts[i])):
        # access disease and non-disease counts for given gene class in this disease class
        disease.append(GeneCounts[i][j][0] / sum(GeneCounts[i][j]))
        nondisease.append(GeneCounts[i][j][1] / sum(GeneCounts[i][j]))
    DisProp.append(disease)    
    NonDisProp.append(nondisease)


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, Title, YRange, XLabel):
    '''
    Returns a ax instance in figure
    '''    
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # check if plot only disease genes or proportions of disease and non-disease genes
    if Proportions == 'both':
        # Create a horizontal bar plot for proportions of disease genes
        ax.bar([0, 0.3, 0.6, 0.9, 1.2], Data[0], width = 0.2, label = 'disease', color= 'black', linewidth = 0.7)
        # Create a horizontal bar plot for proportions of non-disease genes
        ax.bar([0, 0.3, 0.6, 0.9, 1.2], Data[1], width = 0.2, bottom = Data[0], label = 'non-disease', color= 'lightgrey', linewidth = 0.7)
    

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    if Proportions == 'both':
        YLabel = 'Proportion'
    elif Proportions == 'disease':
        YLabel = 'Proportion of\ndisease genes'
    elif Proportions == 'non-disease':    
        YLabel = 'Proportion of\nnon-disease genes'   
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
        
    # add ticks and lebels
    if XLabel == True:
        plt.xticks([0.1, 0.4, 0.7, 1, 1.3], ['NoOvl', 'Nst', 'Int', 'Ext', 'Pgk', 'Con', 'Div'], rotation = 30, size = 7, color = 'black', ha = 'right', **FigFont)
    elif XLabel == False:
        plt.xticks([0.1, 0.4, 0.7, 1, 1.3], [''] * 7, size = 7, color = 'black', ha = 'right', **FigFont)
    
    # edit y axis ticks
    plt.yticks(YRange)    
        
    # add title
    plt.title(Title, color = 'black',  size = 7, ha = 'center', **FigFont)
    # add a range for the Y axis
    plt.ylim([0, YMax])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
    # edit tick parameters    
    if XLabel == True:
        plt.tick_params(axis='both', which='both', bottom='on', top='off',
                        right = 'off', left = 'on', labelbottom='on',
                        colors = 'black', labelsize = 7, direction = 'out')  
    else:
        plt.tick_params(axis='both', which='both', bottom='on', top='off',
                        right = 'off', left = 'on', labelbottom='off',
                        colors = 'black', labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # add margins
    plt.margins(0.1)
    return ax


# make a figure with proportion of disease and non-disease genes

# create figure
fig = plt.figure(1, figsize = (2.5, 6))
# plot data
ax1 = CreateAx(1, 5, 1, fig, [GADDis, GADNonDis], 'complex diseases', 'disease', np.arange(0, 0.71, 0.1), 0.71, False)
ax2 = CreateAx(1, 5, 2, fig, [GWASDis, GWASNonDis], 'GWAS', 'disease', np.arange(0, 0.21, 0.05), 0.2, False)
ax3 = CreateAx(1, 5, 3, fig, [DriversDis, DriversNonDis], 'tumor drivers', 'disease', np.arange(0, 0.041, 0.010), 0.04,  False)
ax4 = CreateAx(1, 5, 4, fig, [OMIMDis, OMIMNonDis], 'medelian diseases', 'disease', np.arange(0, 0.26, 0.05), 0.25, False)
ax5 = CreateAx(1, 5, 5, fig, [AllDis, AllNonDis], 'all diseases', 'disease', np.arange(0, 0.71, 0.1), 0.71, True)






# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]

ypos = [0.55, 0.50, 0.65, 0.55, 0.6, 0.6]
for i in range(len(PValGAD)):
    ax1.text(xpos[i], ypos[i], PValGAD[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
ypos = [0.12, 0.05, 0.16, 0.09, 0.10, 0.10]
for i in range(len(PValGWAS)):
    ax2.text(xpos[i], ypos[i], PValGWAS[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
ypos = [0.027, 0.017, 0.040, 0.015, 0.037, 0.030]
for i in range(len(PValDrivers)):
    ax3.text(xpos[i], ypos[i], PValDrivers[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
ypos = [0.17, 0.12, 0.25, 0.17, 0.22, 0.22]
for i in range(len(PValOMIM)):
    ax4.text(xpos[i], ypos[i], PValOMIM[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
ypos = [0.55, 0.45, 0.7, 0.55, 0.65, 0.65]
for i in range(len(PValAll)):
    ax5.text(xpos[i], ypos[i], PValAll[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)

Proportions = 'disease'
if Proportions == 'both':
    # add legend
    N = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'non-disease')
    D = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'disease')
    ax1.legend(handles = [D, N], loc = (0, 1.1), fontsize = 6, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure to file
fig.savefig('truc.pdf', bbox_inches = 'tight')


#fig.savefig('ProportionDiseaseGenes.pdf', bbox_inches = 'tight')
#fig.savefig('ProportionDiseaseGenes.eps', bbox_inches = 'tight')


