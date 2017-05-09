# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 22:27:34 2017

@author: Richard
"""

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



# use this scipt to test for enrichement of disease genes among external genes

# load dictionaries of host and nested genes 
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanNestedGenes.json') as human_json_data:
    Nested = json.load(human_json_data)

# get GFF file
GFF = 'Homo_sapiens.GRCh38.88.gff3'
 
# find nested and intronic-nested genes 
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
# Map Transcript names to gene names {transcript: gene}
MapTranscriptGene = TranscriptToGene(GFF)
# get the coordinates of all exons    
ExonCoord = GeneExonCoord(GFF)
ExonCoord = CleanGeneFeatureCoord(ExonCoord, MapTranscriptGene)
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
IntronCoord = GeneIntronCoord(ExonCoord)
IntronCoord = CleanGeneFeatureCoord(IntronCoord, MapTranscriptGene)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
TranscriptCoordinates = TranscriptsCoord(GFF)
# map genes to their longest transcript {gene: longest_transcript}
GeneLongestTranscript = LongestTranscript(TranscriptCoordinates, MapGeneTranscript)
# match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
Matches = MatchHostTranscriptWithNestedTranscript(Nested, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# list all host, nested transcript pairs [[host, nested]]
HostNestedPairs = GetHostNestedPairs(Matches)

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


# compare proportions of disease and non-disease genes among external genes
# for which with internal genes that are intronless or intron-containing

# compare proportions of disease and non-disease genes among internal genes
# that are intronless or intron-containing

# make sets of external genes with intronless and intron-containing internal genes
# make sets of internal genes with intron or without intron
ExtWithIntrons, ExtNoIntrons, IntWithIntrons, IntNoIntrons = set(), set(), set(), set()
 
# loop over host-nested transcript pairs
for i in range(len(HostNestedPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
    if HostNestedPairs[i][1] in IntronCoord:
        # internal gene has introns, populate set with external gene name
            ExtWithIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
            IntWithIntrons.add(MapTranscriptGene[HostNestedPairs[i][1]])
    else:
        # internal gene is intronless, populate set with external gene name
        ExtNoIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
        IntNoIntrons.add(MapTranscriptGene[HostNestedPairs[i][1]])
    
# make a list of external genes
ExtGenes = [ExtWithIntrons, ExtNoIntrons]  
# make a list of external genes
IntGenes = [IntWithIntrons, IntNoIntrons]  


# count disease and non-disease genes in external and internal genes [[disease, non_disease],... ]
ExtCounts = []
for DiseaseOrigin in [GAD, GWAS, Drivers, OMIM, DiseaseGenes]:
    Counts = []
    for i in range(len(ExtGenes)):
        disease = len([j for j in ExtGenes[i] if j in DiseaseOrigin])
        nondisease = len([j for j in ExtGenes[i] if j not in DiseaseOrigin])
        Counts.append([disease, nondisease])
    ExtCounts.append(Counts)
IntCounts = []
for DiseaseOrigin in [GAD, GWAS, Drivers, OMIM, DiseaseGenes]:
    Counts = []
    for i in range(len(IntGenes)):
        disease = len([j for j in IntGenes[i] if j in DiseaseOrigin])
        nondisease = len([j for j in IntGenes[i] if j not in DiseaseOrigin])
        Counts.append([disease, nondisease])
    IntCounts.append(Counts)

# compare the proportion of disease genes for each set of disease gene
ExtPVals = []
for i in range(len(ExtCounts)):
    p = stats.fisher_exact([ExtCounts[i][0], ExtCounts[i][1]])[1]
    ExtPVals.append(p)    
IntPVals = []
for i in range(len(IntCounts)):
    p = stats.fisher_exact([IntCounts[i][0], IntCounts[i][1]])[1]
    IntPVals.append(p)

# convert P values to significance
ExtPVals = ConvertPToStars(ExtPVals)
IntPVals = ConvertPToStars(IntPVals)

# compute proportions of disease and non-disease
ExtDisProp, ExtNonDisProp = [], []
for i in range(len(ExtCounts)):
    disease, nondisease = [], []
    for j in range(len(ExtCounts[i])):
        disease.append(ExtCounts[i][j][0] / sum(ExtCounts[i][j]))
        nondisease.append(ExtCounts[i][j][1] / sum(ExtCounts[i][j]))
    ExtDisProp.append(disease)    
    ExtNonDisProp.append(nondisease)
IntDisProp, IntNonDisProp = [], []
for i in range(len(IntCounts)):
    disease, nondisease = [], []
    for j in range(len(IntCounts[i])):
        disease.append(IntCounts[i][j][0] / sum(IntCounts[i][j]))
        nondisease.append(IntCounts[i][j][1] / sum(IntCounts[i][j]))
    IntDisProp.append(disease)
    IntNonDisProp.append(nondisease)
 
# create a list of lists of external and internal disease genes with and without introns for each disease class
DisProp, NonDisProp = [], []
for i in range(len(ExtDisProp)):
    ExtIntGenes = []
    ExtIntGenes.extend(ExtDisProp[i])
    ExtIntGenes.extend(IntDisProp[i])
    DisProp.append(ExtIntGenes)
for i in range(len(ExtNonDisProp)):
    ExtIntGenes = []
    ExtIntGenes.extend(ExtNonDisProp[i])
    ExtIntGenes.extend(IntNonDisProp[i])
    NonDisProp.append(ExtIntGenes)
    
    
# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YRange, isYLabel):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot proportions of disease genes only
    # Create a horizontal bar plot for proportions of disease genes
    ax.bar([0, 0.2, 0.5, 0.7], Data[0], width = 0.2, label = 'disease', color= ['#2b8cbe', '#fd8d3c', '#2b8cbe', '#fd8d3c'], edgecolor = 'white', linewidth = 0.7)
    ax.bar([0, 0.2, 0.5, 0.7], Data[1], width = 0.2, bottom = Data[0], label = 'non-disease', color= ['#88419d', '#e31a1c', '#88419d', '#e31a1c'], edgecolor = 'white', linewidth = 0.7)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write axis labels    
    if isYLabel == True:
        YLabel = 'Proportion of disease genes'
        ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    plt.xticks([0.2, 0.7], ['ext', 'int'], color = 'black',  size = 7, ha = 'center', **FigFont)
    # add x label    
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # edit y axis ticks
    if isYLabel == True:
        plt.yticks(YRange)    
    # add a range for the Y axis
    plt.ylim([0, 1])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    if isYLabel == True:
        ax.spines["left"].set_visible(True)
    else:
        ax.spines["left"].set_visible(False)
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 3))        
    # edit tick parameters    
    if isYLabel == True:
        plt.tick_params(axis='both', which='both', bottom='on', top='off',
                        right = 'off', left = 'on', labelbottom='on', colors = 'black',
                        labelsize = 7, direction = 'out')  
    else:
        plt.tick_params(axis='both', which='both', bottom='on', top='off',
                        right = 'off', left = 'off', labelbottom='on', labelleft = 'off', colors = 'black',
                        labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # add margins
    plt.margins(0.1)
    return ax


# make a figure with proportion of disease and non-disease genes

# create figure
fig = plt.figure(1, figsize = (6, 2))
# plot data
ax1 = CreateAx(5, 1, 1, fig, [DisProp[0], NonDisProp[0]], 'complex', np.arange(0, 1.25, 0.25), True)
ax2 = CreateAx(5, 1, 2, fig, [DisProp[1], NonDisProp[1]], 'GWAS', np.arange(0, 1.25, 0.25), False)
ax3 = CreateAx(5, 1, 3, fig, [DisProp[2], NonDisProp[2]], 'tumors', np.arange(0, 1.25, 0.25), False)
ax4 = CreateAx(5, 1, 4, fig, [DisProp[3], NonDisProp[3]], 'mendelian', np.arange(0, 1.25, 0.25), False)
ax5 = CreateAx(5, 1, 5, fig, [DisProp[4], NonDisProp[4]], 'all', np.arange(0, 1.25, 0.25), False)

# annotate figure to add significance
for i in range(len(ExtPVals)):
    if ExtPVals[i]  != '':
        if i == 0:
            ax1 = AddSignificanceToBars(ax1, ExtPVals[i], 0.1, 0.3, 1, 0.2, 1.05)
        elif i == 1:
            ax2 = AddSignificanceToBars(ax2, ExtPVals[i], 0.1, 0.3, 1, 0.2, 1.05)
        elif i == 2:
            ax3 = AddSignificanceToBars(ax3, ExtPVals[i], 0.1, 0.3, 1, 0.2, 1.05)
        elif i == 3:
            ax4 = AddSignificanceToBars(ax4, ExtPVals[i], 0.1, 0.3, 1, 0.2, 1.05)
        elif i == 4:
            ax5 = AddSignificanceToBars(ax5, ExtPVals[i], 0.1, 0.3, 1, 0.2, 1.05)
for i in range(len(IntPVals)):
    if IntPVals[i]  != '':
        if i == 0:
            ax1 = AddSignificanceToBars(ax1, IntPVals[i], 0.6, 0.8, 1, 0.7, 1.05)
        elif i == 1:
            ax2 = AddSignificanceToBars(ax2, IntPVals[i], 0.6, 0.8, 1, 0.7, 1.05)
        elif i == 2:
            ax3 = AddSignificanceToBars(ax3, IntPVals[i], 0.6, 0.8, 1, 0.7, 1.05)
        elif i == 3:
            ax4 = AddSignificanceToBars(ax4, IntPVals[i], 0.6, 0.8, 1, 0.7, 1.05)
        elif i == 4:
            ax5 = AddSignificanceToBars(ax5, IntPVals[i], 0.6, 0.8, 1, 0.7, 1.05)

# add subplot labels
ax1.text(-0.45, 1.1, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax2.text(0, 1.1, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax3.text(0, 1.1, 'C', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax4.text(0, 1.1, 'D', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax5.text(0, 1.1, 'E', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)

# add legend
DI = mpatches.Patch(facecolor = '#2b8cbe', edgecolor = 'white', linewidth = 0.7, label= 'disease & intron')
DN = mpatches.Patch(facecolor = '#fd8d3c', edgecolor = 'white', linewidth = 0.7, label= 'disease & intronless')
NDI = mpatches.Patch(facecolor = '#88419d', edgecolor = 'white', linewidth = 0.7, label= 'non-disease & intron')
NDN = mpatches.Patch(facecolor = '#e31a1c' , edgecolor = 'white', linewidth = 0.7, label= 'non-disease & intronless')
ax1.legend(handles = [DI, DN, NDI, NDN], loc = (0.05, 1.2), fontsize = 6, frameon = False, ncol = 4)

# make sure subplots do not overlap
plt.tight_layout()

outputfile = 'ProportionDiseaseNested'
# save figure to file
fig.savefig(outputfile + '.pdf', bbox_inches = 'tight')
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
