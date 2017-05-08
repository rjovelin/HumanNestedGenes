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
Drivers = ParseTumorDrivers('driver_genes_per_tumor_syn7314119.csv')
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




# use this function to get gene proportions
def GetProportions(Counts):
    '''
    (list) -> list, list
    Take the list of inner lists with counts of disease and non-disease
    and return 2 lists with proportions of disease and non-disease genes respectively
    '''
    disease, nondisease = [], []
    for i in range(len(Counts)):
        disease.append(Counts[i][0] / sum(Counts[i]))
        nondisease.append(Counts[i][1] / sum(Counts[i]))
    return disease, nondisease




    
 
    
# create lists of proportions for disease and non-disease genes
DisProp, NonDisProp = [], []
for i in range(len(counts)):
    disease, nondisease = GetProportions(counts[i])
    DisProp.append(disease)
    NonDisProp.append(nondisease)


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YRange, YMax):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot proportions of disease genes only
    # Create a horizontal bar plot for proportions of disease genes
    ax.bar([0, 0.3], Data, width = 0.2, label = 'disease', color= ['black'] + ['lightgrey'] * 6, linewidth = 0.7)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write axis labels    
    YLabel = 'Proportion of disease genes'
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # edit y axis ticks
    plt.yticks(YRange)    
    # add a range for the Y axis
    plt.ylim([0, YMax])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right = 'off', left = 'on', labelbottom='off', colors = 'black',
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
ax1 = CreateAx(5, 1, 1, fig, DisProp[0], 'complex diseases', np.arange(0, 1, 0.1), 0.8)
ax2 = CreateAx(5, 1, 2, fig, DisProp[1], 'GWAS', np.arange(0, 0.30, 0.05), 0.25)
ax3 = CreateAx(5, 1, 3, fig, DisProp[2], 'tumor drivers', np.arange(0, 0.1, 0.02), 0.08)
ax4 = CreateAx(5, 1, 4, fig, DisProp[3], 'medelian diseases', np.arange(0, 0.40, 0.05), 0.351)
ax5 = CreateAx(5, 1, 5, fig, DisProp[4], 'all diseases', np.arange(0, 1, 0.1), 0.8)

# use this function to annotate the graph with significance levels
def AddSignificance(ax, SignificanceLevel, XLine1, XLine2, YLine, XText, YText):
    '''
    (ax, str, num, num, num, num, num) -> ax
    Take a matplotlib ax object, the significance level (as stars), the positions
    of the bracket and star and return the ax with annotated significance level
    '''
    ax.annotate("", xy=(XLine1, YLine), xycoords='data', xytext=(XLine2, YLine), textcoords='data',
                 arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 0.7))
    # add stars for significance
    ax.text(XText, YText, SignificanceLevel, horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 7)
    return ax


# annotate figure to add significance
ax1 = AddSignificance(ax1, PVals[0], 0.1, 0.4, 0.73, 0.25, 0.78)
ax2 = AddSignificance(ax2, PVals[1], 0.1, 0.4, 0.22, 0.25, 0.24)
ax3 = AddSignificance(ax3, PVals[2], 0.1, 0.4, 0.07, 0.25, 0.075)
ax4 = AddSignificance(ax4, PVals[3], 0.1, 0.4, 0.32, 0.25, 0.34)
ax5 = AddSignificance(ax5, PVals[4], 0.1, 0.4, 0.75, 0.25, 0.79)


# add subplot labels
ax1.text(-0.45, 0.85, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(0.8, 0.85, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(2.2, 0.85, 'C', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(3.6, 0.85, 'D', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(5, 0.85, 'E', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)


# add legend
N = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'intronless internal genes')
D = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'internal genes with introns')
ax1.legend(handles = [D, N], loc = (1, 1.15), fontsize = 6, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

outputfile = 'ProportionDiseaseExternal'
# save figure to file
fig.savefig('truc.pdf', bbox_inches = 'tight')

############################## $$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# use this function to get gene proportions
def GetProportions(Counts):
    '''
    (list) -> list, list
    Take the list of inner lists with counts of disease and non-disease
    and return 2 lists with proportions of disease and non-disease genes respectively
    '''
    disease, nondisease = [], []
    for i in range(len(Counts)):
        disease.append(Counts[i][0] / sum(Counts[i]))
        nondisease.append(Counts[i][1] / sum(Counts[i]))
    return disease, nondisease
    
# count disease and non-disease genes    
GADCounts = CountDiseaseGenes(IntGenes, GAD)    
GWASCounts = CountDiseaseGenes(IntGenes, GWAS)
DriversCounts = CountDiseaseGenes(IntGenes, Drivers)
OMIMCounts = CountDiseaseGenes(IntGenes, OMIM)
AllCounts = CountDiseaseGenes(IntGenes, DiseaseGenes)

# compare the proportion of disease genes for each set of disease gene
PVals = []
counts = [GADCounts, GWASCounts, DriversCounts, OMIMCounts, AllCounts]
for i in range(len(counts)):
    p = stats.fisher_exact([counts[i][0], counts[i][1]])[1]
    PVals.append(p)    
    
# transform p values in stars
PVals = AssignSignificance(PVals)    
    
# create lists of proportions for disease and non-disease genes
DisProp, NonDisProp = [], []
for i in range(len(counts)):
    disease, nondisease = GetProportions(counts[i])
    DisProp.append(disease)
    NonDisProp.append(nondisease)


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YRange, YMax):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot proportions of disease genes only
    # Create a horizontal bar plot for proportions of disease genes
    ax.bar([0, 0.3], Data, width = 0.2, label = 'disease', color= ['black'] + ['lightgrey'] * 6, linewidth = 0.7)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write axis labels    
    YLabel = 'Proportion of disease genes'
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # edit y axis ticks
    plt.yticks(YRange)    
    # add a range for the Y axis
    plt.ylim([0, YMax])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right = 'off', left = 'on', labelbottom='off', colors = 'black',
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
ax1 = CreateAx(5, 1, 1, fig, DisProp[0], 'complex diseases', np.arange(0, 1, 0.1), 0.8)
ax2 = CreateAx(5, 1, 2, fig, DisProp[1], 'GWAS', np.arange(0, 0.30, 0.05), 0.25)
ax3 = CreateAx(5, 1, 3, fig, DisProp[2], 'tumor drivers', np.arange(0, 0.1, 0.02), 0.08)
ax4 = CreateAx(5, 1, 4, fig, DisProp[3], 'medelian diseases', np.arange(0, 0.40, 0.05), 0.351)
ax5 = CreateAx(5, 1, 5, fig, DisProp[4], 'all diseases', np.arange(0, 1, 0.1), 0.8)


# annotate figure to add significance
ax1 = AddSignificance(ax1, PVals[0], 0.1, 0.4, 0.73, 0.25, 0.78)
ax2 = AddSignificance(ax2, PVals[1], 0.1, 0.4, 0.22, 0.25, 0.24)
ax3 = AddSignificance(ax3, PVals[2], 0.1, 0.4, 0.07, 0.25, 0.075)
ax4 = AddSignificance(ax4, PVals[3], 0.1, 0.4, 0.32, 0.25, 0.34)
ax5 = AddSignificance(ax5, PVals[4], 0.1, 0.4, 0.75, 0.25, 0.79)


# add subplot labels
ax1.text(-0.45, 0.85, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(0.8, 0.85, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(2.2, 0.85, 'C', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(3.6, 0.85, 'D', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
ax1.text(5, 0.85, 'E', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)


# add legend
N = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'intronless internal genes')
D = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'internal genes with introns')
ax1.legend(handles = [D, N], loc = (1, 1.15), fontsize = 6, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()


# save figure to file
fig.savefig('truc.pdf', bbox_inches = 'tight')


## compare proportion of disease genes among external genes with different orientation with their internal genes
#
#AllSame, AllOpposite = set(), set()
#IntronlessSame, IntronlessOpposite = set(), set()
#IntronSame, IntronOpposite =set(), set()
#
## loop over host-nested transcript pairs
#for i in range(len(HostNestedPairs)):
#    # check that both transcripts have coordinates and have corresponding gene names    
#    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
#    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
#    # get the orientation of the external and internal transcripts
#    OrientationPair = GenePairOrientation(HostNestedPairs[i], TranscriptCoordinates)
#    # populate sets with external genes
#    if len(set(OrientationPair)) == 1:
#        assert set(OrientationPair) == {'-'} or set(OrientationPair) == {'+'}
#        AllSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        # check if the internal gene has introns
#        if HostNestedPairs[i][1] in IntronCoord:
#            IntronSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        else:
#            IntronlessSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#    elif len(set(OrientationPair)) == 2:
#        assert set(OrientationPair) == {'-', '+'}
#        AllOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        # check if the internal gene has introns
#        if HostNestedPairs[i][1] in IntronCoord:
#            IntronOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        else:
#            IntronlessOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#       
## make a list of external genes
#ExternalGenes = [AllSame, AllOpposite, IntronSame, IntronOpposite, IntronlessSame, IntronlessOpposite]  
#
## count disease and non-disease genes    
#GADCounts = CountDiseaseGenes(ExternalGenes, GAD)    
#GWASCounts = CountDiseaseGenes(ExternalGenes, GWAS)
#DriversCounts = CountDiseaseGenes(ExternalGenes, Drivers)
#OMIMCounts = CountDiseaseGenes(ExternalGenes, OMIM)
#AllCounts = CountDiseaseGenes(ExternalGenes, DiseaseGenes)
#
## compare the proportion of disease genes for each set of disease gene
#PVals = []
#counts = [GADCounts, GWASCounts, DriversCounts, OMIMCounts, AllCounts]
#for i in range(len(counts)):
#    pvalues = []
#    for j in range(0, len(counts[i]), 2):
#        p = stats.fisher_exact([counts[i][j], counts[i][j+1]])[1]
#        pvalues.append(p)
#    PVals.append(pvalues)    
#
## create lists of proportions for disease and non-disease genes
#DisProp, NonDisProp = [], []
#for i in range(len(counts)):
#    disease, nondisease = GetProportions(counts[i])
#    DisProp.append(disease)
#    NonDisProp.append(nondisease)
#
#Origins = ['GAD', 'GWAS', 'Drivers', 'OMIM', 'All']
#GeneCats = ['External_Same', 'External_opposite', 'WithIntronSame', 'WithIntronOpposite', 'IntronlessSame', 'IntronlessOpposite']
#
#
#newfile = open('truc.txt', 'w')
#newfile.write('Proportion of disease genes among external genes depending on their orientation\n')
#newfile.write('with their intronless or intron-containing internal genes\n\n')
#header = ['Disease genes', 'Gene category', 'N disease genes', 'N non-disease genes', 'Proportion disease genes', 'P']
#newfile.write('\t'.join(header) + '\n')
#
#
## loop over the gene counts for each disease origin
#for i in range(len(counts)):
#    # set up variable to get the index of the pvalue list (the list doesn't have same length)
#    m = 0
#    for j in range(len(counts[i])):
#        line = [Origins[i], GeneCats[j], str(counts[i][j][0]), str(counts[i][j][1]), str(round(DisProp[i][j] * 100, 2))]
#        # add p value on the line of the opposite orientation
#        if j % 2 != 0:
#            # update variable m
#            m += 1
#            # get index of the p value in list
#            k = j -m
#            line.append(str(PVals[i][k]))
#        newfile.write('\t'.join(line) + '\n')
#newfile.close()        






















##########################################
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#fig.savefig(outputfile + '.pdf', bbox_inches = 'tight')
#fig.savefig(outputfile + '.eps', bbox_inches = 'tight')


## compare proportion of disease genes among external genes with different orientation with their internal genes
#
#AllSame, AllOpposite = set(), set()
#IntronlessSame, IntronlessOpposite = set(), set()
#IntronSame, IntronOpposite =set(), set()
#
## loop over host-nested transcript pairs
#for i in range(len(HostNestedPairs)):
#    # check that both transcripts have coordinates and have corresponding gene names    
#    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
#    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
#    # get the orientation of the external and internal transcripts
#    OrientationPair = GenePairOrientation(HostNestedPairs[i], TranscriptCoordinates)
#    # populate sets with external genes
#    if len(set(OrientationPair)) == 1:
#        assert set(OrientationPair) == {'-'} or set(OrientationPair) == {'+'}
#        AllSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        # check if the internal gene has introns
#        if HostNestedPairs[i][1] in IntronCoord:
#            IntronSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        else:
#            IntronlessSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#    elif len(set(OrientationPair)) == 2:
#        assert set(OrientationPair) == {'-', '+'}
#        AllOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        # check if the internal gene has introns
#        if HostNestedPairs[i][1] in IntronCoord:
#            IntronOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        else:
#            IntronlessOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#       
## make a list of external genes
#ExternalGenes = [AllSame, AllOpposite, IntronSame, IntronOpposite, IntronlessSame, IntronlessOpposite]  
#
## count disease and non-disease genes    
#GADCounts = CountDiseaseGenes(ExternalGenes, GAD)    
#GWASCounts = CountDiseaseGenes(ExternalGenes, GWAS)
#DriversCounts = CountDiseaseGenes(ExternalGenes, Drivers)
#OMIMCounts = CountDiseaseGenes(ExternalGenes, OMIM)
#AllCounts = CountDiseaseGenes(ExternalGenes, DiseaseGenes)
#
## compare the proportion of disease genes for each set of disease gene
#PVals = []
#counts = [GADCounts, GWASCounts, DriversCounts, OMIMCounts, AllCounts]
#for i in range(len(counts)):
#    pvalues = []
#    for j in range(0, len(counts[i]), 2):
#        p = stats.fisher_exact([counts[i][j], counts[i][j+1]])[1]
#        pvalues.append(p)
#    PVals.append(pvalues)    
#
## create lists of proportions for disease and non-disease genes
#DisProp, NonDisProp = [], []
#for i in range(len(counts)):
#    disease, nondisease = GetProportions(counts[i])
#    DisProp.append(disease)
#    NonDisProp.append(nondisease)
#
#Origins = ['GAD', 'GWAS', 'Drivers', 'OMIM', 'All']
#GeneCats = ['External_Same', 'External_opposite', 'WithIntronSame', 'WithIntronOpposite', 'IntronlessSame', 'IntronlessOpposite']
#
#
#newfile = open('DiseaseEnrichExternal.txt', 'w')
#newfile.write('Proportion of disease genes among external genes depending on their orientation\n')
#newfile.write('with their intronless or intron-containing internal genes\n\n')
#header = ['Disease genes', 'Gene category', 'N disease genes', 'N non-disease genes', 'Proportion disease genes', 'P']
#newfile.write('\t'.join(header) + '\n')
#
#
## loop over the gene counts for each disease origin
#for i in range(len(counts)):
#    # set up variable to get the index of the pvalue list (the list doesn't have same length)
#    m = 0
#    for j in range(len(counts[i])):
#        line = [Origins[i], GeneCats[j], str(counts[i][j][0]), str(counts[i][j][1]), str(round(DisProp[i][j] * 100, 2))]
#        # add p value on the line of the opposite orientation
#        if j % 2 != 0:
#            # update variable m
#            m += 1
#            # get index of the p value in list
#            k = j -m
#            line.append(str(PVals[i][k]))
#        newfile.write('\t'.join(line) + '\n')
#newfile.close()        
