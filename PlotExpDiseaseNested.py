# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:11:00 2017

@author: RJovelin
"""

# use this script to plot expression divergence between host and nested genes 
# for disease and non-disease genes 

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
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
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
HostNestedTSPairs = GetHostNestedPairs(Matches)
# make a a list of host, nested gene pairs
NestedPairs = GetHostNestedPairs(Nested)

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

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

AllPairs = []
for disease in [GAD, GWAS, Drivers, OMIM, DiseaseGenes]:
    # create list with both, external, internal and none being disease genes
    atleast, none = [], []
    # loop over gene pairs
    for pair in NestedPairs:
        # check that both genes are expressed
        if pair[0] in ExpressionProfile and pair[1] in ExpressionProfile:
            # check which genes are disease genes
            if pair[0] in disease or pair[1] in disease:
                atleast.append(pair)
            elif pair[0] not in disease and pair[1] not in disease:
                none.append(pair)
    AllPairs.append([atleast, none])

# make a list with disease types
DiseaseClasses = ['complex', 'GWAS', 'Drivers', 'mendelian', 'all']

# filter diseases if any pair count < 20
DiseaseToRemove, PairsToRemove = [], []
for i in range(len(AllPairs)):
    # check each pair count for that disease
    for j in range(len(AllPairs[i])):
        if len(AllPairs[i][j]) < 10:
            if AllPairs[i] not in PairsToRemove:
                PairsToRemove.append(AllPairs[i])
            if DiseaseClasses[i] not in DiseaseToRemove:
                DiseaseToRemove.append(DiseaseClasses[i])
#remove lists and diseases      
for i in DiseaseToRemove:
    DiseaseClasses.remove(i)
for pairs in PairsToRemove:
    AllPairs.remove(pairs)
# do QC
for i in range(len(AllPairs)):
    for j in range(len(AllPairs[i])):
        assert len(AllPairs[i][j]) >= 20
   

for i in range(len(AllPairs)):
    print(i, list(map(lambda x: len(x), AllPairs[i])))

    
## generate gene pairs of external and internal genes for intronless and intron-containing internal genes
#PairsWithIntrons, PairsNoIntrons = [], []
## loop over host-nested transcript pairs
#for i in range(len(HostNestedTSPairs)):
#    # check that both transcripts have coordinates and have corresponding gene names    
#    assert HostNestedTSPairs[i][0] in MapTranscriptGene and HostNestedTSPairs[i][0] in TranscriptCoordinates    
#    assert HostNestedTSPairs[i][1] in MapTranscriptGene and HostNestedTSPairs[i][1] in TranscriptCoordinates    
#    if HostNestedTSPairs[i][1] in IntronCoord:
#        # internal gene has introns, add gene pairs to list 
#            PairsWithIntrons.append([MapTranscriptGene[HostNestedTSPairs[i][0]], MapTranscriptGene[HostNestedTSPairs[i][1]]])
#    else:
#        # internal gene is intronless, add gene pair to list
#        PairsNoIntrons.append([MapTranscriptGene[HostNestedTSPairs[i][0]], MapTranscriptGene[HostNestedTSPairs[i][1]]])
#
## make a list of lists of gene pairs
#IntronPairs = [PairsWithIntrons, PairsNoIntrons]
## remove gene pairs if any gene in the pair lacks expression
#for i in range(len(IntronPairs)):
#    IntronPairs[i] = FilterGenePairsWithoutExpression(IntronPairs[i], ExpressionProfile, 'strict')
## add gene pairs defined by distance to list
#for L in [Proximal, Moderate, Intermediate, Distant]:
#    IntronPairs.append(L)

## generate gene pairs of external and internal genes with same and opposite orientation
#Same, Opposite = [], []
#for pair in NestedPairs:
#    orientation = GenePairOrientation(pair, GeneCoord)
#    if len(set(orientation)) == 2:
#        Opposite.append(pair)
#    elif len(set(orientation)) == 1:
#        Same.append(pair)
#
## make a list of lists of gene pairs
#OrientationPairs = [Same, Opposite]
## remove gene pairs if any gene in the pair lacls expression
#for i in range(len(OrientationPairs)):
#    OrientationPairs[i] = FilterGenePairsWithoutExpression(OrientationPairs[i], ExpressionProfile, 'strict')
## add gene pairs defined by distance to list
#for L in [Proximal, Moderate, Intermediate, Distant]:
#    OrientationPairs.append(L)
#



# compute expression divergence between pairs of genes
Divergence = []
for i in range(len(AllPairs)):
    # access lists of genes for each disease
    Div = []
    for j in range(len(AllPairs[i])):
        # access list of gene pairs
        D = ComputeExpressionDivergenceGenePairs(AllPairs[i][j], ExpressionProfile)
        Div.append(D)
    Divergence.append(Div)
    
# make a list of gene category names parallel to the list of gene pairs
GeneCat = ['both', 'ext', 'int', 'none']

# create lists with means and SEM for each gene category
Means, SEMs = [], []
for i in range(len(Divergence)):
    # create lists with means and sem
    average, error = [], []
    for j in range(len(Divergence[i])):
        average.append(np.mean(Divergence[i][j]))
        error.append(np.std(Divergence[i][j]) / math.sqrt(len(Divergence[i][j])))
    Means.append(average)
    SEMs.append(error)        


# perform statistical tests between gene categories
PValues = []
for i in range(len(Divergence)):
    Pvals = []
    # loop over inner lists    
    for j in range(0, len(Divergence[i]) -1):
        for k in range(j+1, len(Divergence[i])):
            P = PermutationResampling(Divergence[i][j], Divergence[i][k], 1000, statistic = np.mean)
            Pvals.append(P)
    PValues.append(Pvals)

for i in range(len(PValues)):
    print(PValues[i])

print(Means)
print(DiseaseClasses)

#
#
#
#
#
#        
## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, figure, Data, GeneCats, Title, YRange, YMax, isYLabel):
#    '''
#    return an ax object part of figure
#    '''
#
#    # add a plot to figure (N row, N column, plot N)
#    ax = fig.add_subplot(Rows, Columns, Position)
#    # set colors
#    colorscheme = ['#225ea8', '#e31a1c', 'lightgrey', 'lightgrey']
#    # plot nucleotide divergence
#    ax.bar([0.05, 0.35, 0.65, 0.95], Data[0], 0.2, yerr = Data[1], color = colorscheme,
#           edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write y axis label
#    if isYLabel == True:
#        ax.set_ylabel('Expression\ndivergence', color = 'black',  size = 7, ha = 'center', **FigFont)
#    # add ticks and lebels
#    plt.xticks([0.15, 0.45, 0.75, 1.05], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
#    # add title
#    ax.set_title(Title, color = 'black', size = 7, ha = 'center', **FigFont)    
#    # add a range for the Y and X axes
#    plt.ylim([0, YMax])
#    plt.xlim([0, 1.2])
#    # edit y axis ticks
#    plt.yticks(YRange) 
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)
#    if isYLabel == True:
#        ax.spines["left"].set_visible(True)
#    else:
#        ax.spines["left"].set_visible(False)
#    # edit tick parameters
#    if isYLabel == True:
#        plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                        right = 'off', left = 'on', labelbottom='on',
#                        colors = 'black', labelsize = 7, direction = 'out')  
#    else:
#        plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                        right = 'off', left = 'off', labelbottom='on', labelleft='off',
#                        colors = 'black', labelsize = 7, direction = 'out')  
#    # Set the tick labels font name
#    for label in ax.get_yticklabels():
#        label.set_fontname('Arial')   
#    return ax  
#
#
## create figure
#fig = plt.figure(1, figsize = (2.5, 2.2))
#
#ax1 = CreateAx(1, 2, 1, fig, [MeanIntron, SEMIntron], GeneCatIntrons, 'Introns in internal genes', np.arange(0, 1.2, 0.2), 1)
#ax2 = CreateAx(1, 2, 2, fig, [MeanOrientation, SEMOrientation], GeneCatOrientation, 'Strand orientation', np.arange(0, 1.2, 0.2), 1)
#
#
#
## add subplot labels
#ax1.text(-0.55, 1.2, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#ax2.text(-0.55, 1.2, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#
## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#DiffIntrons = ['A', 'B', 'C', 'D', 'E', 'A']
#DiffOrientation = ['A,E', 'A', 'B', 'C', 'D', 'E']
#ypos = [0.7] * 6
#xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65]
#for i in range(len(DiffIntrons)):
#     ax1.text(xpos[i], ypos[i], DiffIntrons[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#for i in range(len(DiffOrientation)):
#     ax2.text(xpos[i], ypos[i], DiffOrientation[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#
## make sure subplots do not overlap
#plt.tight_layout()    
#    
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')
#
#
#
#
