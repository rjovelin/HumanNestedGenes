# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:39:25 2017

@author: RJovelin
"""

# use this script to plot conservation of overlapping genes


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

# make a list of json files
files = ['OverlappingGenes.json', 'NestedGenes.json', 'PiggyBackGenes.json', 'ConvergentGenes.json', 'DivergentGenes.json']

# make a list of GFF files
GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
             'Gorilla_gorilla.gorGor3.1.88.gff3', 'Macaca_mulatta.Mmul_8.0.1.88.gff3',
             'Pongo_abelii.PPYG2.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
             'Mus_musculus.GRCm38.88.gff3', 'Bos_taurus.UMD3.1.88.gff3',
             'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
             'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
             'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
             'Monodelphis_domestica.BROADO5.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3',
             'Sorex_araneus.COMMON_SHREW1.88.gff3']
# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
           'Mouse', 'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
           'Cat', 'Opossum', 'Platypus', 'Shrew']

# make parallel lists of json files with each overlapping classes [human, mouse]
NestedFiles = [i + 'NestedGenes.json' for i in ['Human', 'Mouse']]
PiggybackFiles = [i + 'PiggyBackGenes.json' for i in ['Human', 'Mouse']]
ConvergentFiles = [i + 'ConvergentGenes.json' for i in ['Human', 'Mouse']]
DivergentFiles = [i + 'DivergentGenes.json' for i in ['Human', 'Mouse']]
# make a parallel list of ortholog files
OrthoFiles = ['Human' + i + 'Orthologs.txt' for i in Species[1:]]

# make lists of dictionaries for each type of overlapping gene
AllOverlapGenes  = []
# loop over files
AllFiles = [NestedFiles, PiggybackFiles, ConvergentFiles, DivergentFiles]
for i in range(len(AllFiles)):
    OvDicts = []
    for j in range(len(AllFiles[i])):
        # load dictionary from json file
        json_data = open(AllFiles[i][j])
        OvDicts.append(json.load(json_data))
        json_data.close()
    AllOverlapGenes.append(OvDicts)
for i in AllOverlapGenes:
    assert len(AllOverlapGenes[i]) == 2
    
# make a list of dicts with orthologs
AllOrthologs = []
for i in range(len(OrthoFiles)):
    orthologs = MatchOrthologs(OrthoFiles[i])
    AllOrthologs.append(orthologs)
    






# 1) determine overlapping conservation with mouse

# for each overlapping gene class, count the proportion of gene pairs with both genes
# in the same configuration, 1 gene only or none

HumanMouseConservation = []


# get nested pairs 
NestedPairs = [GetHostNestedPairs(AllNestedGenes[i]) for i in range(len(AllNestedGenes))]
PgkPairs = [GetHostNestedPairs(AllPbkGenes[i]) for i in range(len(AllPbkGenes))]
ConPairs = [GetHostNestedPairs(AllConGenes[i]) for i in range(len(AllConGenes))]
DivPairs = [GetHostNestedPairs(AllDivGenes[i]) for i in range(len(AllDivGenes))]

for i in [NestedPairs, PgkPairs, ConPairs, DivPairs]:
    print(list(map(lambda x: len(x), i)))

















# make a list of gene pairs
AllPairs = []
for i in range(len(Overlap)):
    AllPairs.append(GetHostNestedPairs(Overlap[i]))

# make a list of gene sets
AllGenes = []
for i in range(len(Overlap)):
    AllGenes.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)


# make a list of host and testd transcript pairs
HostNestedTSPairs = GetHostNestedPairs(Matches)

# make sets of external genes with intronless internal genes or intron-containing genes
# make sets of internal intronless genes and internal genes with introns
ExternalIntronless, ExternalWithIntron = set(), set()
InternalIntronless, InternalWithIntron = set(), set()
for i in range(len(HostNestedTSPairs)):
    # get the external and internal transcripts of the pair
    external, internal = HostNestedTSPairs[i][0], HostNestedTSPairs[i][1]
    # check if internal transcript has intron
    if internal in IntronCoord:
        ExternalWithIntron.add(MapTranscriptGene[external])
        InternalWithIntron.add(MapTranscriptGene[internal])
    else:
        ExternalIntronless.add(MapTranscriptGene[external])
        InternalIntronless.add(MapTranscriptGene[internal])

# create lists of nested gene pairs with same and opposite directions
same, opposite = [], []
for pair in AllPairs[1]:
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


# check options to select the expression file
if Species == 'human':
    if Breadth == 'full':
        # parse the GTEX expression summary file to obtain the expression profile of each gene
        ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
    elif Breadth == 'restricted':
        # parse the GTEX expression summary file to obtain the expression profile of each gene
        ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
        # match expression profiles between mouse and human 
        ExpressionProfile = MatchHumanToMouseExpressionProfiles(ExpressionProfile)
    elif Breadth == 'narrow':
        # parse expression data in common between chimp and human
        ExpressionProfile = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
elif Species == 'chimp':
    ExpressionProfile = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
elif Species == 'mouse':
    ExpressionProfile = ParseExpressionFile('Mouse_Median_Normalized_FPKM.txt')

# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# make a list of all gene sets
if GenesOfInterest == 'overlapping':
    AllGeneSets = [NonOverlappingGenes, AllGenes[1], AllGenes[2], AllGenes[3], AllGenes[4]]
elif GenesOfInterest == 'external':
    AllGeneSets = [NonOverlappingGenes, InternalIntronless, InternalWithIntron, ExternalIntronless, ExternalWithIntron]
elif GenesOfInterest == 'nested':
    AllGeneSets = [NonOverlappingGenes, InternalSameGenes, InternalOppositeGenes, ExternalSameGenes, ExternalOppositeGenes]

# remove genes that do not have expression profile
for i in range(len(AllGeneSets)):
    # make a list of genes to remove
    to_remove = [gene for gene in AllGeneSets[i] if gene not in ExpressionProfile]
    for gene in to_remove:
        AllGeneSets[i].remove(gene)

# make a list of tissues

if Species == 'human':
    if Breadth == 'full':
        infile = open('GTEX_Median_Normalized_FPKM.txt')
        header = infile.readline().rstrip().split('\t')
        infile.close()
        Tissues = header[1:]
    elif Breadth == 'restricted':
        # get the list of tissues in common between human and mouse
        infile = open('Mouse_Median_Normalized_FPKM.txt')        
        Tissues = infile.readline().rstrip().split('\t')
        infile.close()        
        Tissues = Tissues[1:]
    elif Breadth == 'narrow':
        Tissues = ['Brain', 'Cerebellum', 'Heart', 'Kidney', 'Liver', 'Testis']
elif Species == 'chimp':
    Tissues = ['Brain', 'Cerebellum', 'Heart', 'Kidney', 'Liver', 'Testis']
elif Species == 'mouse':
    infile = open('Mouse_Median_Normalized_FPKM.txt')        
    Tissues = infile.readline().rstrip().split('\t')
    infile.close()    
    Tissues = Tissues[1:]
 
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
if GenesOfInterest == 'overlapping':
    GeneCats = ['NoOv', 'Nst', 'Pbk', 'Conv', 'Div']
elif GenesOfInterest == 'external':
    GeneCats = ['NoOv', 'IntN', 'IntW', 'ExtN', 'ExtW']
elif GenesOfInterest == 'nested':
    # make a list of gene category names parallel to the list of gene sets
    GeneCats = ['NoOv', 'CisInt', 'TransInt', 'CisExt', 'TransExt']

# create a dict to count the number of genes in each category with highest expression in each tissue
HighestExpression = {}
# inititialize dict with list of 0s
for i in GeneCats:
    HighestExpression[i] = [0] * len(Tissues)

# loop over the gene sets
for i in range(len(AllGeneSets)):
    # loop over each gene in given set
    for gene in AllGeneSets[i]:
        # check if gene has expression profile
        assert gene in ExpressionProfile
        # get the index of the maximum expression value
        pos = ExpressionProfile[gene].index(max(ExpressionProfile[gene]))
        # check that there is a single maximum expression value
        assert ExpressionProfile[gene].count(max(ExpressionProfile[gene])) == 1, '> 1 max value'
        # update counter at pos index
        HighestExpression[GeneCats[i]][pos] += 1


# divide by the total number of genes in each category to get proportions
for i in range(len(GeneCats)):
    for j in range(len(HighestExpression[GeneCats[i]])):
        HighestExpression[GeneCats[i]][j] = round(HighestExpression[GeneCats[i]][j] / len(AllGeneSets[i]), 4)

# create a dictionary with tissue as key and a list of gene proportions for each gene category as value
Proportions = {}
for i in range(len(Tissues)):
    Proportions[Tissues[i]] = []
    for j in range(len(GeneCats)):
        Proportions[Tissues[i]].append(HighestExpression[GeneCats[j]][i])





# 2) determine the proportion of gene in each gene category that have orthologs in each other species

# load the dictionary of human overlapping genes
# load dictionary of overlapping gene pairs
json_data = open('HumanOverlappingGenes.json')
HumanOverlapping = json.load(json_data)
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
# Map Transcript names to gene names {transcript: gene}
MapTranscriptGene = TranscriptToGene(GFF)
# make a set of un-nested genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(HumanOverlapping, GeneCoord)

# make a set of external genes




# make a set of internal genes



# make a set of piggyback genes



# make a set of convergent genes



# make a set of divergent genes

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, colorscheme, XLabel, YMax, YLabel):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot nucleotide divergence
    ax.bar([0, 0.2, 0.4, 0.6, 0.8], Data, 0.2, color = colorscheme,
           edgecolor = 'black', linewidth = 0.5)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    if YLabel == True:
        ax.set_ylabel('Proportion', color = 'black',  size = 7, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    
    # add a range for the Y axis
    plt.ylim([0, YMax])
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)
    if YLabel == True:
        ax.spines["left"].set_visible(True)  
    elif YLabel == False:
        ax.spines['left'].set_visible(False)
    
    # edit tick parameters    
    if YLabel == True:
        plt.tick_params(axis='both', which='both', bottom='off', top='off',
                        right = 'off', left = 'on', labelbottom='off', 
                        colors = 'black', labelsize = 7, direction = 'out')  
    elif YLabel == False:
        plt.tick_params(axis='both', which='both', bottom='off', top='off',
                        right = 'off', left = 'off', labelleft='off', labelbottom='off', 
                        colors = 'black', labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    return ax      


# create figure
fig = plt.figure(1, figsize = (5, 3))

# plot data
j = 1
# get y axis range and colors
if GenesOfInterest == 'overlapping':
    YMax = 0.31
elif GenesOfInterest == 'external':
    YMax = 0.61
elif GenesOfInterest == 'nested':
    YMax = 0.41

# create legend
if Species == 'human':
    if Breadth == 'narrow':
        Position = 6
        Xpos, Ypos = -5, 0.7
    else:
        Position = 10
        Xpos, Ypos = -10, 0.6
elif Species == 'chimp':
    Position = 6
    Xpos, Ypos = -5, 0.7
elif Species == 'mouse':
    Position = 10
    Xpos, Ypos = -10, 0.6

# set up colors
colorscheme = ['#fb9a99', '#9ecae1','#3182bd', '#a1d99b','#31a354']

for i in range(len(Tissues)):
    tissue = Tissues[i]
    if j == 1 or j == 11 or j == 21:
        YLabel = True
    else:
        YLabel = False
    ax = CreateAx(10, 3, j, fig, Proportions[tissue], colorscheme, tissue.lower().replace('_', '\n'), 0.61, YLabel)
    j += 1

    if j == Position:
        if GenesOfInterest == 'external':
            Labels = [['#fb9a99', 'NoOv'], ['#9ecae1', 'IntN'], ['#3182bd', 'IntW'], ['#a1d99b', 'ExtN'], ['#31a354', 'ExtW']]
        elif GenesOfInterest == 'overlapping':
            Labels = [['#fb9a99', 'NoOv'], ['#a6cee3', 'Nst'], ['#1f78b4', 'Pbk'], ['#b2df8a', 'Con'], ['#33a02c', 'Div']] 
        elif GenesOfInterest == 'nested':
            Labels = [['#fb9a99', 'NoOv'], ['#9ecae1', 'CisInt'], ['#3182bd', 'TransInt'], ['#a1d99b', 'CisExt'], ['#31a354', 'TransExt']]
        a = mpatches.Patch(facecolor = Labels[0][0], edgecolor = 'black', linewidth = 0.5, label= Labels[0][1])
        b = mpatches.Patch(facecolor = Labels[1][0], edgecolor = 'black', linewidth = 0.5, label= Labels[1][1])
        c = mpatches.Patch(facecolor = Labels[2][0], edgecolor = 'black', linewidth = 0.5, label= Labels[2][1])
        d = mpatches.Patch(facecolor = Labels[3][0], edgecolor = 'black', linewidth = 0.5, label= Labels[3][1])
        e = mpatches.Patch(facecolor = Labels[4][0], edgecolor = 'black', linewidth = 0.5, label= Labels[4][1])
        ax.legend(handles = [a, b, c, d, e], bbox_to_anchor=(Xpos, Ypos), loc = 3, fontsize = 6, frameon = False, ncol = 5)
    
# adjust padding between subplots
# pad controls the padding around the figure border
# hpad and wpad control the padding between subplots
plt.tight_layout(pad=0.2, w_pad=0, h_pad=0.5)



fig.savefig('truc.pdf', bbox_inches = 'tight')




