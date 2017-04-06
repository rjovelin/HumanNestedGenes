# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 21:24:45 2017

@author: Richard
"""

# use this script to plot differences in functional similarities between overlapping genes

# usage PlotFunctionalSimilarity
#- [all/molfunc/celcomp/biolproc]: use all GO terns, molecular function, cellular compartments or biological processes

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


#get option from command
GOClass = sys.argv[1]
assert GOClass in ['all', 'molfunc', 'celcomp', 'biolproc']

# load dictionaries of overlapping genes
JsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',
             'HumanPiggyBackGenes.json', 'HumanConvergentGenes.json',
             'HumanDivergentGenes.json']
# get GFF file
GFF = 'Homo_sapiens.GRCh38.88.gff3'

# make a list of dictionaries
Overlap = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    Overlap.append(overlapping)

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

# make pairs of overlapping genes
OverlappingPairs = []
for i in range(len(Overlap)):
    pairs = GetHostNestedPairs(Overlap[i])
    OverlappingPairs.append(pairs)

# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)

# map gene names with gene ID {gene ID: name}
Names = MapNametoID(GFF)
# parse file with GO annotations 
GOAnnotations = ParseGOFile('goa_human.gaf')
# map gene IDs to GO annotations            
GeneOntology = MapEnsemblGenesToGOTerms(GOAnnotations, Names)
# check option to see GO ids are filtered or not
if GOClass == 'molfunc':
    GOSubSet = FilterGOTerms('GOMolecularFunction.txt')
elif GOClass == 'celcomp':
    GOSubSet = FilterGOTerms('GOCellularCompartments.txt')
elif GOClass == 'biolproc':
    GOSubSet = FilterGOTerms('GOBiologicalProcesses.txt')
# keep only terms for specific GO class
if GOClass != 'all':
    for gene in GeneOntology:
         GeneOntology[gene] = set(filter(lambda x: x in GOSubSet, GeneOntology[gene]))
# remove genes if gene has no GO term
to_remove = [gene for gene in GeneOntology if GeneOntology[gene] == 0]
if len(to_remove) != 0:
    for gene in to_remove:
        del GeneOntology[gene]

# remove pairs if genes do not have GO terms
for i in range(len(OverlappingPairs)):
    to_remove = [pair for pair in OverlappingPairs[i] if pair[0] not in GeneOntology or pair[1] not in GeneOntology]
    for pair in to_remove:
        OverlappingPairs[i].remove(pair)
    
# create a list of lists with functional similarities between genes of the same pair
FunctionalSimilarity = []
# include only the 4 overlapping types
for i in range(1, len(OverlappingPairs)):
    # compute jaccard similarity index between each pair
    GOOverlap = []
    for pair in OverlappingPairs[i]:
        if pair[0] in GeneOntology and pair[1] in GeneOntology:
            JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
            GOOverlap.append(JI)
    FunctionalSimilarity.append(GOOverlap)
 
# generate random pairs of non-overlapping genes as a baseline random JI
# make a list of non-overlapping genes convert set of non-overlapping genes to list
NonOverlappingGenes = list(NonOverlappingGenes)
BaseLine = []
# draw 10000 random pairs
replicates = 10000
while replicates != 0:
    # pick 2 random genes
    j = random.randint(0, len(NonOverlappingGenes) -1)
    k = random.randint(0, len(NonOverlappingGenes) -1)
    gene1, gene2 = NonOverlappingGenes[j], NonOverlappingGenes[k]
    if gene1 in GeneOntology and gene2 in GeneOntology:
        JI = JaccardIndex(GeneOntology[gene1], GeneOntology[gene2])
        BaseLine.append(JI)
        replicates -= 1

OverlapTypes = ['Nested', 'PiggyBack', 'Convergent', 'Divergent']
for i in range(len(FunctionalSimilarity)):
    P = PermutationResampling(BaseLine, FunctionalSimilarity[i], 1000, statistic = np.mean)
    print('Baseline', OverlapTypes[i], len(BaseLine), len(FunctionalSimilarity[i]), np.mean(BaseLine), np.mean(FunctionalSimilarity[i]), P)

# make a set of overlapping genes
Overlapping = MakeFullPartialOverlapGeneSet(Overlap[0])

# generate a dict to draw genes {chromo: {num: gene}}    
ToDrawGenesFrom = {}
# loop over chromosomes
for chromo in OrderedGenes:
    # set up counter
    k = 0
    # add chromo as key and intialize inner dict
    ToDrawGenesFrom[chromo] = {}
    # loop over the list of ordered genes
    for i in range(len(OrderedGenes[chromo])):
        # check that gene does not overlap with any other gene
        if OrderedGenes[chromo][i] not in Overlapping and OrderedGenes[chromo][i] in GeneOntology:
            # add gene pair and update counter
            ToDrawGenesFrom[chromo][k] = OrderedGenes[chromo][i]
            k += 1

# generate matching control pairs for each type of overllaping class
MatchingControl = []
# loop over list of overlapping pairs
for i in range(1, len(OverlappingPairs)):
    control  = []
    for pair in OverlappingPairs[i]:
        # make a list of matching gene pairs (orientation, chromosome, distance)
        PairPool = GenerateMatchingPoolPairs(pair, ToDrawGenesFrom, GeneCoord, 2000)
        if len(PairPool) != 0:
            # draw a matching gene pair at random
            j = random.randint(0, len(PairPool) -1)
            assert PairPool[j][0] not in Overlapping and PairPool[j][1] not in Overlapping
            assert PairPool[j][0] in NonOverlappingGenes and PairPool[j][1] in NonOverlappingGenes
            control.append(PairPool[j])
    MatchingControl.append(control)
# do some QC
for item in MatchingControl:
    assert len(item) != 0

# compute the jaccard similarity index for each pair of each controls
SimilarityControls = []
for i in range(len(MatchingControl)):
    controlJI = []
    for pair in MatchingControl[i]:
        JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
        controlJI.append(JI)
    SimilarityControls.append(controlJI)

# make a list of P values of permutation test between each control and each ovverlapping type
PVals = []
for i in range(len(FunctionalSimilarity)):
    P = PermutationResampling(SimilarityControls[i], FunctionalSimilarity[i], 1000, statistic = np.mean)
    print(OverlapTypes[i], len(SimilarityControls[i]), len(FunctionalSimilarity[i]), np.mean(SimilarityControls[i]), np.mean(FunctionalSimilarity[i]), P)
    PVals.append(P)
    
# create lists with means JI and SEM for each gene category and its control [[control_nested], [nested]...]
MeanExpDiv, SEMExpDiv = [], []
for i in range(len(FunctionalSimilarity)):
    MeanExpDiv.append(np.mean(SimilarityControls[i]))
    MeanExpDiv.append(np.mean(FunctionalSimilarity[i]))
    SEMExpDiv.append(np.std(SimilarityControls[i]) / math.sqrt(len(SimilarityControls[i])))
    SEMExpDiv.append(np.std(FunctionalSimilarity[i]) / math.sqrt(len(FunctionalSimilarity[i])))


# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# set colors
colorscheme = ['black','lightgrey'] * 4
# plot functional similarity
xpos = [0.05, 0.25, 0.65, 0.85, 1.25, 1.45, 1.85, 2.05]
# plot nucleotide divergence
ax.bar(xpos, MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
       edgecolor = 'black', linewidth = 0.5, error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Functional similarity', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.25, 0.85, 1.45, 2.05], OverlapTypes, size = 7, color = 'black', ha = 'center', **FigFont)

# add a range for the Y and X axes
plt.ylim([0, 1])
#plt.xlim([0, 2.45])

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

# plot the baseline
plt.plot((0, 2.05), (np.mean(BaseLine), np.mean(BaseLine)), color = 'grey', linestyle = '--', linewidth = 0.5)
      
# convert p-values to star significance level
Significance = []
for pvalue in PVals:
    if pvalue >= 0.05:
        Significance.append('')
    elif pvalue < 0.05 and pvalue >= 0.01:
        Significance.append('*')
    elif pvalue < 0.01 and pvalue >= 0.001:
        Significance.append('**')
    elif pvalue < 0.001:
        Significance.append('***')

# make a list of x, y axis positions for brackets
XLinePos = list(map(lambda x: x + 0.1, xpos))
YLinePos = [0.12, 0.50, 0.12, 0.12]
# make a list of x, y axis positions for stars
XStarPos = [0.25, 0.85, 1.45, 1.85]
YStarPos = [0.14, 0.52, 0.14, 0.14]
j = 0
# annotate graph with significance
for i in range(len(Significance)):
    if Significance[i] != '':
        AddSignificanceToBars(ax, Significance[i], XLinePos[j], XLinePos[j+1], YLinePos[i], XStarPos[i], YStarPos[i])        
    j += 2        
        
# save figure
fig.savefig('test.pdf', bbox_inches = 'tight')
