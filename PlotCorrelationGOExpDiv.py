# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:51:38 2017

@author: RJovelin
"""

# use this script to plot the correlation between expression divergence between gene pairs and functional similarity

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
 
# generate random pairs of non-overlapping genes as a control
# make a list of non-overlapping genes convert set of non-overlapping genes to list
NonOverlappingGenes = list(NonOverlappingGenes)
ControlPairs = []
# draw 10000 random pairs
replicates = 10000
while replicates != 0:
    # pick 2 random genes
    j = random.randint(0, len(NonOverlappingGenes) -1)
    k = random.randint(0, len(NonOverlappingGenes) -1)
    gene1, gene2 = NonOverlappingGenes[j], NonOverlappingGenes[k]
    if gene1 in GeneOntology and gene2 in GeneOntology:
        JI = JaccardIndex(GeneOntology[gene1], GeneOntology[gene2])
        ControlPairs.append(JI)
        replicates -= 1
  
# insert pairs of JI in list
FunctionalSimilarity.insert(0, ControlPairs)

DataSets = ['Control', 'Nested', 'PiggyBack', 'Convergent', 'Divergent']

for i in range(1, len(FunctionalSimilarity)):
    P = PermutationResampling(FunctionalSimilarity[0], FunctionalSimilarity[i], 1000, statistic = np.mean)
    print(DataSets[0], DataSets[i], len(FunctionalSimilarity[0]), len(FunctionalSimilarity[i]), np.mean(FunctionalSimilarity[0]), np.mean(FunctionalSimilarity[i]), P)


###############################

    

# make sets of host and nested nested genes
OverlapSets = MakeFullPartialOverlapGeneSet(Overlap[0])
HumanExpression = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
# compute expression specificity
HumanSpecificity = ExpressionSpecificity(HumanExpression)
# generate a dict to draw genes in human    
HumanRandomGenes = GenerateAllUnNestedGenes(OverlapSets, OrderedGenes, HumanExpression)
    
    
NestedPairs = copy.deepcopy(OverlappingPairs[1])
# remove human pairs if genes are not expressed
to_remove = [pair for pair in NestedPairs if pair[0] not in HumanExpression or pair[1] not in HumanExpression]
for pair in to_remove:
    NestedPairs.remove(pair)
print('expressed nested', len(NestedPairs), len(OverlappingPairs[1]))
  



  
# make a list of control un-nested pairs in sister species
HumanControlPairs = []    
for pair in NestedPairs:
    # make a list of matching gene pairs (orientation, chromosome, distance)
    PairPool = GenerateMatchingPoolPairs(pair, HumanRandomGenes, GeneCoord, 2000)
    # draw a matching gene pair at random
    i = random.randint(0, len(PairPool) -1)
    assert PairPool[i][0] not in OverlapSets and PairPool[i][1] not in OverlapSets
    HumanControlPairs.append(PairPool[i])


NestedJI = []
for pair in NestedPairs:
    if pair[0] in GeneOntology and pair[1] in GeneOntology:
        JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
        NestedJI.append(JI)
ControlJI = []
for pair in HumanControlPairs:
    if pair[0] in GeneOntology and pair[1] in GeneOntology:
        JI = JaccardIndex(GeneOntology[pair[0]], GeneOntology[pair[1]])
        ControlJI.append(JI)    

P = PermutationResampling(NestedJI, ControlJI, 1000, statistic = np.mean)
print(len(NestedJI), len(ControlJI), np.mean(NestedJI), np.mean(ControlJI), P)








######################






# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 18:38:17 2016

@author: Richard
"""

# use this script to plot expression divergence between gene pairs as a function
# of the distance between the 2 genes



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


# create lists of gene pairs
OverlappingPairs = GetHostNestedPairs(Overlapping)
NestedPairs = GetHostNestedPairs(Nested)
PiggybackPairs = GetHostNestedPairs(Piggyback)
ConvergentPairs = GetHostNestedPairs(Convergent)
DivergentPairs = GetHostNestedPairs(Divergent)


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

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# generate lists of gene pairs separated by distance 
Proximal, Moderate, Intermediate, Distant = GenerateSetsGenePairsDistance(GeneCoord, OrderedGenes, ExpressionProfile)

# filter gene pairs lacking expression
AllPairs = [NestedPairs, PiggybackPairs, ConvergentPairs, DivergentPairs]
for i in range(len(AllPairs)):
    AllPairs[i] = FilterGenePairsWithoutExpression(AllPairs[i], ExpressionProfile)

# add gene pairs to Allpairs list
AllPairs.append(Proximal)
AllPairs.append(Moderate)
AllPairs.append(Intermediate)
AllPairs.append(Distant)

# compute expression divergence between pairs of genes
Divergence = []
for i in range(len(AllPairs)):
    Div = ComputeExpressionDivergenceGenePairs(AllPairs[i], ExpressionProfile)
    Divergence.append(Div)

# make a list of gene category names parallel to the list of gene pairs
GeneCats = ['Nst', 'Pbk', 'Conv', 'Div', 'Prox', 'Mod', 'Int', 'Dist']

# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(Divergence)):
    MeanExpDiv.append(np.mean(Divergence[i]))
    SEMExpDiv.append(np.std(Divergence[i]) / math.sqrt(len(Divergence[i])))

# create figure
fig = plt.figure(1, figsize = (3, 2))

# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# set colors
colorscheme = ['grey','grey','grey','grey', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
# plot nucleotide divergence
ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
       edgecolor = 'black', linewidth = 0.5,
       error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
# add a range for the Y and X axes
plt.ylim([0, 0.61])
plt.xlim([0, 2.45])
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
      

# perform statistical tests between gene categories
# create list to store the p-values
PValues = []
# loop over inner list, compare gene categories
for i in range(0, len(Divergence) -1):
    for j in range(i+1, len(Divergence)):
        P = stats.ranksums(Divergence[i], Divergence[j])[1]
        PValues.append(P)
# print p values
for p in PValues:
    print(p)

# convert p-values to star significance level
Significance = []
for pvalue in PValues:
    if pvalue >= 0.05:
        Significance.append('')
    elif pvalue < 0.05 and pvalue >= 0.01:
        Significance.append('*')
    elif pvalue < 0.01 and pvalue >= 0.001:
        Significance.append('**')
    elif pvalue < 0.001:
        Significance.append('***')


# annotate figure to add significance
# significant comparisons were already determined, add letters to show significance
Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
ypos = [0.55] * 8
xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
for i in range(len(Diff)):
    ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
    
# save figure
fig.savefig('ExpressionDivergenceDistance.pdf', bbox_inches = 'tight')
fig.savefig('ExpressionDivergenceDistance.eps', bbox_inches = 'tight')
