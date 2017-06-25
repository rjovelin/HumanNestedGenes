# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""

# use this script to plot expression divergence between host-nested pairs and their un-nested orthologs


# usage python3 PlotExpDivergYoungOld.py [options]
# [human/mouse]: nested genes in human or mouse


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


# get option from command
FocalSp = sys.argv[1]
assert FocalSp in ['mouse', 'human']
if FocalSp == 'human':
    SisterSp = 'mouse'
elif FocalSp == 'mouse':
    SisterSp = 'human'



if FocalSp == 'human':
    # make a list of GFF files
    GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Mus_musculus.GRCm38.88.gff3',
                 'Pan_troglodytes.CHIMP2.1.4.88.gff3', 'Gorilla_gorilla.gorGor3.1.88.gff3',
                 'Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Pongo_abelii.PPYG2.88.gff3',
                 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3', 'Bos_taurus.UMD3.1.88.gff3',
                 'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
                 'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
                 'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
                 'Monodelphis_domestica.BROADO5.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3',
                 'Sorex_araneus.COMMON_SHREW1.88.gff3']
    # make a parallel list of Species names
    Species = ['Human', 'Mouse', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
               'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
               'Cat', 'Opossum', 'Platypus', 'Shrew']    
    
elif FocalSp == 'mouse':
    # make a list of GFF files
    GFF_Files = ['Mus_musculus.GRCm38.88.gff3', 'Homo_sapiens.GRCh38.88.gff3',
                 'Pan_troglodytes.CHIMP2.1.4.88.gff3', 'Gorilla_gorilla.gorGor3.1.88.gff3',
                 'Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Pongo_abelii.PPYG2.88.gff3',
                 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3', 'Bos_taurus.UMD3.1.88.gff3',
                 'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
                 'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
                 'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
                 'Monodelphis_domestica.BROADO5.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3',
                 'Sorex_araneus.COMMON_SHREW1.88.gff3']
    # make a parallel list of Species names
    Species = ['Mouse', 'Human', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
               'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
               'Cat', 'Opossum', 'Platypus', 'Shrew']    
    
    
# make a parallel list of json files with nested genes
JsonFiles = [i + 'NestedGenes.json' for i in Species]

# make a parallel list of ortholog files
OrthoFiles = [FocalSp.title() + i + 'Orthologs.txt' for i in Species[1:]]


# make a list of dictionaries with nested genes
AllNestedGenes = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    nested = json.load(json_data)
    json_data.close()
    AllNestedGenes.append(nested)
# make a list of dicts with orthologs
AllOrthologs = []
for i in range(len(OrthoFiles)):
    orthologs = MatchOrthologs(OrthoFiles[i])
    AllOrthologs.append(orthologs)


# make a set of overlapping genes in the focal species
json_data = open(FocalSp.title() + 'OverlappingGenes.json')
FocalSpOverlap = json.load(json_data)
json_data.close()
FocalSpOverlapGenes = MakeFullPartialOverlapGeneSet(FocalSpOverlap)
# make a set of overlapping genes in the sisterspecies
json_data = open(SisterSp.title() + 'OverlappingGenes.json')
SisterSpOverlap = json.load(json_data)
json_data.close()
SisterSpOverlapGenes = MakeFullPartialOverlapGeneSet(SisterSpOverlap)
print('generated sets of overlapping genes in {0}  ({1}) and {2} ({3})'.format(FocalSp, len(FocalSpOverlapGenes), SisterSp, len(SisterSpOverlapGenes)))


# generate a set of nested genes in focal and sister species
FocalSpNestedGenes = MakeFullPartialOverlapGeneSet(AllNestedGenes[0])
SisterSpNestedGenes = MakeFullPartialOverlapGeneSet(AllNestedGenes[1])
print('nested genes in {0} and {1}: {2}, {3}: '.format(FocalSp, SisterSp, len(FocalSpNestedGenes), len(SisterSpNestedGenes)))


# get nested pairs in focal sp and other species (sister species is 1st in list)
FocalSpNestedPairs = GetHostNestedPairs(AllNestedGenes[0])
# get nested pairs in other species
SpeciesNestedPairs = [GetHostNestedPairs(AllNestedGenes[i]) for i in range(1, len(AllNestedGenes))]
# infer old and young nested events in focal sp
OldNested, YoungNested = InferYoungOldNestingEvents(FocalSpNestedPairs, SpeciesNestedPairs, AllOrthologs) 
print('inferred young and old nesting events: ', 'old:', len(OldNested), 'new:', len(YoungNested))



# get expression profiles of focal and sister species
if FocalSp == 'human':
    FocalSpExpression = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
    SisterSpExpression = ParseExpressionFile('Mouse_Median_Normalized_FPKM.txt')
    # match expression profiles between mouse and human 
    FocalSpExpression = MatchHumanToMouseExpressionProfiles(FocalSpExpression)
elif FocalSp == 'mouse':
    SisterSpExpression = ParseExpressionFile('GTEX_Median_Normalized_FPKM.txt')
    FocalSpExpression = ParseExpressionFile('Mouse_Median_Normalized_FPKM.txt')
    # match expression profiles between mouse and human 
    SisterSpExpression = MatchHumanToMouseExpressionProfiles(SisterSpExpression)
# remove genes without expression
FocalSpExpression = RemoveGenesLackingExpression(FocalSpExpression)
SisterSpExpression = RemoveGenesLackingExpression(SisterSpExpression)
# scale expression values to the median in each species
FocalSpExpression = ScaleExpression(FocalSpExpression, 'level_scaling')
SisterSpExpression = ScaleExpression(SisterSpExpression, 'level_scaling')
# get relative expression
FocalSpExpression = TransformRelativeExpression(FocalSpExpression)
SisterSpExpression = TransformRelativeExpression(SisterSpExpression)
# compute expression specificity
FocalSpspecificity = ExpressionSpecificity(FocalSpExpression)
SisterSpSpecificity = ExpressionSpecificity(SisterSpExpression)


# make a list of gene coordinates       
AllCoordinates, AllOrdered = [], []
# loop over GFF files
for i in range(len(GFF_Files)):
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF_Files[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF_Files[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    AllCoordinates.append(GeneCoord)
    AllOrdered.append(OrderedGenes)
print('extracted gene coordinates and ordered genes')

# remove orthologs without gene coordinates
SpeciesCoordinates = AllCoordinates[1:]
for i in range(len(SpeciesNestedPairs)):
    to_remove = []
    for pair in SpeciesNestedPairs[i]:
        if pair[0] not in SpeciesCoordinates[i] or pair[1] not in SpeciesCoordinates[i]:
            to_remove.append(pair)
    print('remove {0} pairs without coordinates in {1}'.format(len(to_remove), Species[i+1]))    
    for pair in to_remove:
        SpeciesNestedPairs[i].remove(pair)



# compare expression divergence between gene pairs for young nested pairs in
# human or mouse with their ancestral non-overlapping pairs in mouse or human respectively
    
  
# remove pairs for which both genes are not expressed
YoungNested = FilterGenePairsWithoutExpression(YoungNested, FocalSpExpression, 'strict')
OldNested = FilterGenePairsWithoutExpression(OldNested, FocalSpExpression, 'strict')
print('{0} old and {1} new nested pairs after removing genes lacking expression'.format(len(OldNested), len(YoungNested)))
    
    
# generate gene pairs with orthologs of young nested in sister species
AncestralPairs = []
for pair in YoungNested:
    for ortho1 in AllOrthologs[0][pair[0]]:
        for ortho2 in AllOrthologs[0][pair[1]]:
            AncestralPairs.append([ortho1, ortho2])
            # check that pair not in sister species
            assert [ortho1, ortho2] not in SpeciesNestedPairs[0]
print('generated {0} ancestral pairs in {1}'.format(len(AncestralPairs), SisterSp))


# remove pairs for which both members are not expressed in sister species            
AncestralPairs = FilterGenePairsWithoutExpression(AncestralPairs, SisterSpExpression, 'strict')
print('{0} ancestral pairs after removing pairs without expression'. format(len(AncestralPairs)))    
# remove pairs for which members are not valid genes (eg pseudogenes, etc)        
to_remove = [pair for pair in AncestralPairs if pair[0] not in AllCoordinates[1] or pair[1] not in AllCoordinates[1]]
for pair in to_remove:
    AncestralPairs.remove(pair)        
print('{0} ancestral pairs after removing pairs without coordinates'. format(len(AncestralPairs)))    

# remove pairs if any gene is nested
to_remove = [pair for pair in AncestralPairs if pair[0] in SisterSpNestedGenes or pair[1] in SisterSpNestedGenes]
for pair in to_remove:
    AncestralPairs.remove(pair)
print('{0} ancestral pairs after removing nested genes'.format(len(AncestralPairs)))    
   
    
# pick random pairs if size ancestralapirs > youngnested
L = []
while len(L) != len(YoungNested):
    # pick pair and populate new list
   i = random.randint(0, len(AncestralPairs) -1) 
   pair = AncestralPairs[i]
   L.append(pair)
   # remove pair from list
   AncestralPairs.remove(pair)
# reassign variable name
AncestralPairs = copy.deepcopy(L)
print('{0} ancestral pairs after randomly picking up pairs of orthologs'.format(len(AncestralPairs)))    


# generate a list of control un-nested pairs
SisterSpRandomGenes = GenerateAllUnNestedGenes(SisterSpNestedGenes, AllOrdered[1], SisterSpExpression)
# make a list of control un-nested pairs in sister species
SisterSpControlPairs = []
for pair in AncestralPairs:
    # make a list of matching gene pairs (orientation, chromosome, distance)
    PairPool = GenerateMatchingPoolPairs(pair, SisterSpRandomGenes, AllCoordinates[1], 1000)
    # draw a matching gene pair at random
    i = random.randint(0, len(PairPool) -1)
    SisterSpControlPairs.append(PairPool[i])
    
    
# generate a list of control genes in focal species
FocalSpRandomGenes = GenerateAllUnNestedGenes(FocalSpNestedGenes, AllOrdered[0], FocalSpExpression)
# make a list of control un-nested (non-overlapping) pars in focal species
FocalSpControlPairs = []
for pair in YoungNested:
    # make a list of matching gene pairs (orientation, chromosome, distance)
    PairPool = GenerateMatchingPoolPairs(pair, FocalSpRandomGenes, AllCoordinates[0], 1000)
    # draw a matching gene pair at random
    i = random.randint(0, len(PairPool)-1)
    FocalSpControlPairs.append(PairPool[i])
    assert PairPool[i][0] not in FocalSpNestedGenes
    assert PairPool[i][1] not in FocalSpNestedGenes
print('generated control pairs')    
    
# compute divergence in young nested pairs
FocalSpYoungDiv = ComputeExpressionDivergenceGenePairs(YoungNested, FocalSpExpression)    
# compute divergence in ancestral un-nested pairs
SisterSpAncestralDiv = ComputeExpressionDivergenceGenePairs(AncestralPairs, SisterSpExpression)
# compute expression divergence between un-nested control pairs in sister species    
SisterSpControlDiv = ComputeExpressionDivergenceGenePairs(SisterSpControlPairs, SisterSpExpression)
# compute expression divergence between un-nested control pairs in focal species
FocalSpControlDiv = ComputeExpressionDivergenceGenePairs(FocalSpControlPairs, FocalSpExpression)


# compute P values using permutation tests
P = PermutationResampling(FocalSpYoungDiv, SisterSpAncestralDiv, 1000, statistic = np.mean)
# add P to list
PValues = [P]
# convert P to star significance
PValues = ConvertPToStars(PValues)[0]
print('young vs ancestral', len(FocalSpYoungDiv), len(SisterSpAncestralDiv), np.mean(FocalSpYoungDiv), np.mean(SisterSpAncestralDiv), P)
P = PermutationResampling(SisterSpAncestralDiv, SisterSpControlDiv , 1000, statistic = np.mean)
print('ancestral vs control', len(SisterSpAncestralDiv), len(SisterSpControlDiv), np.mean(SisterSpAncestralDiv), np.mean(SisterSpControlDiv), P)
P = PermutationResampling(FocalSpYoungDiv, FocalSpControlDiv, 1000, statistic = np.mean)
print('young vs contol', len(FocalSpYoungDiv), len(FocalSpControlDiv), np.mean(FocalSpYoungDiv), np.mean(FocalSpControlDiv), P)


# plot results to file

# make a list of gene categories
GeneCatOrientation = ['Derived nested in ' + FocalSp, 'Ancestral non-nested in ' + SisterSp]
# create lists with means and SEM for each gene category
Means, SEM = [], []
for L in [FocalSpYoungDiv, SisterSpAncestralDiv]:
    Means.append(np.mean(L))
    SEM.append(np.std(L) / math.sqrt(len(L)))

# create figure
fig = plt.figure(1, figsize = (1.5, 2))

# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# set colors
colorscheme = ['#225ea8', '#e31a1c']
# plot nucleotide divergence
ax.bar([0.05, 0.35], Means, 0.2, yerr = SEM, color = colorscheme,
        edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.15, 0.45], GeneCatOrientation, size = 7, color = 'black', ha = 'center', **FigFont)
# add title
ax.set_xlabel('Orientation', color = 'black', size = 7, ha = 'center', **FigFont)    
# add a range for the Y and X axes
MaxVal = Means[0] + SEM[0] + 0.2
print(Means[0] + SEM[0], MaxVal)

#plt.ylim([0, 1.5])
plt.xlim([0, 0.6])
# edit y axis ticks
plt.yticks(np.arange(0, 1.5, 0.1)) 
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


# annotate figure to add significance
if PValues != '':
    if FocalSp == 'human':
        ax = AddSignificanceToBars(ax, PValues, 0.15, 0.45, 1.5, 0.3, 1.6)
    elif FocalSp == 'mouse':
        ax = AddSignificanceToBars(ax, PValues, 0.15, 0.45, 1.1, 0.3, 1.2)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
  
## save figure
#for extension in ['.pdf', '.eps', '.png']:
#    fig.savefig('truc' + extension, bbox_inches = 'tight')



