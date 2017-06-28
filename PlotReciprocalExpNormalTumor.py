# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:00:08 2017

@author: RJovelin
"""


# use this script to plot the median number of pairs per individuals with reciprocal
# expression in tumor and normal tissue from public TCGA data


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
print('gene coord', len(GeneCoord))

# load dictionary of overlapping genes
json_data = open('HumanOverlappingGenes.json')
Overlap = json.load(json_data)
json_data.close()
# make a list of gene pairs
OverlapPairs = GetHostNestedPairs(Overlap)
print('generated gene pairs', len(OverlapPairs))


# create a dictionary with {sample_ID: [cancer, disease_status]} 
Samples = {}
# create a dictionary with sample_ID: participant_ID pairs
Participants = {}
infile = open('DGE_paired_PCAWG_TCGA_metadata_FLAMAZE_for_methods_Dec_8_2016.tsv')
header = infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get sample ID, cancer type and disease status
        ID, tumor, status = line[0], line[5], line[6]
        # populate dicts        
        assert ID not in Samples
        Samples[ID] = [tumor, status]
        # get participant ID
        participant = line[3]
        assert ID not in Participants
        Participants[ID] = participant
infile.close()
print('paired IDs to cancer and disease status')

# make a list of tumors
Tissues = ['BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA', 'HNSC', 'KICH', 'KIRC',
           'KIRP', 'LIHC', 'LIRI', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'PRAD', 'READ',
           'RECA', 'SARC', 'STAD', 'THCA', 'THYM', 'UCEC']



# make a set with aliquot IDs from TCGA
TCGA = set()
infile = open('DGE_paired_TCGA_meta_1_to_9_raw_count_FLAMAZE_for_methods_Dec_8_2016.tsv')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        TCGA.add(line[4])
infile.close()
print('extracted valid participant IDs from TCGA')


# make a dict with aliquotID: expression pairs for each gene {gene: aliquot: expression}
Genes = {}
infile = open('DGE_paired_PCAWG_TCGA_rowENSGids_colAliquotIds_FLAMAZE_for_methods_Dec_8_2016.tsv')
# make a list of aliquot IDs
aliquots = infile.readline().rstrip().split('\t')
aliquots = aliquots[1:]
# loop over genes in file
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # get gene and expression values
        gene, vals = line[0], line[1:]
        # initialize dict
        assert gene not in Genes
        Genes[gene] = {}
        # loop over the expression values
        for i in range(len(vals)):
            # grab the corresponding aliquot ID
            ID = aliquots[i]
            assert ID not in Genes[gene]
            # check that aliquot ID is in TCGA data
            if ID in TCGA:
                # get expression data for that gene in given sample
                Genes[gene][ID] = int(vals[i])
infile.close()
print('recorded expression for each gene: ', len(Genes))


# remove genes without coordinates
to_remove = [gene for gene in Genes if gene not in GeneCoord]
for gene in to_remove:
    del Genes[gene]
print('removed genes without coordinates: ', len(Genes))


# create 2 separate dicts for expression in normal and tumor {gene: {tissue: {participant: expression}
NormalTissue, TumorTissue = {}, {}
# loop over genes
for gene in Genes:
    NormalTissue[gene], TumorTissue[gene] = {}, {}
    # loop over aliquot IDs
    for ID in Genes[gene]:
        # determine cancer and disease status and participant ID
        cancer, status = Samples[ID][0], Samples[ID][1]
        participant = Participants[ID]
        if status == 'normal':
            # check if cancer (tissue) is recorded in dict
            if cancer not in NormalTissue[gene]:
                # initialize list
                NormalTissue[gene][cancer] = {}
            # populate dict with participant: expression pairs
            assert participant not in NormalTissue[gene][cancer]
            NormalTissue[gene][cancer][participant] = Genes[gene][ID]    
        elif status == 'tumour':
            # check if cancer (tissue) is recorded in dict
            if cancer not in TumorTissue[gene]:
                # initialize dict
                TumorTissue[gene][cancer] = {}
            # populate dict with participant: expression pairs
            assert participant not in TumorTissue[gene][cancer]
            TumorTissue[gene][cancer][participant] = Genes[gene][ID]

print('recorded expression for cancer: ', 'normal', len(NormalTissue), 'tissue', len(TumorTissue))


# make a dict {tissue: {participant: {gene: expression}}}
TissueHealthy, TissueCancer = {}, {}
for gene in NormalTissue:
    for tissue in NormalTissue[gene]:
        # initialize dict if tissue not in dict
        if tissue not in TissueHealthy:
            TissueHealthy[tissue] = {}
        for participant in NormalTissue[gene][tissue]:
            # initialize dict if participant not in dict
            if participant not in TissueHealthy[tissue]:
                TissueHealthy[tissue][participant] = {}
            # populate dict
            assert gene not in TissueHealthy[tissue][participant]
            TissueHealthy[tissue][participant][gene] = NormalTissue[gene][tissue][participant]
for gene in TumorTissue:
    for tissue in TumorTissue[gene]:
        # initialize dict if tissue not in dict
        if tissue not in TissueCancer:
            TissueCancer[tissue] = {}
        for participant in TumorTissue[gene][tissue]:
            # initialize dict if participant not in dict
            if participant not in TissueCancer[tissue]:
                TissueCancer[tissue][participant] = {}
            # populate dict
            assert gene not in TissueCancer[tissue][participant]
            TissueCancer[tissue][participant][gene] = TumorTissue[gene][tissue][participant]
print('recorded expression per participant and tissue')

print('number of healthy and tumor tissues', len(TissueHealthy), len(TissueCancer))
for tissue in TissueHealthy:
    print(tissue, len(TissueHealthy[tissue]), len(TissueCancer[tissue]))


# remove genes without expression
for tissue in TissueHealthy:
    for participant in TissueHealthy[tissue]:
        to_remove = [gene for gene in TissueHealthy[tissue][participant] if TissueHealthy[tissue][participant][gene] == 0]
        for gene in to_remove:
            del TissueHealthy[tissue][participant][gene]
for tissue in TissueCancer:
    for participant in TissueCancer[tissue]:
        to_remove = [gene for gene in TissueCancer[tissue][participant] if TissueCancer[tissue][participant][gene] == 0]
        for gene in to_remove:
            del TissueCancer[tissue][participant][gene]
print('removed genes without expression')

     
# count pairs with reciprocal expression arrangements for normal and tumor
# qet counts for oiverlapping gene pairs
PairCountsOverlap = ReciprocalExpressionTumorNormal(OverlapPairs, Tissues, TissueCancer, TissueHealthy)

# generate distributions of reciprocal expression by resampling random pairs
# make a list with all genes
AllGenes = list(GeneCoord.keys())

# create a list to store the counts of each for each reciprocal arrangement
PairCountsRandom = []
for i in range(9):
    PairCountsRandom.append([])

for i in range(1000):
    # make a list of random gene pairs
    RandomPairs = []
    while len(RandomPairs) != len(OverlapPairs):
        # pick 2 genes at random
        j, k = random.randint(0, len(AllGenes) -1), random.randint(0, len(AllGenes) -1)
        # do not form pairs of the same gene
        if j != k:
            gene1, gene2 = AllGenes[j], AllGenes[k]
            RandomPairs.append([gene1, gene2])
    # get counts for random pairs
    CountsRun = ReciprocalExpressionTumorNormal(RandomPairs, Tissues, TissueCancer, TissueHealthy)
    for m in range(len(CountsRun)):
        PairCountsRandom[m].extend(CountsRun[m])
print('generated random distribution of median reciprocal expression')


# compare reciprocal expressions between overlap and random
Pattern = ['CoExpBoth', 'CoExpCancerDiscordNorm', 'CoExpCancerNotNorm',
            'DiscordCancerCoExpNorm', 'DiscordBoth', 'DiscordCancerNotNorm',
            'NotCancerCoExpNorm', 'NotCancerDiscordNorm', 'NotBoth']        
for i in range(len(PairCountsOverlap)):
    print(Pattern[i], np.median(PairCountsOverlap[i]), np.median(PairCountsRandom[i]),
          np.mean(PairCountsOverlap[i]), np.mean(PairCountsRandom[i]),
          stats.ranksums(PairCountsOverlap[i], PairCountsRandom[i])[1])




# create plotting function
def CreateAx(Columns, Rows, Position, figure, XCoord, YCoord, Data, Title, MaxVal):
    '''
    return an ax in figure
    '''    
        
    # create subplot in figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # get colors    
    #cm = plt.cm.get_cmap('Blues')
    
    # compute median and standard deviation of the data, use as bubble size
    area1 = list(map(lambda x: np.median(x), Data))
    area2 = list(map(lambda x: np.std(x), Data))
    
    # plot data 
    sc1 = ax.scatter(XCoord, YCoord, s = area1, c = area1, cmap = 'Purples', linewidths = 0.7, edgecolor = 'black', alpha = 1, vmin = 0, vmax = MaxVal)      
    sc2 = ax.scatter(XCoord, YCoord, s = area2, c = area2, cmap = 'Greens', linewidths = 0.7, edgecolor = 'black', alpha = 1)      

    # add color scale
    cbar = plt.colorbar(sc1)
    # edit tcik parameters of the scale
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.tick_params(direction = 'out')
    
    # add grid behind plot
    ax.yaxis.grid(color='gray', linestyle='dotted', linewidth = 0.7)
    ax.xaxis.grid(color='gray', linestyle='dotted', linewidth = 0.7)
    ax.set_axisbelow(True)

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # set axis labels
    ax.set_ylabel('Expression in normal', size = 8, color = 'black', ha = 'center', **FigFont)
    ax.set_xlabel('Expression in tumor', size = 8, color = 'black', ha = 'center', **FigFont )        
    # add title
    ax.set_title(Title, size = 8, color = 'black', ha = 'center', **FigFont )     
    # set x axis ticks
    plt.xticks([1, 2, 3], ['Co-expressed', 'Discordant', 'Not expressed'])
    plt.yticks([1, 2, 3], ['Co-expressed', 'Discordant', 'Not expressed'])
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False) 
    ax.spines["left"].set_visible(False)  
    # edit tick paramters
    plt.tick_params(axis='both', which='both', bottom='on', top='off', 
                    right = 'off', left = 'on', labelbottom='on', colors = 'black',
                    labelsize = 8, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    return ax




# create figure
fig = plt.figure(1, figsize = (8.8, 4))

# make a list with x and y coordinates
x = [1, 1, 1, 2, 2, 2, 3, 3, 3]
y = [1, 2, 3, 1, 2, 3, 1, 2, 3]

# compute maximum value for color scale
A, B = list(map(lambda x: np.median(x), PairCountsOverlap)), list(map(lambda x: np.std(x), PairCountsOverlap))
C, D = list(map(lambda x: np.median(x), PairCountsRandom)), list(map(lambda x: np.std(x), PairCountsRandom))
MaxVal = 0
for L in [A, B, C, D]:
    for val in L:
        if val >= MaxVal:
            MaxVal = val
# generate subplots
ax1 = CreateAx(2, 1, 1, fig, x, y, PairCountsOverlap, 'Overlapping pairs', MaxVal)
ax2 = CreateAx(2, 1, 2, fig, x, y, PairCountsRandom, 'Random pairs', MaxVal)

## add subplot labels
ax1.text(-0.9, 3.8, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 8)
ax2.text(-0.9, 3.8, 'B', ha='center', va='center', color = 'black', fontname = 'Arial', size = 8)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
for extension in ['.pdf', '.eps', '.png']:
    fig.savefig('ReciprocalExpNormalTumor' + extension, bbox_inches = 'tight')
