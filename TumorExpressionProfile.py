# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:00:08 2017

@author: RJovelin
"""

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
            # get expression data for that gene in given sample
            Genes[gene][ID] = int(vals[i])
infile.close()
print('recorded expression for each gene')
print(len(Genes))

# create 2 separate dicts for expression in normal and tumor {gene: tissue: [expression]}
NormalTissue, TumorTissue = {}, {}
# loop over genes
for gene in Genes:
    NormalTissue[gene], TumorTissue[gene] = {}, {}
    # loop over aliquot IDs
    for ID in Genes[gene]:
        # determine cancer and disease status and participant ID
        cancer, status = Samples[ID][0], Samples[ID][1]
        if status == 'normal':
            # check if cancer (tissue) is recorded in dict
            if cancer not in NormalTissue[gene]:
                # initialize list
                NormalTissue[gene][cancer] = [Genes[gene][ID]]
            else:
                NormalTissue[gene][cancer].append(Genes[gene][ID])
        elif status == 'tumour':
            # check if cancer (tissue) is recorded in dict
            if cancer not in TumorTissue[gene]:
                # initialize list
                TumorTissue[gene][cancer] = [Genes[gene][ID]]
            else:
                TumorTissue[gene][cancer].append(Genes[gene][ID])
print('recorded expression for cancer: ', 'normal', len(NormalTissue), 'tissue', len(TumorTissue))

# compute the median expression across samples for a given tissue
for gene in NormalTissue:
    for cancer in NormalTissue[gene]:
        NormalTissue[gene][cancer] = np.median(NormalTissue[gene][cancer])
for gene in TumorTissue:
    for cancer in TumorTissue[gene]:
        TumorTissue[gene][cancer] = np.median(TumorTissue[gene][cancer])
print('computed median expression: ', 'normal', len(NormalTissue), 'tissue', len(TumorTissue))

# make a dict with expression profile for each gene for normal and tumor
NormalExp, TumorExp = {}, {}
# loop for genes
for gene in NormalTissue:
    # make a list of tissue
    NormalExp[gene] = []
    for cancer in Tissues:
        NormalExp[gene].append(NormalTissue[gene][cancer])
for gene in TumorTissue:
    # make a list of tissue
    TumorExp[gene] = []
    for cancer in Tissues:
        TumorExp[gene].append(TumorTissue[gene][cancer])
print('recorded expression profile: ', 'normal', len(NormalExp), 'tumor', len(TumorExp))



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

# record only genes with coordinates
to_remove = [gene for gene in NormalExp if gene not in GeneCoord]
for gene in to_remove:
    del NormalExp[gene]
to_remove = [gene for gene in TumorExp if gene not in GeneCoord]
for gene in to_remove:
    del TumorExp[gene]
print('removed genes without coordinates: ', 'normal', len(NormalExp), 'tumor', len(TumorExp))

# remove genes without expression
NormalExp = RemoveGenesLackingExpression(NormalExp)
TumorExp = RemoveGenesLackingExpression(TumorExp)
print('removed genes without expression: ', 'normal', len(NormalExp), 'tumor', len(TumorExp))

# remove genes that do not have expression in both normal and tumor
to_remove = [gene for gene in NormalExp if gene not in TumorExp]
for gene in TumorExp:
    if gene not in NormalExp:
        to_remove.append(gene)
for gene in to_remove:
    if gene in NormalExp:
        del NormalExp[gene]
    if gene in TumorExp:
        del TumorExp[gene]
print('removed genes without paired expression: ', 'normal', len(NormalExp), 'tumor', len(TumorExp))       

# sacle expression data
NormalExp = ScaleExpression(NormalExp, 'level_scaling')
TumorExp = ScaleExpression(TumorExp, 'level_scaling')
print('scaled expression using median: ', 'normal', len(NormalExp), 'tumor', len(TumorExp))

# get relative expression
NormalExp = TransformRelativeExpression(NormalExp)
TumorExp = TransformRelativeExpression(TumorExp)
print('transformed relative expression: ', 'normal', len(NormalExp), 'tumor', len(TumorExp))
        

# load dictionary of nested genes
json_data = open('HumanNestedGenes.json')
Nested = json.load(json_data)
json_data.close()

# make a list of gene pairs
NestedPairs = GetHostNestedPairs(Nested)
print('generated gene pairs', len(NestedPairs))

# remove pairs for which both genes are not expressed
NestedPairs = FilterGenePairsWithoutExpression(NestedPairs, NormalExp, 'strict')
print('filtered gene pairs lacking expression', len(NestedPairs))

# make parallel lists with euclidian distances between external and internal genes for normal and tumor 
ExpDivNormal = ComputeExpressionDivergenceGenePairs(NestedPairs, NormalExp)
ExpDivTumor = ComputeExpressionDivergenceGenePairs(NestedPairs, TumorExp)
print('computed expression divergence')
print('normal', np.mean(ExpDivNormal))
print('tumor', np.mean(ExpDivTumor))


# create figure
fig = plt.figure(1, figsize = (6, 2.5))
# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# plot all gene or non-overlapping gene density first
ax.plot(ExpDivNormal, ExpDivTumor, marker = 'o', markeredgecolor = 'black',
        markeredgewidth = 1, markerfacecolor = 'red',  markersize = 4)      
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# set axis labels
ax.set_ylabel('Expression divergence in tumors', size = 8, color = 'black', ha = 'center', **FigFont)
ax.set_xlabel('Expression divergence in normal', size = 8, color = 'black', ha = 'center', **FigFont )        
# add a range for the Y axis
#plt.ylim([0, YMax])
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False) 
ax.spines["left"].set_visible(True)  
    
# edit tick paramters
plt.tick_params(axis='both', which='both', bottom='on', top='off', 
                right = 'off', left = 'on', labelbottom='on', colors = 'black',
                labelsize = 8, direction = 'out')  
# set x axis ticks
#plt.xticks([], [])
    
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   
          
# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
