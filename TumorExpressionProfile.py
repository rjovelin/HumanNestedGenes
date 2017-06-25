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

# load dictionary of nested genes
#json_data = open('HumanNestedGenes.json')


json_data = open('HumanOverlappingGenes.json')


Nested = json.load(json_data)
json_data.close()
# make a list of gene pairs
NestedPairs = GetHostNestedPairs(Nested)
print('generated gene pairs', len(NestedPairs))


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
print('recorded expression for each gene: ', len(Genes))


# remove genes without coordinates
to_remove = [gene for gene in Genes if gene not in GeneCoord]
for gene in to_remove:
    del Genes[gene]
print('removed genes without coordinates: ', len(Genes))


## remove pairs for which genes have missing expression
#NestedPairs = FilterGenePairsWithoutExpression(NestedPairs, Genes, 'strict')
#print('deleted pairs lacking expression', len(NestedPairs))



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


print(len(TissueHealthy), len(TissueCancer))
for tissue in TissueHealthy:
    print(tissue, len(TissueHealthy[tissue]), len(TissueCancer[tissue]))


# remove genes without expression
for tissue in TissueHealthy:
    for participant in TissueHealthy[tissue]:
        to_remove = [gene for gene in TissueHealthy[tissue][participant] if TissueHealthy[tissue][participant][gene] in [0, 1]]
        for gene in to_remove:
            del TissueHealthy[tissue][participant][gene]
for tissue in TissueCancer:
    for participant in TissueCancer[tissue]:
        to_remove = [gene for gene in TissueCancer[tissue][participant] if TissueCancer[tissue][participant][gene] in [0,1]]
        for gene in to_remove:
            del TissueCancer[tissue][participant][gene]
print('removed genes without expression')


# count the number of pairs in which genes are coexpressed, discordantly expressed or not expressed
CoExpBoth, CoExpCancerDiscordNorm, CoExpCancerNotNorm = 0, 0, 0
DiscordCancerCoExpNorm, DiscordBoth, DiscordCancerNotNorm = 0, 0, 0
NotCancerCoExpNorm, NotCancerDiscordNorm, NotBoth = 0, 0, 0

          
#a, b, c, d, e, f, g, h, i = [], [], [], [], [], [], [], [], []



total = 0

# loop over tissue
for tissue in TissueCancer:
    if tissue in Tissues:
        # loop over participant
        for participant in TissueCancer[tissue]:
#            CoExpBoth, CoExpCancerDiscordNorm, CoExpCancerNotNorm = 0, 0, 0
#            DiscordCancerCoExpNorm, DiscordBoth, DiscordCancerNotNorm = 0, 0, 0
#            NotCancerCoExpNorm, NotCancerDiscordNorm, NotBoth = 0, 0, 0
            for pair in NestedPairs:
                total += 1
                # check expression in cancer and normal
                if pair[0] in TissueCancer[tissue][participant] and pair[1] in TissueCancer[tissue][participant]:
                    # gene coexpressed in cancer
                    if pair[0] in TissueHealthy[tissue][participant] and pair[1] in TissueHealthy[tissue][participant]:
                        # gene coexpressed in normal
                        CoExpBoth += 1
                    elif (pair[0] in TissueHealthy[tissue][participant] and pair[1] not in TissueHealthy[tissue][participant]) or (pair[0] not in TissueHealthy[tissue][participant] and pair[1] in TissueHealthy[tissue][participant]):
                        # gene discordant in normal
                        CoExpCancerDiscordNorm += 1
                    elif pair[0] not in TissueHealthy[tissue][participant] and pair[1] not in TissueHealthy[tissue][participant]:
                        # gene not expressed in normal
                        CoExpCancerNotNorm += 1
                elif (pair[0] in TissueCancer[tissue][participant] and pair[1] not in TissueCancer[tissue][participant]) or (pair[0] not in TissueCancer[tissue][participant] and pair[1] in TissueCancer[tissue][participant]):
                    # gene discordant in cancer
                    if pair[0] in TissueHealthy[tissue][participant] and pair[1] in TissueHealthy[tissue][participant]:
                        # gene coexpressed in normal
                        DiscordCancerCoExpNorm += 1
                    elif (pair[0] in TissueHealthy[tissue][participant] and pair[1] not in TissueHealthy[tissue][participant]) or (pair[0] not in TissueHealthy[tissue][participant] and pair[1] in TissueHealthy[tissue][participant]):
                        # gene discordant in normal
                        DiscordBoth += 1
                    elif pair[0] not in TissueHealthy[tissue][participant] and pair[1] not in TissueHealthy[tissue][participant]:
                        # gene not expressed in normal
                        DiscordCancerNotNorm += 1
                elif pair[0] not in TissueCancer[tissue][participant] and pair[1] not in TissueCancer[tissue][participant]:
                    # gene not expressed in cancer
                    if pair[0] in TissueHealthy[tissue][participant] and pair[1] in TissueHealthy[tissue][participant]:
                        # gene coexpressed in normal
                        NotCancerCoExpNorm += 1
                    elif (pair[0] in TissueHealthy[tissue][participant] and pair[1] not in TissueHealthy[tissue][participant]) or (pair[0] not in TissueHealthy[tissue][participant] and pair[1] in TissueHealthy[tissue][participant]):
                        # gene discordant in normal
                        NotCancerDiscordNorm += 1
                    elif pair[0] not in TissueHealthy[tissue][participant] and pair[1] not in TissueHealthy[tissue][participant]:
                        # gene not expressed in normal
                        NotBoth += 1
#            a.append(CoExpBoth)
#            b.append(CoExpCancerDiscordNorm)
#            c.append(CoExpCancerNotNorm)
#            d.append(DiscordCancerCoExpNorm)
#            e.append(DiscordBoth)
#            f.append(DiscordCancerNotNorm)
#            g.append(NotCancerCoExpNorm)
#            h.append(NotCancerDiscordNorm)
#            i.append(NotBoth)        

















       
print('CoExpBoth', CoExpBoth)
print('CoExpCancerDiscordNorm', CoExpCancerDiscordNorm)
print('CoExpCancerNotNorm', CoExpCancerNotNorm)
print('DiscordCancerCoExpNorm', DiscordCancerCoExpNorm)
print('DiscordBoth', DiscordBoth)
print('DiscordCancerNotNorm', DiscordCancerNotNorm)
print('NotCancerCoExpNorm', NotCancerCoExpNorm) 
print('NotCancerDiscordNorm', NotCancerDiscordNorm)
print('NotBoth', NotBoth)

            


x = [1, 1, 1, 2, 2, 2, 3, 3, 3]
y = [1, 2, 3, 1, 2, 3, 1, 2, 3]
area = [CoExpBoth, CoExpCancerDiscordNorm, CoExpCancerNotNorm, DiscordCancerCoExpNorm, DiscordBoth, DiscordCancerNotNorm, NotCancerCoExpNorm, NotCancerDiscordNorm, NotBoth]

assert sum(area) == total


area = list(map(lambda x: x / 1000, area))

#area = list(map(lambda x: np.median(x) * 10, [a, b, c, d, e, f, g, h, i]))




##area = np.pi * (15 * np.random.rand(N))**2 # 0 to 15 point radiuses

#cm = plt.cm.get_cmap('jet')
#sc = ax.scatter(x,y,s=z*500,c=z,cmap=cm,linewidth=0,alpha=0.5)
#ax.grid()
#fig.colorbar(sc)

# color.append(data[6]) # larceny_theft 


#area = list(map(lambda x: math.sqrt(x) * 10, area))



#color = ['#2b8cbe'] * len(area)

#color = area



cm = plt.cm.get_cmap('Blues')


# create figure
fig = plt.figure(1, figsize = (5, 5))
# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# plot all gene or non-overlapping gene density first
#sc = ax.scatter(x, y, c = color, s = area, linewidths = 0.7, edgecolor = 'black', alpha = 1)      


sc = ax.scatter(x, y, s = area, c = area, cmap = 'Blues', linewidths = 0.7, edgecolor = 'black', alpha = 1)      



# add color scale
plt.colorbar(sc)


ax.yaxis.grid(color='gray', linestyle='dotted', linewidth = 0.7)
ax.xaxis.grid(color='gray', linestyle='dotted', linewidth = 0.7)
ax.set_axisbelow(True)


# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# set axis labels
ax.set_ylabel('Expression in normal', size = 8, color = 'black', ha = 'center', **FigFont)
ax.set_xlabel('Expression in tumor', size = 8, color = 'black', ha = 'center', **FigFont )        
# set x axis ticks
plt.xticks([1, 2, 3], ['Co-expressed', 'Discordant', 'Not expressed'])
plt.yticks([1, 2, 3], ['Co-expressed', 'Discordant', 'Not expressed'])
# add a range for the Y axis
#plt.ylim([0, YMax])
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
# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')


