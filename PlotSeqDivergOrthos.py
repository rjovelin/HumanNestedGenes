# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 11:58:03 2017

@author: Richard
"""

# use this script to plot sequence divergence between orthologs 

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

# generate gene sets
NestedGenes  = MakeFullPartialOverlapGeneSet(Nested)
OverlappingGenes = MakeFullPartialOverlapGeneSet(Overlapping)
ConvergentGenes = MakeFullPartialOverlapGeneSet(Convergent)
DivergentGenes = MakeFullPartialOverlapGeneSet(Divergent)
PiggyBackGenes = MakeFullPartialOverlapGeneSet(Piggyback)

# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)

# create sets of internal and external nested gene pairs
NestedPairs = GetHostNestedPairs(Nested)
InternalGenes, ExternalGenes = set(), set()
for pair in NestedPairs:
    ExternalGenes.add(pair[0])
    InternalGenes.add(pair[1])
    
## get 1:1 orthologs between human and chimp
#Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')


Orthos = {}
files = [i for i in os.listdir('HumanChimpCodeml/') if '.txt' in i and 'ENSG' in i and 'ENSPTR' in i and 'codeml' not in i and '.out' not in i]
# loop over files, get the aligned coding sequences
for filename in files:
    genes = filename[:filename.index('.txt')]
    genes = genes.split('_')
    Orthos[genes[0]] = genes[1]



# create lists of orthologous pairs for each gene category 
GeneCats = ['NoOv', 'Nst', 'Int', 'Ext', 'Pbk', 'Conv', 'Div']
AllPairs = []
AllGenes = [NonOverlappingGenes, NestedGenes, InternalGenes, ExternalGenes,
            PiggyBackGenes, ConvergentGenes, DivergentGenes] 
# loop over gene sets
for i in range(len(AllGenes)):
    # create a list of gene pairs
    orthologs = []
    # loop over genes in given gene set
    for gene in AllGenes[i]:
        # check that gene has ortholog
        if gene in Orthos:
            orthologs.append([gene, Orthos[gene]])
    AllPairs.append(orthologs)

    
# compute AA differences between orthologs
Proteins = {}
# make a list of codon based alignment files
files = [i for i in os.list('HumanChimpCodeml/') if '.txt' in i and 'ENSG' in i and 'ENSPTR' in i and 'codeml' not in i and '.out' not in i]
# loop over files, get the aligned coding sequences
for filename in files:
    infile = open('HumanChimpCodeml/' + filename)
    # skip the first line
    infile.readline()
    content = infile.read().rstrip()
    infile.close()
    content = content.split('>')
    content.remove('')
    for i in range(len(content)):
        content[i] = content[i].split('\n')
    hsagene, hsaseq = content[0][0], content[0][1]
    ptrgene, ptrseq = content[1][0], content[1][1]
    assert 'ENSG' in hsagene and 'ENSPTR' in ptrgene
    assert hsagene not in Proteins and ptrgene not in Proteins
    Proteins[hsagene] = hsaseq
    Proteins[ptrgene] = ptrseq
# translate sequences
for gene in Proteins:
    Proteins[gene] = TranslateCDS(Proteins[gene])



weird = 0
    
# loop over ortholog pairs for each gene category and compute the protein distance
AADiffs = []
for i in range(len(AllPairs)):
    protdiff = []
    # loop of pairs in given list
    for pair in AllPairs[i]:
        if pair[0] in Proteins:
            D = ProteinDistance(Proteins(pair[0]), Proteins(pair[1]))
            if D != 'NA':
                protdiff.append(D)
            else:
                weird += 1
    AADiffs.append(protdiff)
print('undefined prot', weird)    
    
# create a dict with divergence values
SeqDiv = {}
infile = open('HumanChimpSeqDiverg.txt')
infile.readline()
for line in infile:
    if line.startswith('ENSG'):
        line = line.rstrip().split('\t')
        SeqDiv[line[0]] = [float(line[2]), float(line[3])]
        if line[-1] != 'NA':
            SeqDiv[line[0]].append(float(line[-1]))
        else:
            SeqDiv[line[0]].append(line[-1])
infile.close()            
    
# make lists of dN, dS and omega values for each gene category
dN, dS, Omega = [], [], []
for i in range(len(AllPairs)):
    if pair[0] in SeqDiv:
        nonsyn, syn, ratio = [], [], []
        for pair in AllPairs[i]:
             nonsyn.append(SeqDiv[pair[0]][0])
             syn.append(SeqDiv[pair[0]][1])
             if SeqDiv[pair[0]][-1] != 'NA':
                 ratio.append(SeqDiv[pair[0]][-1])
    dN.append(nonsyn)
    dS.append(syn)
    Omega.append(ratio)
    
# create lists with means and SEM for each gene category
MeanAADiff, SEMAADiff = [], []
for i in range(len(AADiffs)):
    MeanExpDiv.append(np.mean(AADiffs[i]))
    SEMExpDiv.append(np.std(AADiffs[i]) / math.sqrt(len(AADiffs[i])))

 # create lists with means and SEM for dN for each gene category
MeandN, SEMdN
for i in range(len(dN)):
    MeandN.append(np.mean(dN[i]))
    SEMdN.append(np.std(dN[i]) / math.sqrt(len(dN[i])))
    
for i in range(len(GeneCats)):
    print(i, Genecats[i], MeanAADiff[i], MeandN[i], sep = '\t')




   
## create figure
#fig = plt.figure(1, figsize = (2.1, 2))
#
## add a plot to figure (N row, N column, plot N)
#ax = fig.add_subplot(1, 1, 1)
## set colors
#colorscheme = ['black','lightgrey','lightgrey','lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
## plot nucleotide divergence
#ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
#       edgecolor = 'black', linewidth = 0.7,
#       error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
## set font for all text in figure
#FigFont = {'fontname':'Arial'}   
## write y axis label
#ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
## add ticks and lebels
#plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
## add a range for the Y and X axes
#plt.ylim([0, 0.305])
#plt.xlim([0, 2.10])
#
## add margins
#plt.margins(0.1)
#
## do not show lines around figure  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(True)  
## edit tick parameters    
#plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                right = 'off', left = 'on', labelbottom='on',
#                colors = 'black', labelsize = 7, direction = 'out')  
## Set the tick labels font name
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')   
#      
#
## perform statistical tests between gene categories using Kolmogorov-Smirnof test
## create list to store the p-values
#PValues = []
#for i in range(1, len(ExpressionDivergence)):
#    # compare each gene category to non-overlapping genes
#    val, P = stats.ks_2samp(ExpressionDivergence[0], ExpressionDivergence[i])
#    PValues.append(P)
#
## convert p-values to star significance level
#Significance = []
#for pvalue in PValues:
#    if pvalue >= 0.05:
#        Significance.append('')
#    elif pvalue < 0.05 and pvalue >= 0.01:
#        Significance.append('*')
#    elif pvalue < 0.01 and pvalue >= 0.001:
#        Significance.append('**')
#    elif pvalue < 0.001:
#        Significance.append('***')
#
## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#ypos = [0.28] * 4 + [0.22] * 2
#xpos = [0.45, 0.75, 1.05, 1.35, 1.65, 1.95]
#for i in range(len(Significance)):
#    ax.text(xpos[i], ypos[i], Significance[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')
