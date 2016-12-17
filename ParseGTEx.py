# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 18:37:02 2016

@author: RJovelin
"""


# use this script to write a file with median normalized FPKM in each tissue for each gene

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


# create a dictionary with sample_ID: tissue type pairs 
Samples = {}
infile = open('GTEx_Data_V6_Annotations_SampleAttributesDS.txt')
header = infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get tissue type (code SMTS)
        tissue = line[5]
        # get sample ID, replace '-' with '.'
        ID = line[0].replace('-', '.')
        assert ID not in Samples
        Samples[ID] = tissue
infile.close()


# create a dictionary with tissue type: [sample ids] 
Tissues = {}
infile = open('GTEx_Data_V6_Annotations_SampleAttributesDS.txt')
header = infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get tissue type (code SMTS)
        tissue = line[5]
        # get sample ID, replace '-' with '.'
        ID = line[0].replace('-', '.')
        if tissue not in Tissues:
            Tissues[tissue] = [ID]
        else:
            Tissues[tissue].append(ID)
infile.close()

# print sample size for each tissue
for i in Tissues:
    print(i, len(Tissues[i]))

# create a dict with sample_index: sample_id in normalized FPKM file
infile = open('Matrix.gtex.rnaseq.fpkmuq.v6p.txt')
header = infile.readline().rstrip().split('\t')
# do not record the gene ID
header = header[1:]
SampleIndex = {}
for i in range(len(header)):
    # QC check that all indices have a corresponding tissue
    assert header[i] in Samples
    SampleIndex[i] = header[i]
# loop over file, collect expression information for each gene
Genes = {}
for line in infile:
    if line.startswith('ENS'):
        line = line.rstrip().split('\t')
        # get gene ID
        gene = line[0]
        # remove extra information that is not part of the valid ensembl gene ID
        if '.' in gene:
            gene = gene[:gene.index('.')]
        # check that gene is not already recorded (ie, if information specifies transcript)
        assert gene not in Genes
        # record only genes present in GeneCoord
        if gene in GeneCoord:
            # record expression information
            Genes[gene] = line[1:]
            # check that length of header and expression match
            assert len(header) == len(Genes[gene])
infile.close()


# create a dict {gene: {tissue:[expression values]}}
Expression = {}
# loop over each gene
for gene in Genes:
    # initialize outer dict
    Expression[gene] = {}
    # loop over indices of expression list
    for i in range(len(Genes[gene])):
        # get the corresponding ID
        ID = SampleIndex[i]
        # get the corresponding tissue
        tissue = Samples[ID]
        # check if tissue is key in dict
        if tissue not in Expression[gene]:
            # initialize list value
            Expression[gene][tissue] = []
        # populate list with expression values
        Expression[gene][tissue].append(float(Genes[gene][i]))

# QC that all genes have the same expression profile (ie, same tissue list)
# make a list of expressed genes
ExpressedGenes = list(Expression.keys())
# get the list o tissues for one of the expressed genes
TissueExpression = list(Expression[ExpressedGenes[0]].keys())
TissueExpression.sort()
for gene in Expression:
    a = list(Expression[gene].keys())
    a.sort()
    assert a == TissueExpression
# make a list of expressed tissues not in tissue types from metada
MissingInTissueExp = [i for i in Tissues if i not in TissueExpression]
print('missing in expressed tissues', MissingInTissueExp)
# make a list of metadata tissues not in expressed tissues
MissingInTissueType = [i for i in TissueExpression if i not in Tissues]
print('missing in matadata tissues', MissingInTissueType)

# an empty string has 10 sample IDs in metadata
# remove empty strings from expression dict
for gene in Expression:
    assert '' in Expression[gene]
    del Expression[gene]['']
# remove empty string from list of expressed tissues
TissueExpression.remove('')
TissueExpression.sort()

# create a dict with expression profile for each gene
# expression profile is the median expression of given gene in each of the tissues 
# in the same order as the list of tissue types
# create a dict {gene: [median_tissue1, median_tissue2]
MedianExp = {}
# loop over genes with expression
for gene in Expression:
    # initialize list value
    MedianExp[gene] = []
    # loop over tissue expression
    for tissue in TissueExpression:
        # take the median expression
        MedianExp[gene].append(np.median(Expression[gene][tissue]))
print(len(MedianExp))


# remove genes with no expression in any tissues
to_remove = []
for gene in MedianExp:
    if sum(MedianExp[gene]) == 0:
        to_remove.append(gene)
for gene in to_remove:
    del MedianExp[gene]
print(len(MedianExp))

# create a dict with sample size per tissue
SampleCounts = {}
# loop over gene with expression
for gene in Expression:
    # loop over tissue
    for tissue in Expression[gene]:
        if tissue not in SampleCounts:
            SampleCounts[tissue] = len(Expression[gene][tissue])
        else:
            assert SampleCounts[tissue] == len(Expression[gene][tissue])
# print expressed tissues and sample sizes
for tissue in SampleCounts:
    print(tissue, SampleCounts[tissue])


# write median normalized fpkm to file
newfile = open('GTEX_Median_Normalized_FPKM.txt', 'w')
header = ['Gene'] + TissueExpression
# write header to file
newfile.write('\t'.join(header) + '\n')
# write gene and expression profile
for gene in MedianExp:
    profile = list(map(lambda x: str(x), MedianExp[gene]))
    line = [gene] + [profile]
    newfile.write('\t'.join(line) + '\n')
newfile.close()
 

