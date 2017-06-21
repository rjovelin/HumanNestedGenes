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

## get GFF file
#GFF = 'Homo_sapiens.GRCh38.88.gff3'
## get the coordinates of genes on each chromo
## {chromo: {gene:[chromosome, start, end, sense]}}
#GeneChromoCoord = ChromoGenesCoord(GFF)
## map each gene to its mRNA transcripts
#MapGeneTranscript = GeneToTranscripts(GFF)
## remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
## get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)


# create a dictionary with sample_ID: cancer pair 
Samples = {}
infile = open('DGE_paired_PCAWG_TCGA_metadata_FLAMAZE_for_methods_Dec_8_2016.tsv')
header = infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get sample ID, cancer type and disease status
        ID, tumor, status = line[0], line[5], line[6]
        assert ID not in Samples
        Samples[ID] = [tumor, status]
infile.close()

# create a dictionary with tumor type {tumor: {ID: status}}
Tumors = {}
infile = open('DGE_paired_PCAWG_TCGA_metadata_FLAMAZE_for_methods_Dec_8_2016.tsv')
header = infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get sample ID, cancer type and disease status
        ID, tumor, status = line[0], line[5], line[6]
        if tumor not in Tumors:
            Tumors[tumor] = {}
        else:
            Tumors[tumor][ID] = status
infile.close()

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


# create a dictionary with tissue as key and a list of dictionaries with gene expression as value
Expression = {}



#XXXXXXXXXXXXXXXXXXXX
#
#
## loop over tissues, compute expression for each sample, store the dictionary
#for tissue in Tissues:
#    L = ComputeTissueExpression(tissue)
#    Expression[tissue] = L
#    
#		
## create a dict with the tissue names matching human tissues and list with mouse gene expression
#MouseExpression = {}
#for tissue in HumanMatches:
#    if len(HumanMatches[tissue]) == 1:
#        # map list of mouse expression to human tissue name
#        MouseExpression[tissue] = Expression[HumanMatches[tissue][0]]
#    else:
#        # merge lists of dictionaries in a single list
#        MouseExpression[tissue] = []
#        for subtissue in HumanMatches[tissue]:
#            MouseExpression[tissue].extend(Expression[subtissue])
#
## Merge all samples for a given tissue
#for tissue in MouseExpression:
#    MouseExpression[tissue] = MergeSamples(MouseExpression[tissue])
#
## take the median expression for a given gene across the different samples
#for tissue in MouseExpression:
#    for gene in MouseExpression[tissue]:
#        MouseExpression[tissue][gene] = np.median(MouseExpression[tissue][gene])
#
## make a list with human tissue names in the same order as the tissues in the human expression file
#HumanTissueNames = ['Adipose Tissue', 'Adrenal Gland', 'Bladder',
#                    'Brain', 'Breast', 'Colon', 'Heart', 'Kidney', 
#                    'Liver', 'Lung', 'Ovary', 'Pancreas', 'Small Intestine',
#                    'Spleen', 'Stomach', 'Testis']	
#
## create a dict to store the expression profile of gene 
## (ie, list of expression values for each tissue in the same order as the list of HumanTissueNames)
#GeneExpression = {}
## initiate dict with empty list for each gene
#for gene in MouseExpression[HumanTissueNames[0]]:
#    GeneExpression[gene] = []
## loop over tissue, add expression for that tissue to the gene expression profile
#for tissue in HumanTissueNames:
#    for gene in MouseExpression[tissue]:
#        assert gene in GeneExpression
#        GeneExpression[gene].append(MouseExpression[tissue][gene])
#
## move back to parent directory
#os.chdir('../') 
#
## check that there are no gene doublons when gene name is parsed to match the GFF gene ID
#genes = [i[:i.index('.')] for i in GeneExpression]
#for i in genes:
#    assert genes.count(i) == 1
#
## make a table with median UP-FPKM in each tissue for each gene
#header = ['Gene'] + HumanTissueNames
#header = '\t'.join(header)
#
#
#
#
#
#
#
#
#
#
#
#
#
#$$$$$XXXXXXXXXXXXXXXX



## create a dict with sample_index: sample_id in normalized FPKM file
#infile = open('Matrix.gtex.rnaseq.fpkmuq.v6p.txt')
#header = infile.readline().rstrip().split('\t')
## do not record the gene ID
#header = header[1:]
#SampleIndex = {}
#for i in range(len(header)):
#    # QC check that all indices have a corresponding tissue
#    assert header[i] in Samples
#    SampleIndex[i] = header[i]
## loop over file, collect expression information for each gene
#Genes = {}
#for line in infile:
#    if line.startswith('ENS'):
#        line = line.rstrip().split('\t')
#        # get gene ID
#        gene = line[0]
#        # remove extra information that is not part of the valid ensembl gene ID
#        if '.' in gene:
#            gene = gene[:gene.index('.')]
#        # check that gene is not already recorded (ie, if information specifies transcript)
#        assert gene not in Genes
#        # record only genes present in GeneCoord
#        if gene in GeneCoord:
#            # record expression information
#            Genes[gene] = line[1:]
#            # check that length of header and expression match
#            assert len(header) == len(Genes[gene])
#infile.close()
#
#
## create a dict {gene: {tissue:[expression values]}}
#Expression = {}
## loop over each gene
#for gene in Genes:
#    # initialize outer dict
#    Expression[gene] = {}
#    # loop over indices of expression list
#    for i in range(len(Genes[gene])):
#        # get the corresponding ID
#        ID = SampleIndex[i]
#        # get the corresponding tissue
#        tissue = Samples[ID]
#        # check if tissue is key in dict
#        if tissue not in Expression[gene]:
#            # initialize list value
#            Expression[gene][tissue] = []
#        # populate list with expression values
#        Expression[gene][tissue].append(float(Genes[gene][i]))
#
## QC that all genes have the same expression profile (ie, same tissue list)
## make a list of expressed genes
#ExpressedGenes = list(Expression.keys())
## get the list o tissues for one of the expressed genes
#TissueExpression = list(Expression[ExpressedGenes[0]].keys())
#TissueExpression.sort()
#for gene in Expression:
#    a = list(Expression[gene].keys())
#    a.sort()
#    assert a == TissueExpression
## make a list of expressed tissues not in tissue types from metada
#MissingInTissueExp = [i for i in Tissues if i not in TissueExpression]
#print('missing in expressed tissues', MissingInTissueExp)
## make a list of metadata tissues not in expressed tissues
#MissingInTissueType = [i for i in TissueExpression if i not in Tissues]
#print('missing in matadata tissues', MissingInTissueType)
#
## an empty string has 10 sample IDs in metadata
## remove empty strings from expression dict
#for gene in Expression:
#    assert '' in Expression[gene]
#    del Expression[gene]['']
## remove empty string from list of expressed tissues
#TissueExpression.remove('')
#TissueExpression.sort()
#
## create a dict with expression profile for each gene
## expression profile is the median expression of given gene in each of the tissues 
## in the same order as the list of tissue types
## create a dict {gene: [median_tissue1, median_tissue2]
#MedianExp = {}
## loop over genes with expression
#for gene in Expression:
#    # initialize list value
#    MedianExp[gene] = []
#    # loop over tissue expression
#    for tissue in TissueExpression:
#        # take the median expression
#        MedianExp[gene].append(np.median(Expression[gene][tissue]))
#print(len(MedianExp))
#
#
## remove genes with no expression in any tissues
#to_remove = []
#for gene in MedianExp:
#    if sum(MedianExp[gene]) == 0:
#        to_remove.append(gene)
#for gene in to_remove:
#    del MedianExp[gene]
#print(len(MedianExp))
#
## create a dict with sample size per tissue
#SampleCounts = {}
## loop over gene with expression
#for gene in Expression:
#    # loop over tissue
#    for tissue in Expression[gene]:
#        if tissue not in SampleCounts:
#            SampleCounts[tissue] = len(Expression[gene][tissue])
#        else:
#            assert SampleCounts[tissue] == len(Expression[gene][tissue])
## print expressed tissues and sample sizes
#for tissue in SampleCounts:
#    print(tissue, SampleCounts[tissue])


## write median normalized fpkm to file
#newfile = open('truc', 'w')
#header = ['Gene'] + TissueExpression
## write header to file
#newfile.write('\t'.join(header) + '\n')
## write gene and expression profile
#for gene in MedianExp:
#    profile = list(map(lambda x: str(x), MedianExp[gene]))
#    line = [gene] + profile
#    newfile.write('\t'.join(line) + '\n')
#newfile.close()
 














#$$$$$$$$$$$$$$$$$$$$$$



            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            