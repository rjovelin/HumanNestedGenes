# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:54:24 2017

@author: RJovelin
"""

# use this script to test for enrichement of overlapping genes among disease genes


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


# map ensembl gene IDs to gene names
GeneNames = {}
infile = open('Homo_sapiens.GRCh38.86.gff3')
# consider only protein coding genes
for line in infile:
    if 'gene' in line and not line.startswith('#'):
        line = line.rstrip().split('\t')
        if line[2] == 'gene':
            biotype = line[8][line[8].index('biotype=')+8: line[8].index(';', line[8].index('biotype=')+8)]
            if biotype == 'protein_coding':
                # get gene ID
                gene = line[8][line[8].index('ID=gene:')+8: line[8].index(';')]
                # get gene name
                name = line[8][line[8].index('Name=')+len('Name='): line[8].index(';biotype')]  
                assert gene not in GeneNames
                assert name not in GeneNames.values()
                GeneNames[name] = gene
infile.close()

print(len(GeneNames))


# make a set of complex disease genes
GAD = set()
infile = open('GADCDC_data.tsv')
header = infile.readline().rstrip().split('\t')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # get gene name
        gene = line[5]
        if '"' in gene:
            assert gene.count('"') == 2 and gene[0] == '"' and gene[-1] == '"'
            gene = gene[1:-1]
        # remove space in gene
        gene = gene.replace(' ', '')
        # check if multiple genes are listed
        if ',' in gene:
            assert ':' not in gene and ';' not in gene
            gene = gene.split(',')
        elif ';' in gene:
            assert ':' not in gene and ',' not in gene
            gene = gene.split(';')
        elif ':' in gene:
            assert ';' not in gene and ',' not in gene
            gene = gene.split(':')
        if type(gene) == str:
            GAD.add(gene)
        elif type(gene) == list:
            for item in gene:
                GAD.add(item)
#        # get alternative gene names
#        alternative = line[1].split('|')
#        while '"' in alternative:
#            alternative.remove('"')
#        while '' in alternative:
#            alternative.remove('')
#        for name in alternative:
#            GAD.add(name)
infile.close()
print(len(GAD))

GADID = set()
for name in GAD:
    if name in GeneNames:
        GADID.add(GeneNames[name])
print('GAD', len(GADID))


# make a set of GWAS disease genes
GWAS = set()
# exclude non-disease traits
ExcludeTraits = set()
infile = open('TraitsToRemove.txt')
for line in infile:
    if line.rstrip() != '':
        line = line.strip()
        ExcludeTraits.add(line)
infile.close()
# open GWAS catalog
infile = open('gwas_catalog_v1.0-associations_e87_r2017-01-09.tsv', encoding='utf8')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # check that trait is disease only        
        if line[7] not in ExcludeTraits:
            GWAS.add(line[13])
infile.close()
print(len(GWAS))


GWASID = set()
for name in GWAS:
    if name in GeneNames:
        GWASID.add(GeneNames[name])
print('GWAS', len(GWASID))





AllGenes = [NonOverlappingGenes, NestedGenes, InternalGenes, ExternalGenes,
            PiggyBackGenes, ConvergentGenes, DivergentGenes] 

print('GAD genes')
for i in range(1, len(AllGenes)):
    # compute the number of disease and non-disease genes
    DiseaseNonOv = len(AllGenes[0].intersection(GADID))
    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in GADID])
    disease = len([j for j in AllGenes[i] if j in GADID])
    nondisease = len([j for j in AllGenes[i] if j not in GADID])
    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
    print(i, round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), P)
     


print('GWAS genes')
for i in range(1, len(AllGenes)):
    # compute the number of disease and non-disease genes
    DiseaseNonOv = len(AllGenes[0].intersection(GWASID))
    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in GWASID])
    disease = len([j for j in AllGenes[i] if j in GWASID])
    nondisease = len([j for j in AllGenes[i] if j not in GWASID])
    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
    print(i, round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), P)







    