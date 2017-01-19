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

# create lists of nested gene pairs with same and opposite directions
same, opposite = [], []
for pair in NestedPairs:
    orientation = GenePairOrientation(pair, GeneCoord)
    if len(set(orientation)) == 2:
        opposite.append(pair)
    elif len(set(orientation)) == 1:
        same.append(pair)
# create sets of internal and external nested genes depending on orientation 
InternalSameGenes, InternalOppositeGenes, ExternalSameGenes, ExternalOppositeGenes = set(), set(), set(), set()
for pair in same:
    ExternalSameGenes.add(pair[0])
    InternalSameGenes.add(pair[1])
for pair in opposite:
    ExternalOppositeGenes.add(pair[0])
    InternalOppositeGenes.add(pair[1])

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


# make a set of cancer driver genes
Drivers = set()
infile = open('driver_genes_per_tumor_syn7314119.csv')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        Drivers.add(line[2][:line[2].index('.')])
infile.close()


# make a set of mendelian disease genes
OMIM = set()
discat = set()
# get the mim IDs corresponding to associations between phenotypes and genes
mimIDs = set()
infile = open('mimTitles.txt')
for line in infile:
    if (not line.startswith('#')) and line.rstrip() != '':
        line = line.rstrip().split('\t')
        if line[0] == 'Number Sign' or line[0] == 'Percent' or line[0] == 'Plus':
            mimIDs.add(line[1])
            
infile.close()

print(len(mimIDs))

truc = set()

# get the set of phenotype associated genes
infile = open('morbidmap.txt')
for line in infile:
    if (not line.startswith('#')) and line.rstrip() != '':
        line = line.rstrip().split('\t')
        pheno = line[0].replace(' ', '')
        pheno = pheno.split(',')
        pheno = pheno[-1]
        pheno = pheno[:pheno.index('(')]
        if pheno in mimIDs:
            for item in line[1:-2]:
                if item in GeneNames:
                    OMIM.add(GeneNames[item])
                    truc.add(item)
infile.close()

print('omim', len(OMIM), len(truc))          
print(truc)


AllGenes = [NonOverlappingGenes, NestedGenes, InternalGenes, ExternalGenes,
            PiggyBackGenes, ConvergentGenes, DivergentGenes] 

GeneCats = ['NoOvl', 'Nst', 'Int', 'Ext', 'Pgk', 'Con', 'Div'] 



print('GAD genes')
for i in range(1, len(AllGenes)):
    # compute the number of disease and non-disease genes
    DiseaseNonOv = len([j for j in AllGenes[0] if j in GADID])
    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in GADID])
    disease = len([j for j in AllGenes[i] if j in GADID])
    nondisease = len([j for j in AllGenes[i] if j not in GADID])
    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
          disease, nondisease, round(disease / (disease + nondisease), 4), p)
     
print('GWAS genes')
for i in range(1, len(AllGenes)):
    # compute the number of disease and non-disease genes
    DiseaseNonOv = len([j for j in AllGenes[0] if j in GWASID])
    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in GWASID])
    disease = len([j for j in AllGenes[i] if j in GWASID])
    nondisease = len([j for j in AllGenes[i] if j not in GWASID])
    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
    #print(i, GeneCats[i], round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), p)
    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
          disease, nondisease, round(disease / (disease + nondisease), 4), p)


print('driver genes')
for i in range(1, len(AllGenes)):
    # compute the number of disease and non-disease genes
    DiseaseNonOv = len([j for j in AllGenes[0] if j in Drivers])
    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in Drivers])
    disease = len([j for j in AllGenes[i] if j in Drivers])
    nondisease = len([j for j in AllGenes[i] if j not in Drivers])
    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
    #print(i, GeneCats[i], round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), p)
    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
          disease, nondisease, round(disease / (disease + nondisease), 4), p)








# create a set with all disease genes
DiseaseGenes = set()
for i in Drivers:
    DiseaseGenes.add(i)
for i in GWASID:
    DiseaseGenes.add(i)
for i in GADID:
    DiseaseGenes.add(i)
for i in OMIM:
    DiseaseGenes.add(i)



print('all genes')
for i in range(1, len(AllGenes)):
    # compute the number of disease and non-disease genes
    DiseaseNonOv = len([j for j in AllGenes[0] if j in DiseaseGenes])
    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in DiseaseGenes])
    disease = len([j for j in AllGenes[i] if j in DiseaseGenes])
    nondisease = len([j for j in AllGenes[i] if j not in DiseaseGenes])
    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
    #print(i, GeneCats[i], round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), p)
    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
          disease, nondisease, round(disease / (disease + nondisease), 4), p)
    