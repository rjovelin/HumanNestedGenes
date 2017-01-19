# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 14:54:24 2017

@author: RJovelin
"""

# use this script to test for enrichement of overlapping genes among disease genes
# save data as table and figure

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
            if gene in GeneNames:
                GAD.add(GeneNames[gene])
        elif type(gene) == list:
            for item in gene:
                if item in GeneNames:
                    GAD.add(item)
        # get alternative gene names
        alternative = line[1].split('|')
        while '"' in alternative:
            alternative.remove('"')
        while '' in alternative:
            alternative.remove('')
        for name in alternative:
            if name in GeneNames:
                GAD.add(GeneNames[name])
infile.close()
print(len(GAD))


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
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # check that trait is disease only        
        if line[7] not in ExcludeTraits:
            genes = line[13].replace(' ', '')
            # check if multiple genes are listed
            if ', ' in genes:
                genes = genes.split(',')
            if type(genes) == str:
                if genes in GeneNames:
                    GWAS.add(GeneNames[genes])
            elif type(genes) == list:
                for item in genes:
                    if item in GeneNames:
                        GWAS.add(GeneNames[item])       
infile.close()
print(len(GWAS))


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
            # extract the genes
            genes = line[-3].replace(' ', '')
            # check if multiple genes are listed
            if ',' in genes:
                genes = genes.split(',')
            if type(genes) == str:
                if genes in GeneNames:
                    OMIM.add(GeneNames[genes])
            elif type(genes) == list:
                for item in genes:
                    if item in GeneNames:
                        OMIM.add(GeneNames[item])
infile.close()

   


AllGenes = [NonOverlappingGenes, NestedGenes, InternalGenes, ExternalGenes,
            PiggyBackGenes, ConvergentGenes, DivergentGenes] 

GeneCats = ['NoOvl', 'Nst', 'Int', 'Ext', 'Pgk', 'Con', 'Div'] 







#print('GAD genes')
#for i in range(1, len(AllGenes)):
#    # compute the number of disease and non-disease genes
#    DiseaseNonOv = len([j for j in AllGenes[0] if j in GADID])
#    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in GADID])
#    disease = len([j for j in AllGenes[i] if j in GADID])
#    nondisease = len([j for j in AllGenes[i] if j not in GADID])
#    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
#    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
#          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
#          disease, nondisease, round(disease / (disease + nondisease), 4), p)
#     
#print('GWAS genes')
#for i in range(1, len(AllGenes)):
#    # compute the number of disease and non-disease genes
#    DiseaseNonOv = len([j for j in AllGenes[0] if j in GWASID])
#    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in GWASID])
#    disease = len([j for j in AllGenes[i] if j in GWASID])
#    nondisease = len([j for j in AllGenes[i] if j not in GWASID])
#    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
#    #print(i, GeneCats[i], round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), p)
#    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
#          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
#          disease, nondisease, round(disease / (disease + nondisease), 4), p)
#
#
#print('driver genes')
#for i in range(1, len(AllGenes)):
#    # compute the number of disease and non-disease genes
#    DiseaseNonOv = len([j for j in AllGenes[0] if j in Drivers])
#    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in Drivers])
#    disease = len([j for j in AllGenes[i] if j in Drivers])
#    nondisease = len([j for j in AllGenes[i] if j not in Drivers])
#    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
#    #print(i, GeneCats[i], round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), p)
#    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
#          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
#          disease, nondisease, round(disease / (disease + nondisease), 4), p)
#
#
#print('OMIM genes')
#for i in range(1, len(AllGenes)):
#    # compute the number of disease and non-disease genes
#    DiseaseNonOv = len([j for j in AllGenes[0] if j in OMIM])
#    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in OMIM])
#    disease = len([j for j in AllGenes[i] if j in OMIM])
#    nondisease = len([j for j in AllGenes[i] if j not in OMIM])
#    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
#    #print(i, GeneCats[i], round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), p)
#    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
#          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
#          disease, nondisease, round(disease / (disease + nondisease), 4), p)
#
## create a set with all disease genes
#DiseaseGenes = set()
#for i in Drivers:
#    DiseaseGenes.add(i)
#for i in GWASID:
#    DiseaseGenes.add(i)
#for i in GADID:
#    DiseaseGenes.add(i)
#for i in OMIM:
#    DiseaseGenes.add(i)
#
#
#
#print('all genes')
#for i in range(1, len(AllGenes)):
#    # compute the number of disease and non-disease genes
#    DiseaseNonOv = len([j for j in AllGenes[0] if j in DiseaseGenes])
#    NonDiseaseNonOV = len([j for j in AllGenes[0] if j not in DiseaseGenes])
#    disease = len([j for j in AllGenes[i] if j in DiseaseGenes])
#    nondisease = len([j for j in AllGenes[i] if j not in DiseaseGenes])
#    p = stats.fisher_exact([[NonDiseaseNonOV, DiseaseNonOv], [nondisease, disease]])[1]
#    #print(i, GeneCats[i], round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3), round(disease / (disease + nondisease), 4), p)
#    print(i, GeneCats[i], DiseaseNonOv, NonDiseaseNonOV, 
#          round(DiseaseNonOv / (DiseaseNonOv + NonDiseaseNonOV), 3),
#          disease, nondisease, round(disease / (disease + nondisease), 4), p)


# make a list of counts for each 


# use this function to count disease and non-disease genes for each gene class
def CountDiseaseGenes(L, DiseaseGenes):
    '''
    (list, set) -> list
    Take a list of gene sets and a set of disease genes and return a parallel
    list with lists of counts of diease and non disease genes for each gene set
    '''
    Counts = []
    for i in range(len(L)):
        disease = len([j for j in L[i] if j in DiseaseGenes])
        nondisease = len([j for j in L[i] if j not in DiseaseGenes])
        Counts.append([disease, nondisease])
    return Counts

# use this function to test for enrichement of disease genes among gene groups
def TestDiseaseEnrichement(Counts):
    '''
    (list) -> list
    Take the list of disease and non-disease gene counts for each gene group and 
    returns a list of p-values from FET comparing each overlapping gene group to 
    non-overlapping genes
    Precondition: the non-overlapping gene counts are first in the list    
    '''
    PVals = []    
    for i in range(1, len(Counts)):
        p = stats.fisher_exact([Counts[0], Counts[i]])[1]
        PVals.append(p)
    return PVals
    
# count disease and non-disease genes    
GADCounts = CountDiseaseGenes(AllGenes, GAD)    
GWASCounts = CountDiseaseGenes(AllGenes, GWAS)
DriversCounts = CountDiseaseGenes(AllGenes, Drivers)
OMIMCounts = CountDiseaseGenes(AllGenes, OMIM)

# create a set with all disease genes
DiseaseGenes = set()
for i in Drivers:
    DiseaseGenes.add(i)
for i in GWAS:
    DiseaseGenes.add(i)
for i in GAD:
    DiseaseGenes.add(i)
for i in OMIM:
    DiseaseGenes.add(i)

AllCounts = CountDiseaseGenes(AllGenes, DiseaseGenes)


# test for enrichement of disease genes between non-overlapping genes and overlapping genes
PValGAD = TestDiseaseEnrichement(GADCounts)
PValGWAS = TestDiseaseEnrichement(GWASCounts)
PValDrivers = TestDiseaseEnrichement(DriversCounts)
PvalOMIM = TestDiseaseEnrichement(OMIMCounts)
PValAll = TestDiseaseEnrichement(AllCounts)





# test for enrichement of disease genes between non-overlapping genes and overlapping genes
PValGAD = TestDiseaseEnrichement(GADCounts)
PValGWAS = TestDiseaseEnrichement(GWASCounts)
PValDrivers = TestDiseaseEnrichement(DriversCounts)
PValOMIM = TestDiseaseEnrichement(OMIMCounts)
PValAll = TestDiseaseEnrichement(AllCounts)


for i in range(1, len(GADCounts)):
    print('GAD', GeneCats[0], GeneCats[i], GADCounts[0][0]/sum(GADCounts[0]), GADCounts[i][0] / sum(GADCounts[i]), PValGAD[i-1])
    print('GWAS', GeneCats[0], GeneCats[i], GWASCounts[0][0] / sum(GWASCounts[0]), GWASCounts[i][0] / sum(GWASCounts[i]), PValGWAS[i-1])
    print('divers', GeneCats[0], GeneCats[i], DriversCounts[0][0] / sum(DriversCounts[0]), DriversCounts[i][0] / sum(DriversCounts[i]), PValDrivers[i-1])
    print('OMIM', GeneCats[0], GeneCats[i], OMIMCounts[0][0] / sum(OMIMCounts[0]), OMIMCounts[i][0] / sum(OMIMCounts[i]), PValOMIM[i-1])
    print('all', GeneCats[0], GeneCats[i], AllCounts[0][0] / sum(AllCounts[0]), AllCounts[i][0] / sum(AllCounts[i]), PValAll[i-1])





#########################

            

## count genes with and without homologs
#GeneCounts = []
## loop over gene sets
#for i in range(len(AllGenes)):
#    # initialize counters
#    homo, nohomo = 0, 0
#    # loop over genes in given set
#    for gene in AllGenes[i]:
#        if gene in Homologs:
#            homo += 1
#        else:
#            nohomo += 1
#    GeneCounts.append([homo, nohomo])    
#
## compare the proportions of gene with and without homologs
## create a list to store the P-values
#PProp = []
#for i in range(1, len(GeneCounts)):
#    p = stats.fisher_exact([GeneCounts[0], GeneCounts[i]])[1]
#    PProp.append(p)
## replace P values by significance
#for i in range(len(PProp)):
#    if PProp[i] >= 0.05:
#        PProp[i] = ''
#    elif PProp[i] < 0.05 and PProp[i] >= 0.01:
#        PProp[i] = '*'
#    elif PProp[i] < 0.01 and PProp[i] >= 0.001:
#        PProp[i] = '**'
#    elif PProp[i] < 0.001:
#        PProp[i] = '***'
#
## get the proportions of genes with and without homologs
#WithHomolog, NoHomolog = [], []
#for i in range(len(GeneCounts)):
#    WithHomolog.append(GeneCounts[i][0] / sum(GeneCounts[i]))
#    NoHomolog.append(GeneCounts[i][1] / sum(GeneCounts[i]))
#    assert sum(GeneCounts[i]) == len(AllGenes[i])
#
#
## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YLabel, DataType, YMax):
#    '''
#    Returns a ax instance in figure
#    '''    
#
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#    # check type of graphic    
#    if DataType == 'divergence':
#        # set colors
#        colorscheme = ['black','lightgrey','lightgrey','lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
#        # plot nucleotide divergence
#        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[0], 0.2, yerr = Data[1], color = colorscheme,
#               edgecolor = 'black', linewidth = 0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))
#    elif DataType == 'proportion':
#        ## Create a horizontal bar plot for proportions of genes with homologs
#        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[0], width = 0.2, label = 'homolog', color= 'black', linewidth = 0.7)
#        # Create a horizontal bar plot for proportions of same strand pairs
#        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[1], width = 0.2, bottom = Data[0], label = 'no homolog', color= 'lightgrey', linewidth = 0.7)
#
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write y axis label
#    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
#    # add ticks and lebels
#    plt.xticks([0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9], XLabel, rotation = 30, size = 7, color = 'black', ha = 'right', **FigFont)
#    # add a range for the Y and X axes
#    plt.ylim([0, YMax])    
#    
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)
#    ax.spines["left"].set_visible(True)  
#    # edit tick parameters    
#    plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                    right = 'off', left = 'on', labelbottom='on',
#                    colors = 'black', labelsize = 7, direction = 'out')  
#    # Set the tick labels font name
#    for label in ax.get_yticklabels():
#        label.set_fontname('Arial')   
#      
#    # add margins
#    plt.margins(0.1)
#    
#    return ax
#
#
#
## make a figure with mean dN/dS and with proportion of gene with homologs
#
## create figure
#fig = plt.figure(1, figsize = (4.5, 2))
## plot data
#ax1 = CreateAx(2, 1, 1, fig, [MeanOmega, SEMOmega], GeneCats, 'Nucleotide divergence (dN/dS)', 'divergence', 0.50)
#ax2 = CreateAx(2, 1, 2, fig, [WithHomolog, NoHomolog], GeneCats, 'Proportion', 'proportion', 1)
#
## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#ypos = [0.47, 0.50, 0.47, 0.47, 0.45, 0.45]
#xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
#for i in range(len(PValOmega)):
#    ax1.text(xpos[i], ypos[i], PValOmega[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
#for i in range(len(PProp)):
#    ax2.text(xpos[i], 1.02, PValOmega[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
#
## add legend
#NoH = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'no homolog')
#WiH = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'homolog')
#ax2.legend(handles = [WiH, NoH], loc = (0, 1.1), fontsize = 6, frameon = False, ncol = 2)
#
## make sure subplots do not overlap
#plt.tight_layout()
#
#
#
