# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 22:27:34 2017

@author: Richard
"""

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



# use this scipt to test for enrichement of disease genes among external genes

# load dictionaries of host and nested genes 
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanNestedGenes.json') as human_json_data:
    Nested = json.load(human_json_data)

# get GFF file
GFF = 'Homo_sapiens.GRCh38.86.gff3'
 
# find nested and intronic-nested genes 
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
# Map Transcript names to gene names {transcript: gene}
MapTranscriptGene = TranscriptToGene(GFF)
# get the coordinates of all exons    
ExonCoord = GeneExonCoord(GFF)
ExonCoord = CleanGeneFeatureCoord(ExonCoord, MapTranscriptGene)
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
IntronCoord = GeneIntronCoord(ExonCoord)
IntronCoord = CleanGeneFeatureCoord(IntronCoord, MapTranscriptGene)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
TranscriptCoordinates = TranscriptsCoord(GFF)
# map genes to their longest transcript {gene: longest_transcript}
GeneLongestTranscript = LongestTranscript(TranscriptCoordinates, MapGeneTranscript)
# match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
Matches = MatchHostTranscriptWithNestedTranscript(Nested, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# list all host, nested transcript pairs [[host, nested]]
HostNestedPairs = GetHostNestedPairs(Matches)



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

# create a set with all disease genes
DiseaseGenes = GAD.union(GWAS).union(Drivers).union(OMIM)


# generate a figure comparing proportions of disease and non-disease genes among
# external genes with internal genes that are intronless or intron-containing

# make sets of external genes with intronless and intron-containing internal genes
WithIntrons, NoIntrons = set(), set()

# loop over host-nested transcript pairs
for i in range(len(HostNestedPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
    if HostNestedPairs[i][1] in IntronCoord:
        # internal gene has introns, populate set with external gene name
            WithIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
    else:
        # internal gene is intronless, populate set with external gene name
        NoIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
    
    

ExtGenes = [WithIntrons, NoIntrons]  
GeneCats = ['With', 'Without'] 



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


# use this function to assign significance level
def AssignSignificance(L):
    '''
    (list) -> list
    Take a list of p-values and return a modfied list with significance levels
    represented by stars
    '''
    # replace P values by significance
    for i in range(len(L)):
        if L[i] >= 0.05:
            L[i] = ''
        elif L[i] < 0.05 and L[i] >= 0.01:
            L[i] = '*'
        elif L[i] < 0.01 and L[i] >= 0.001:
            L[i] = '**'
        elif L[i] < 0.001:
            L[i] = '***'
    return L
    
 
# use this function to get gene proportions
def GetProportions(Counts):
    '''
    (list) -> list, list
    Take the list of inner lists with counts of disease and non-disease
    and return 2 lists with proportions of disease and non-disease genes respectively
    '''
    disease, nondisease = [], []
    for i in range(len(Counts)):
        disease.append(Counts[i][0] / sum(Counts[i]))
        nondisease.append(Counts[i][1] / sum(Counts[i]))
    return disease, nondisease
    
# count disease and non-disease genes    
GADCounts = CountDiseaseGenes(ExtGenes, GAD)    
GWASCounts = CountDiseaseGenes(ExtGenes, GWAS)
DriversCounts = CountDiseaseGenes(ExtGenes, Drivers)
OMIMCounts = CountDiseaseGenes(ExtGenes, OMIM)
AllCounts = CountDiseaseGenes(ExtGenes, DiseaseGenes)

# compare the proportion of disease genes for each set of disease gene
PVals = []
counts = [GADCounts, GWASCounts, DriversCounts, OMIMCounts, AllCounts]
for i in range(len(counts)):
    p = stats.fisher_exact([counts[i][0], counts[i][1]])[1]
    PVals.append(p)    
    print(i, counts[i][0][0]/sum(counts[i][0]), counts[i][1][0] / sum(counts[i][1]), p)
    
# transform p values in stars
PVals = AssignSignificance(PVals)    
    
# create lists of proportions for disease and non-disease genes
DisProp, NonDisProp = [], []
for i in range(len(counts)):
    disease, nondisease = GetProportions(counts[i])
    DisProp.append(disease)
    NonDisProp.append(nondisease)


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, XLabel, YRange, YMax):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot proportions of disease genes only
    # Create a horizontal bar plot for proportions of disease genes
    ax.bar([0, 0.3], Data, width = 0.2, label = 'disease', color= ['black'] + ['lightgrey'] * 6, linewidth = 0.7)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write axis labels    
    YLabel = 'Proportion of\ndisease genes'
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    ax.set_xlabel(XLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # edit y axis ticks
    plt.yticks(YRange)    
    # add a range for the Y axis
    plt.ylim([0, YMax])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right = 'off', left = 'on', labelbottom='off', colors = 'black',
                    labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # add margins
    plt.margins(0.1)
    return ax


# make a figure with proportion of disease and non-disease genes






# create figure
fig = plt.figure(1, figsize = (6, 2.5))
# plot data
ax1 = CreateAx(5, 1, 1, fig, DisProp[0], 'complex diseases', np.arange(0, 1, 0.1), 1)
ax2 = CreateAx(5, 1, 2, fig, DisProp[1], 'GWAS', np.arange(0, 1, 0.1), 1)
ax3 = CreateAx(5, 1, 3, fig, DisProp[2], 'tumor drivers', np.arange(0, 1, 0.1), 1)
ax4 = CreateAx(5, 1, 4, fig, DisProp[3], 'medelian diseases', np.arange(0, 1, 0.1), 1)
ax5 = CreateAx(5, 1, 5, fig, DisProp[4], 'all diseases', np.arange(0, 1, 0.1), 1)

## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
#
#ypos = [0.55, 0.50, 0.65, 0.55, 0.6, 0.6]
#for i in range(len(PValGAD)):
#    ax1.text(xpos[i], ypos[i], PValGAD[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
#ypos = [0.12, 0.05, 0.16, 0.09, 0.10, 0.10]
#for i in range(len(PValGWAS)):
#    ax2.text(xpos[i], ypos[i], PValGWAS[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
#ypos = [0.027, 0.017, 0.040, 0.015, 0.037, 0.030]
#for i in range(len(PValDrivers)):
#    ax3.text(xpos[i], ypos[i], PValDrivers[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
#ypos = [0.17, 0.12, 0.25, 0.17, 0.22, 0.22]
#for i in range(len(PValOMIM)):
#    ax4.text(xpos[i], ypos[i], PValOMIM[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
#ypos = [0.55, 0.45, 0.7, 0.55, 0.65, 0.65]
#for i in range(len(PValAll)):
#    ax5.text(xpos[i], ypos[i], PValAll[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
#
#Proportions = 'disease'
#if Proportions == 'both':
#    # add legend
#    N = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'non-disease')
#    D = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'disease')
#    ax1.legend(handles = [D, N], loc = (0, 1.1), fontsize = 6, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# save figure to file
fig.savefig('truc.pdf', bbox_inches = 'tight')


## make a table with counts of disease and non-disease genes
#GeneCounts = [GADCounts, GWASCounts, DriversCounts, OMIMCounts, AllCounts]
#Origins = ['GAD', 'GWAS', 'Drivers', 'OMIM', 'All']
#GeneCats = ['Non-overlapping', 'Nested', 'Internal', 'External', 'Piggyback', 'Convergent', 'Divergent'] 
#
#newfile = open('truc.txt', 'w')
#header = ['Disease genes', 'Gene category', 'N disease genes', 'N non-disease genes', 'Proportion disease genes', 'P']
#newfile.write('\t'.join(header) + '\n')
#for i in range(len(GeneCounts)):
#    # get the list of P values
#    Pvals = TestDiseaseEnrichement(GeneCounts[i])
#    # add empty  string for the non-overlapping genes
#    Pvals.insert(0, '')
#    for j in range(len(GeneCounts[i])):
#        line = [Origins[i], GeneCats[j], str(GeneCounts[i][j][0]), str(GeneCounts[i][j][1]), str(round(GeneCounts[i][j][0] / sum(GeneCounts[i][j]), 4) * 100), str(Pvals[j])]
#        newfile.write('\t'.join(line) + '\n')
#newfile.close()        





Same, Opposite = set(), set()
# make sets of external genes with same and opposite direction as their internal genes that be intronless or not
WithIntronsSame, NoIntronsSame, WithIntronsOpposite, NoIntronsOpposite = set(), set(), set(), set()




############################## edit below




## loop over host-nested transcript pairs
#for i in range(len(HostNestedPairs)):
#    # check that both transcripts have coordinates and have corresponding gene names    
#    assert HostNestedPairs[i][0] in MapTranscriptGene
#    assert HostNestedPairs[i][0] in TranscriptCoordinates    
#    assert HostNestedPairs[i][1] in MapTranscriptGene
#    assert HostNestedPairs[i][1] in TranscriptCoordinates    
#    if HostNestedPairs[i][1] in IntronCoord:
#            WithIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
#    else:
#        NoIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
#    
#    
#    
#    
#    # get the orientation of each transcript [+,+]
#    OrientationPair = GenePairOrientation(HostNestedPairs[i], TranscriptCoordinates)
#    # populate sets with external genes
#    if len(set(OrientationPair)) == 1:
#        assert set(OrientationPair) == {'-'} or set(OrientationPair) == {'+'}
#        Same.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        # check if the internal gene has introns
#        if HostNestedPairs[i][1] in IntronCoord:
#            WithIntronsSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        else:
#            NoIntronsSame.add(MapTranscriptGene[HostNestedPairs[i][0]])
#    elif len(set(OrientationPair)) == 2:
#        assert set(OrientationPair) == {'-', '+'}
#        Opposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        # check if the internal gene has introns
#        if HostNestedPairs[i][1] in IntronCoord:
#            WithIntronsOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
#        else:
#            NoIntronsOpposite.add(MapTranscriptGene[HostNestedPairs[i][0]])
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
#
#
#       
## load dictionary of overlapping gene pairs
#json_data = open('HumanOverlappingGenes.json')
#Overlapping = json.load(json_data)
#json_data.close()
## make a set of non-overlapping genes
#OverlappingGenes = MakeFullPartialOverlapGeneSet(Overlapping)
#NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlapping, GeneCoord)
#
#
#
#
#
#
#   
##AllGenes = [NonOverlappingGenes, Same, Opposite, WithIntronsSame, WithIntronsOpposite,
##            NoIntronsSame, NoIntronsOpposite]  
##GeneCats = ['NoOvl', 'Same', 'Opposite', 'IntronCis', 'IntronTrans', 'NoIntronCis', 'NoIntronTrans'] 
#
#
#AllGenes = [NonOverlappingGenes, WithIntrons, NoIntrons]  
#GeneCats = ['NoOvl', 'With', 'Without'] 
#
#
#
#
#
## make a list of counts for each 
#
#
## use this function to count disease and non-disease genes for each gene class
#def CountDiseaseGenes(L, DiseaseGenes):
#    '''
#    (list, set) -> list
#    Take a list of gene sets and a set of disease genes and return a parallel
#    list with lists of counts of diease and non disease genes for each gene set
#    '''
#    Counts = []
#    for i in range(len(L)):
#        disease = len([j for j in L[i] if j in DiseaseGenes])
#        nondisease = len([j for j in L[i] if j not in DiseaseGenes])
#        Counts.append([disease, nondisease])
#    return Counts
#
## use this function to test for enrichement of disease genes among gene groups
#def TestDiseaseEnrichement(Counts):
#    '''
#    (list) -> list
#    Take the list of disease and non-disease gene counts for each gene group and 
#    returns a list of p-values from FET comparing each overlapping gene group to 
#    non-overlapping genes
#    Precondition: the non-overlapping gene counts are first in the list    
#    '''
#    PVals = []    
#    for i in range(1, len(Counts)):
#        p = stats.fisher_exact([Counts[0], Counts[i]])[1]
#        PVals.append(p)
#    return PVals
#
## use this function to assign significance level
#def AssignSignificance(L):
#    '''
#    (list) -> list
#    Take a list of p-values and return a modfied list with significance levels
#    represented by stars
#    '''
#    # replace P values by significance
#    for i in range(len(L)):
#        if L[i] >= 0.05:
#            L[i] = ''
#        elif L[i] < 0.05 and L[i] >= 0.01:
#            L[i] = '*'
#        elif L[i] < 0.01 and L[i] >= 0.001:
#            L[i] = '**'
#        elif L[i] < 0.001:
#            L[i] = '***'
#    return L
#    
# 
## use this function to get gene proportions
#def GetProportions(Counts):
#    '''
#    (list) -> list, list
#    Take the list of inner lists with counts of disease and non-disease
#    and return 2 lists with proportions of disease and non-disease genes respectively
#    '''
#    disease, nondisease = [], []
#    for i in range(len(Counts)):
#        disease.append(Counts[i][0] / sum(Counts[i]))
#        nondisease.append(Counts[i][1] / sum(Counts[i]))
#    return disease, nondisease
#    
## count disease and non-disease genes    
#GADCounts = CountDiseaseGenes(AllGenes, GAD)    
#GWASCounts = CountDiseaseGenes(AllGenes, GWAS)
#DriversCounts = CountDiseaseGenes(AllGenes, Drivers)
#OMIMCounts = CountDiseaseGenes(AllGenes, OMIM)
#
## create a set with all disease genes
#DiseaseGenes = set()
#for i in Drivers:
#    DiseaseGenes.add(i)
#for i in GWAS:
#    DiseaseGenes.add(i)
#for i in GAD:
#    DiseaseGenes.add(i)
#for i in OMIM:
#    DiseaseGenes.add(i)
#
#AllCounts = CountDiseaseGenes(AllGenes, DiseaseGenes)
#
#
#############
#a = [GADCounts, GWASCounts, DriversCounts, OMIMCounts, AllCounts]
#for i in range(0, len(a)):
#    for j in range(1, len(a[i])-1):
#        p = stats.fisher_exact([a[i][j], a[i][j+1]])[1]
#        print(i, j, GeneCats[j], GeneCats[j+1], a[i][j][0]/sum(a[i][j]), a[i][j+1][0] / sum(a[i][j+1]), p)
#    
#
#
#
#
#
#
##########
#
#
#
#
#
## test for enrichement of disease genes between non-overlapping genes and overlapping genes
#PValGAD = TestDiseaseEnrichement(GADCounts)
#PValGWAS = TestDiseaseEnrichement(GWASCounts)
#PValDrivers = TestDiseaseEnrichement(DriversCounts)
#PValOMIM = TestDiseaseEnrichement(OMIMCounts)
#PValAll = TestDiseaseEnrichement(AllCounts)
#
#
#for i in a:
#    print(i)
#
#
#
#print(PValGAD)
#print(PValGWAS)
#print(PValDrivers)
#print(PValOMIM)
#print(PValAll)
#
#
#
#
#PValGAD = AssignSignificance(PValGAD)
#PValGWAS = AssignSignificance(PValGWAS)
#PValDrivers = AssignSignificance(PValDrivers)
#PValOMIM = AssignSignificance(PValOMIM)
#PValAll = AssignSignificance(PValAll)
#
#
## get proportions
#GADDis, GADNonDis = GetProportions(GADCounts)    
#GWASDis, GWASNonDis = GetProportions(GWASCounts)
#DriversDis, DriversNonDis = GetProportions(DriversCounts)
#OMIMDis, OMIMNonDis = GetProportions(OMIMCounts)
#AllDis, AllNonDis = GetProportions(AllCounts)
#
#
#
## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, figure, Data, Title, Proportions, YRange, YMax, XLabel):
#    '''
#    Returns a ax instance in figure
#    '''    
#    
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#    # check if plot only disease genes or proportions of disease and non-disease genes
#    if Proportions == 'both':
#        # Create a horizontal bar plot for proportions of disease genes
#        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[0], width = 0.2, label = 'disease', color= 'black', linewidth = 0.7)
#        # Create a horizontal bar plot for proportions of non-disease genes
#        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[1], width = 0.2, bottom = Data[0], label = 'non-disease', color= 'lightgrey', linewidth = 0.7)
#    elif Proportions == 'disease':
#        # plot proportions of disease genes only
#        # Create a horizontal bar plot for proportions of disease genes
#        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[0], width = 0.2, label = 'disease', color= ['black'] + ['lightgrey'] * 6, linewidth = 0.7)
#    elif Proportions == 'non-disease':
#        # plot proportions of non-disease genes only
#        # Create a horizontal bar plot for proportions of non-disease genes
#        ax.bar([0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8], Data[1], width = 0.2, label = 'disease', color= 'lightgrey', linewidth = 0.7)
#
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    # write y axis label
#    if Proportions == 'both':
#        YLabel = 'Proportion'
#    elif Proportions == 'disease':
#        YLabel = 'Proportion of\ndisease genes'
#    elif Proportions == 'non-disease':    
#        YLabel = 'Proportion of\nnon-disease genes'   
#    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
#        
#    # add ticks and lebels
#    if XLabel == True:
#        plt.xticks([0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9], ['NoOvl', 'Nst', 'Int', 'Ext', 'Pgk', 'Con', 'Div'], rotation = 30, size = 7, color = 'black', ha = 'right', **FigFont)
#    elif XLabel == False:
#        plt.xticks([0.1, 0.4, 0.7, 1, 1.3, 1.6, 1.9], [''] * 7, size = 7, color = 'black', ha = 'right', **FigFont)
#    
#    # edit y axis ticks
#    plt.yticks(YRange)    
#        
#    # add title
#    plt.title(Title, color = 'black',  size = 7, ha = 'center', **FigFont)
#    # add a range for the Y axis
#    plt.ylim([0, YMax])    
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)
#    ax.spines["left"].set_visible(True)  
#    # edit tick parameters    
#    if XLabel == True:
#        plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                        right = 'off', left = 'on', labelbottom='on',
#                        colors = 'black', labelsize = 7, direction = 'out')  
#    else:
#        plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                        right = 'off', left = 'on', labelbottom='off',
#                        colors = 'black', labelsize = 7, direction = 'out')  
#    # Set the tick labels font name
#    for label in ax.get_yticklabels():
#        label.set_fontname('Arial')   
#    # add margins
#    plt.margins(0.1)
#    return ax
#
#
### make a figure with proportion of disease and non-disease genes
##
### create figure
##fig = plt.figure(1, figsize = (2.5, 6))
### plot data
##ax1 = CreateAx(1, 5, 1, fig, [GADDis, GADNonDis], 'complex diseases', 'disease', np.arange(0, 0.71, 0.1), 0.71, False)
##ax2 = CreateAx(1, 5, 2, fig, [GWASDis, GWASNonDis], 'GWAS', 'disease', np.arange(0, 0.21, 0.05), 0.2, False)
##ax3 = CreateAx(1, 5, 3, fig, [DriversDis, DriversNonDis], 'tumor drivers', 'disease', np.arange(0, 0.041, 0.010), 0.04,  False)
##ax4 = CreateAx(1, 5, 4, fig, [OMIMDis, OMIMNonDis], 'medelian diseases', 'disease', np.arange(0, 0.26, 0.05), 0.25, False)
##ax5 = CreateAx(1, 5, 5, fig, [AllDis, AllNonDis], 'all diseases', 'disease', np.arange(0, 0.71, 0.1), 0.71, True)
##
### annotate figure to add significance
### significant comparisons were already determined, add letters to show significance
##xpos = [0.4, 0.7, 1, 1.3, 1.6, 1.9]
##
##ypos = [0.55, 0.50, 0.65, 0.55, 0.6, 0.6]
##for i in range(len(PValGAD)):
##    ax1.text(xpos[i], ypos[i], PValGAD[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
##ypos = [0.12, 0.05, 0.16, 0.09, 0.10, 0.10]
##for i in range(len(PValGWAS)):
##    ax2.text(xpos[i], ypos[i], PValGWAS[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
##ypos = [0.027, 0.017, 0.040, 0.015, 0.037, 0.030]
##for i in range(len(PValDrivers)):
##    ax3.text(xpos[i], ypos[i], PValDrivers[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
##ypos = [0.17, 0.12, 0.25, 0.17, 0.22, 0.22]
##for i in range(len(PValOMIM)):
##    ax4.text(xpos[i], ypos[i], PValOMIM[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
##ypos = [0.55, 0.45, 0.7, 0.55, 0.65, 0.65]
##for i in range(len(PValAll)):
##    ax5.text(xpos[i], ypos[i], PValAll[i], ha='center', va='center', color = 'grey', fontname = 'Arial', size = 7)
##
##Proportions = 'disease'
##if Proportions == 'both':
##    # add legend
##    N = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'non-disease')
##    D = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'disease')
##    ax1.legend(handles = [D, N], loc = (0, 1.1), fontsize = 6, frameon = False, ncol = 2)
##
### make sure subplots do not overlap
##plt.tight_layout()
##
### save figure to file
##fig.savefig('truc.pdf', bbox_inches = 'tight')
#
#
## make a table with counts of disease and non-disease genes
#
#
##GeneCounts = [GADCounts, GWASCounts, DriversCounts, OMIMCounts, AllCounts]
##Origins = ['GAD', 'GWAS', 'Drivers', 'OMIM', 'All']
##GeneCats = ['Non-overlapping', 'Nested', 'Internal', 'External', 'Piggyback', 'Convergent', 'Divergent'] 
##
##newfile = open('truc.txt', 'w')
##header = ['Disease genes', 'Gene category', 'N disease genes', 'N non-disease genes', 'Proportion disease genes', 'P']
##newfile.write('\t'.join(header) + '\n')
##for i in range(len(GeneCounts)):
##    # get the list of P values
##    Pvals = TestDiseaseEnrichement(GeneCounts[i])
##    # add empty  string for the non-overlapping genes
##    Pvals.insert(0, '')
##    for j in range(len(GeneCounts[i])):
##        line = [Origins[i], GeneCats[j], str(GeneCounts[i][j][0]), str(GeneCounts[i][j][1]), str(round(GeneCounts[i][j][0] / sum(GeneCounts[i][j]), 4) * 100), str(Pvals[j])]
##        newfile.write('\t'.join(line) + '\n')
##newfile.close()        
#
#
#
#
#
#Same, Opposite = set(), set()
## make sets of external genes with same and opposite direction as their internal genes that be intronless or not
#WithIntronsSame, NoIntronsSame, WithIntronsOpposite, NoIntronsOpposite = set(), set(), set(), set()
