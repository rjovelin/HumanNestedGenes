# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 21:26:13 2017

@author: Richard
"""



# use this script to plot expression divergence between host and nested genes 
# separately for intron-containing and intronless nested genes 


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


# 1) compare expression divergence between host and nested genes for intronless and intron-containing internal genes

# map genes to their longest transcript {gene: longest_transcript}
GeneLongestTranscript = LongestTranscript(TranscriptCoordinates, MapGeneTranscript)
# match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
Matches = MatchHostTranscriptWithNestedTranscript(Nested, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# list all host, nested transcript pairs [[host, nested]]
HostNestedPairs = GetHostNestedPairs(Matches)

# make pairs of host-nested genes with intronless and intron-containing internal genes
PairsWithIntrons, PairsNoIntrons = [], []
# loop over host-nested transcript pairs
for i in range(len(HostNestedPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
    if HostNestedPairs[i][1] in IntronCoord:
        # internal gene has introns, add gene pairs to list 
            PairsWithIntrons.append([MapTranscriptGene[HostNestedPairs[i][0]], MapTranscriptGene[HostNestedPairs[i][1]]])
    else:
        # internal gene is intronless, add gene pair to list
        PairsNoIntrons.append([MapTranscriptGene[HostNestedPairs[i][0]], MapTranscriptGene[HostNestedPairs[i][1]]])

# parse the GTEX expression summary file to obtain the expression profile of each gene
ExpressionProfile = ParseGTEXExpression('GTEX_Median_Normalized_FPKM.txt')
# remove genes without any expression
ExpressionProfile = RemoveGenesLackingExpression(ExpressionProfile)
# transform absulte expression in relative expression
ExpressionProfile = TransformRelativeExpression(ExpressionProfile)

# filter gene pairs lacking expression
AllPairs = [PairsWithIntrons, PairsNoIntrons]
for i in range(len(AllPairs)):
    AllPairs[i] = FilterGenePairsWithoutExpression(AllPairs[i], ExpressionProfile)

# compute expression divergence between pairs of genes
ExpressDivergence = []
for i in range(len(AllPairs)):
    Div = ComputeExpressionDivergenceGenePairs(AllPairs[i], ExpressionProfile)
    ExpressDivergence.append(Div)

# make a list of gene category names parallel to the list of gene pairs
PairsCats = ['WithIntrons', 'Intronless']
# create lists with means and SEM for each gene category
MeanExpDiv, SEMExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(ExpressDivergence)):
    MeanExpDiv.append(np.mean(ExpressDivergence[i]))
    SEMExpDiv.append(np.std(ExpressDivergence[i]) / math.sqrt(len(ExpressDivergence[i])))

print(MeanExpDiv)

# perform statistical tests between gene categories
PExpressDiv= PermutationResampling(Divergence[0], Divergence[1], 10000, np.mean)
print(PExpressDiv)




# 2) compare sequence divergence between internal genes with and witout introns
# and between external genes for which internal genes have introns or not

# make sets of host-nested genes with intronless and intron-containing internal genes
ExtWithIntrons, ExtNoIntrons = set(), set()
IntWithIntrons, IntNoIntrons = set(), set()
# loop over host-nested transcript pairs
for i in range(len(HostNestedPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
    if HostNestedPairs[i][1] in IntronCoord:
        # internal gene has introns, add genes to sets 
        ExtWithIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
        IntWithIntrons.add(MapTranscriptGene[HostNestedPairs[i][1]])
    else:
        # internal gene is intronless, add genes to sets
        ExtNoIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
        IntNoIntrons.add(MapTranscriptGene[HostNestedPairs[i][1]])

# get 1:1 orthologs between human and chimp
Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')

# create lists of orthologous pairs for each gene category 
GeneCats = ['ExtWith', 'ExtNo', 'IntWith', 'IntNo']
AllGenes = [ExtWithIntrons, ExtNoIntrons , IntWithIntrons, IntNoIntrons]
# create pairs of orthologs
Orthologs = []
# loop over gene sets
for i in range(len(AllGenes)):
    # create a list of gene pairs
    orthologs = []
    # loop over genes in given gene set
    for gene in AllGenes[i]:
        # check that gene has ortholog
        if gene in Orthos:
            orthologs.append([gene, Orthos[gene]])
    Orthologs.append(orthologs)

# create a dict with divergence values
SeqDiv = {}
infile = open('HumanChimpSeqDiverg.txt')
infile.readline()
for line in infile:
    if line.startswith('ENSG'):
        line = line.rstrip().split('\t')
        if line[-1] != 'NA':
            SeqDiv[line[0]] = float(line[-1])
infile.close()            
    
# make list of dN/dS for each gene category
Omega = []
for i in range(len(Orthologs)):
    ratio = [] 
    for pair in Orthologs[i]:
        if pair[0] in SeqDiv:
            ratio.append(SeqDiv[pair[0]])
    Omega.append(ratio)
    
 # create lists with means and SEM for dN/dS for each gene category
MeanOmega, SEMOmega = [], []
for i in range(len(Omega)):
    MeanOmega.append(np.mean(Omega[i]))
    SEMOmega.append(np.std(Omega[i]) / math.sqrt(len(Omega[i])))

print(GeneCats)
print(MeanOmega)

# test differences between external genes and between internal genes
POmega = []
for i in range(0, len(Omega), 2):
    P = PermutationResampling(Omega[i], Omega[i+1], 10000, np.mean)
    POmega.append(P)
for i in POmega:
    print(i)


###############################

# 3) Compare expression divergence between orthologs for internal genes
# and for external genes 

# create sets of internal and external nested gene pairs
NestedPairs = GetHostNestedPairs(Nested)
InternalGenes, ExternalGenes = set(), set()
for pair in NestedPairs:
    ExternalGenes.add(pair[0])
    InternalGenes.add(pair[1])
    
# get expression profile of human genes
HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
# remove genes wuthout expression
HumanExpression = RemoveGenesLackingExpression(HumanExpression)
# get relative expression
HumanExpression = TransformRelativeExpression(HumanExpression)
# get expression profile of chimp genes
ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
# remove genes wuthout expression
ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
# get relative expression
ChimpExpression = TransformRelativeExpression(ChimpExpression)

# generate lists of ortholog pairs for each gene category
ExtWithIntronsOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, ExtWithIntrons, Orthos)
ExtNoIntronsOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, ExtNoIntrons, Orthos)
IntWithIntronsOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, IntWithIntrons, Orthos)
IntNoIntronsOrthos = ExpressedOrthologousPairs(HumanExpression, ChimpExpression, IntNoIntrons, Orthos) 

# make a list of orthologous expressed gene pairs
OrthoPairs = [ExtWithIntronsOrthos, ExtNoIntronsOrthos, IntWithIntronsOrthos, IntNoIntronsOrthos] 

# compute expression divergence between orthologs for each gene category
OrthoExprxDiverg = []
for i in range(len(OrthoPairs)):
    Div = ComputeExpressionDivergenceOrthologs(OrthoPairs[i], HumanExpression, ChimpExpression)
    OrthoExprxDiverg.append(Div)
    
# create lists with means and SEM for each gene category
MeanOrthoExpDiv, SEMOrthoExpDiv = [], []
# loop over lists in Divergence list
for i in range(len(OrthoExprxDiverg)):
    MeanOrthoExpDiv.append(np.mean(OrthoExprxDiverg[i]))
    SEMOrthoExpDiv.append(np.std(OrthoExprxDiverg[i]) / math.sqrt(len(OrthoExprxDiverg[i])))

# test differences between external genes and between internal genes
POrthosExprx = []
for i in range(0, len(OrthoExprxDiverg), 2):
    P = PermutationResampling(OrthoExprxDiverg[i], OrthoExprxDiverg[i+1], 10000, np.mean)
    POrthosExprx.append(P)
for i in POrthosExprx:
    print(i)

    
  




###############################

## create figure
#fig = plt.figure(1, figsize = (3, 2))
#
## add a plot to figure (N row, N column, plot N)
#ax = fig.add_subplot(1, 1, 1)
## set colors
#colorscheme = ['grey','grey','grey','grey', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
## plot nucleotide divergence
#ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], MeanExpDiv, 0.2, yerr = SEMExpDiv, color = colorscheme,
#       edgecolor = 'black', linewidth = 0.5,
#       error_kw=dict(elinewidth=0.5, ecolor='black', markeredgewidth = 0.5))
## set font for all text in figure
#FigFont = {'fontname':'Arial'}   
## write y axis label
#ax.set_ylabel('Expression divergence', color = 'black',  size = 7, ha = 'center', **FigFont)
## add ticks and lebels
#plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
## add a range for the Y and X axes
#plt.ylim([0, 0.61])
#plt.xlim([0, 2.45])
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
#
#
## annotate figure to add significance
## significant comparisons were already determined, add letters to show significance
#Diff = ['A', 'B', 'C', 'B', 'D', 'E', 'F', 'G']
#ypos = [0.55] * 8
#xpos = [0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25]
#for i in range(len(Diff)):
#    ax.text(xpos[i], ypos[i], Diff[i], ha='center', va='center', color = 'black', fontname = 'Arial', size = 7)
#    
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')



###########################


    
    

## perform statistical tests between gene categories using Kolmogorov-Smirnof test
## create list to store the p-values
#PValOmega = []
#for i in range(1, len(Omega)):
#    # compare each gene category to non-overlapping genes
#    val, P = stats.ks_2samp(Omega[0], Omega[i])
#    PValOmega.append(P)
## replace P values with significance
#for i in range(len(PValOmega)):
#    if PValOmega[i] >= 0.05:
#        PValOmega[i] = ''
#    elif PValOmega[i] < 0.05 and PValOmega[i] >= 0.01:
#        PValOmega[i] = '*'
#    elif PValOmega[i] < 0.01 and PValOmega[i] >= 0.001:
#        PValOmega[i] = '**'
#    elif PValOmega[i] < 0.001:
#        PValOmega[i] = '***'
#
## create a set of human genes that have homologs
## include 1:1 and 1 to many or many to many orthologs
#Homologs = set()
#infile = open('HumanChimpOrthologs.txt')
#infile.readline()
#for line in infile:
#    if 'ortholog' in line:
#        line = line.rstrip().split('\t')
#        assert 'ENS' in line[0], 'gene id is not valid'
#        assert 'ortholog' in line[4], 'ortholog should be in homology type'
#        Homologs.add(line[0])
#infile.close()            
#
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
