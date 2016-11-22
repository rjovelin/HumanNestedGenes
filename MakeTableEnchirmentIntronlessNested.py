# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 11:57:29 2016

@author: RJovelin
"""


# use this script to test for over-representation of intronless nested genes in same sense orientation
# generate a contagency table with number of nested intronless-same sense, intronless opposite sense,
# with introns same sense, and with introns opposite sense

# usage python3 MakeTableEnrichmentIntronlessNested.py


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


# load dictionaries of host and nested genes 
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('ChimpHostNestedGenes.json') as chimp_json_data:
    ChimpHostGenes = json.load(chimp_json_data)
with open('GorillaHostNestedGenes.json') as gorilla_json_data:
    GorillaHostGenes = json.load(gorilla_json_data)
with open('OrangOutanHostNestedGenes.json') as orangoutan_json_data:
    OrangOutanHostGenes = json.load(orangoutan_json_data)
with open('MacaqueHostNestedGenes.json') as macaque_json_data:
    MacaqueHostGenes = json.load(macaque_json_data)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
PabGFF = 'Pongo_abelii.PPYG2.86.gff3'
MmlGFF = 'Macaca_mulatta.Mmul_8.0.1.86.gff3' 


# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF, PabGFF, MmlGFF]
# make a list of species names
SpeciesNames = ['human', 'chimp', 'gorilla', 'orangoutan', 'macaque']
# make a list of host:nested genes dictionaries
HostGenes = [HumanHostGenes, ChimpHostGenes, GorillaHostGenes, OrangOutanHostGenes, MacaqueHostGenes]


# create lists to store the number of introns for the longest transcript of the host,
# nested and un-nested genes in each species [[human], [chimp], [gorilla], [orang-outan], [macaque]]
HostIntron, NestedIntron, OthersIntron = [], [], []


# loop over GFF files, find nested and intronic=nested genes in each species 
for i in range(len(GFFs)):
    print(GFFs[i][:GFFs[i].index('.')], SpeciesNames[i])
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    SpGeneChromoCoord = ChromoGenesCoord(GFFs[i])
    # map each gene to its mRNA transcripts
    SpMapGeneTranscript = GeneToTranscripts(GFFs[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    SpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(SpGeneChromoCoord, SpMapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    SpGeneCoord = FromChromoCoordToGeneCoord(SpGeneChromoCoord)
    # Map Transcript names to gene names {transcript: gene}
    SpMapTranscriptGene = TranscriptToGene(GFFs[i])
    # get the coordinates of all exons    
    SpExonCoord = GeneExonCoord(GFFs[i])
    SpExonCoord = CleanGeneFeatureCoord(SpExonCoord, SpMapTranscriptGene)
    # get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
    SpIntronCoord = GeneIntronCoord(SpExonCoord)
    SpIntronCoord = CleanGeneFeatureCoord(SpIntronCoord, SpMapTranscriptGene)
    # get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
    SpTranscriptCoordinates = TranscriptsCoord(GFFs[i])
    # map genes to their longest transcript {gene: longest_transcript}
    SpGeneLongestTranscript = LongestTranscript(SpTranscriptCoordinates, SpMapGeneTranscript)
    # match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
    SpMatches = MatchHostTranscriptWithNestedTranscript(HostGenes[i], SpMapGeneTranscript, SpGeneLongestTranscript, SpTranscriptCoordinates, SpIntronCoord)
    # list all host, nested transcript pairs [[host, nested]]
    SpHostNestedPairs = GetHostNestedPairs(SpMatches)


# count nested intronless-same, intronless-opposite, intron-same, intron-opposite
CelIntronlessSame, CelIntronlessOpposite, CelIntronSame, CelIntronOpposite = 0, 0, 0, 0
# loop over host-nested transcript pairs
for i in range(len(CelHostNestedPairs)):
    # get the orientation of each transcript [+,+]
    OrientationPair = GenePairOrientation(CelHostNestedPairs[i], CelTranscriptCoordinates)
    # determine if the nested transcript has introns
    if CelHostNestedPairs[i][1] in CelIntronCoordinates:
        IntronLess = False
    else:
        IntronLess = True
    if IntronLess == True and len(set(OrientationPair)) == 1:
        # intronless, same orientation
        CelIntronlessSame += 1
    elif IntronLess == True and len(set(OrientationPair)) == 2:
        # intronless, differente orientation
        CelIntronlessOpposite += 1
    elif IntronLess == False and len(set(OrientationPair)) == 1:
        # with introns and same orientation
        CelIntronSame += 1
    elif IntronLess == False and len(set(OrientationPair)) == 2:
        # with introns and opposite orientation
        CelIntronOpposite += 1


# get the coordinates of briggsae genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
CbrGeneCoordChromo = ChromoGenesCoord(CbrGFF)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
CbrGeneCoord = FromChromoCoordToGeneCoord(CbrGeneCoordChromo)
# get briggsae expression {wormbase ID gene: [expression at 10 stages]}
CbrExpressionStage = ExpressionDevelopemtStages('WBPaper00041190.cbg.mr.csv')
# get the coordinates of introns # {chromo: {transcript:[(intron_start, intron_end), ...]}}
CbrIntronCoordChromo = CDSIntronCoord(CbrGFF, 'intron')
# create a dict {transcript: [(intron_start, intron_end), ...]}
CbrIntronCoordinates = {}
for i in CbrIntronCoordChromo:
    for j in CbrIntronCoordChromo[i]:
        CbrIntronCoordinates[j] = copy.deepcopy(CbrIntronCoordChromo[i][j])
# Map Transcript names to gene names {transcript: gene}
CbrMapTranscriptGene = TranscriptToGene(CbrGFF)
# List all transcripts for each gene {gene: [transcripts]}
CbrMapGeneTranscript = GeneToTranscripts(CbrGFF)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
CbrTranscriptCoordinates = TranscriptsCoord(CbrGFF)
# map genes to their longest transcript {gene: longest_transcript}
CbrGeneLongestTranscript = LongestTranscript(CbrTranscriptCoordinates, CbrMapGeneTranscript)
# match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
CbrMatches = MatchHostTranscriptWithNestedTranscript(CbrHostGenes, CbrMapGeneTranscript, CbrGeneLongestTranscript, CbrTranscriptCoordinates, CbrIntronCoordinates)
# list all host, nested transcript pairs [[host, nested]]
CbrHostNestedPairs = GetHostNestedPairs(CbrMatches)

# count nested intronless-same, intronless-opposite, intron-same, intron-opposite
CbrIntronlessSame, CbrIntronlessOpposite, CbrIntronSame, CbrIntronOpposite = 0, 0, 0, 0
# loop over host-nested transcript pairs
for i in range(len(CbrHostNestedPairs)):
    # get the orientation of each transcript [+,+]
    OrientationPair = GenePairOrientation(CbrHostNestedPairs[i], CbrTranscriptCoordinates)
    # determine if the nested transcript has introns
    if CbrHostNestedPairs[i][1] in CbrIntronCoordinates:
        IntronLess = False
    else:
        IntronLess = True
    if IntronLess == True and len(set(OrientationPair)) == 1:
        # intronless, same orientation
        CbrIntronlessSame += 1
    elif IntronLess == True and len(set(OrientationPair)) == 2:
        # intronless, differente orientation
        CbrIntronlessOpposite += 1
    elif IntronLess == False and len(set(OrientationPair)) == 1:
        # with introns and same orientation
        CbrIntronSame += 1
    elif IntronLess == False and len(set(OrientationPair)) == 2:
        # with introns and opposite orientation
        CbrIntronOpposite += 1


# get the coordinates of remanei gene on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
CrmGeneCoordChromo = CrmClaChromoGenesCoord(CrmGFF)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
CrmGeneCoord = FromChromoCoordToGeneCoord(CrmGeneCoordChromo)
# get remanei expression {wormbase ID gene: [expression at 10 stages]}
CrmExpression = ExpressionDevelopemtStages('WBPaper00041190.cre.mr.csv')
# get expression for genes with gene names e dict {CRE_gene_name: [expression values]}
CrmExpressionStage = MapCrmGeneNamesIDs(CrmExpression, 'c_remanei.PRJNA53967.WS256.geneIDs.txt')
# get the coordinates of introns {chromo: {transcript:[(intron_start, intron_end), ...]}}
CrmIntronCoordChromo = CrmClaCDSIntronCoord(CrmGFF, 'intron')
# create a dict {transcript: [(intron_start, intron_end), ...]}
CrmIntronCoordinates = {}
for i in CrmIntronCoordChromo:
    for j in CrmIntronCoordChromo[i]:
        CrmIntronCoordinates[j] = copy.deepcopy(CrmIntronCoordChromo[i][j])
# Map transcripts to their parent gene {transcript: gene}
CrmMapTranscriptGene = CrmClaTranscriptToGene(CrmGFF)
# Get the list of transcripts for each parent gene {gene: [transcripts]}
CrmMapGeneTranscript = CrmClaGeneToTranscripts(CrmMapTranscriptGene)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
CrmTranscriptCoordinates = CrmClaTranscriptCoord(CrmGFF)
# map genes to their longest transcript {gene: longest_transcript}
CrmGeneLongestTranscript = LongestTranscript(CrmTranscriptCoordinates, CrmMapGeneTranscript)
# match host transcripts with nested transcripts 
CrmMatches = CrmClaMatchHostTSNestedTS(CrmHostGenes, CrmTranscriptCoordinates, CrmIntronCoordinates)
# list all host, nested transcript pairs [[host, nested]]
CrmHostNestedPairs = GetHostNestedPairs(CrmHostGenes)

# count nested intronless-same, intronless-opposite, intron-same, intron-opposite
CrmIntronlessSame, CrmIntronlessOpposite, CrmIntronSame, CrmIntronOpposite = 0, 0, 0, 0
# loop over host-nested transcript pairs
for i in range(len(CrmHostNestedPairs)):
    # get the orientation of each transcript [+,+]
    OrientationPair = GenePairOrientation(CrmHostNestedPairs[i], CrmTranscriptCoordinates)
    # determine if the nested transcript has introns
    if CrmHostNestedPairs[i][1] in CrmIntronCoordinates:
        IntronLess = False
    else:
        IntronLess = True
    if IntronLess == True and len(set(OrientationPair)) == 1:
        # intronless, same orientation
        CrmIntronlessSame += 1
    elif IntronLess == True and len(set(OrientationPair)) == 2:
        # intronless, differente orientation
        CrmIntronlessOpposite += 1
    elif IntronLess == False and len(set(OrientationPair)) == 1:
        # with introns and same orientation
        CrmIntronSame += 1
    elif IntronLess == False and len(set(OrientationPair)) == 2:
        # with introns and opposite orientation
        CrmIntronOpposite += 1

# perform statistical tests to find out if there is an over-representation of intronless genes in same orientation
CelP = stats.fisher_exact([[CelIntronlessSame, CelIntronlessOpposite], [CelIntronSame, CelIntronOpposite]])[1]
CbrP = stats.fisher_exact([[CbrIntronlessSame, CbrIntronlessOpposite], [CbrIntronSame, CbrIntronOpposite]])[1]
CrmP = stats.fisher_exact([[CrmIntronlessSame, CrmIntronlessOpposite], [CrmIntronSame, CrmIntronOpposite]])[1]

# create table
# species, intronless-same, intronless-opposite, intron-same, intron-opposite, Pval
newfile = open('ContingencyTableIntronlessGenes.txt', 'w')
header = '\t'.join(['Species', 'Intronless-Same', 'Intronless-Opposite', 'Intron-Same', 'Intron-Opposite', 'P'])
newfile.write('\t'.join(['C. elegans', str(CelIntronlessSame), str(CelIntronlessOpposite), str(CelIntronSame), str(CelIntronOpposite), str(round(CelP, 4))])+'\n')
newfile.write('\t'.join(['C. briggsae', str(CbrIntronlessSame), str(CbrIntronlessOpposite), str(CbrIntronSame), str(CbrIntronOpposite), str(round(CbrP, 4))])+'\n')
newfile.write('\t'.join(['C. remanei', str(CrmIntronlessSame), str(CrmIntronlessOpposite), str(CrmIntronSame), str(CrmIntronOpposite), str(round(CrmP, 4))])+'\n')
newfile.close()