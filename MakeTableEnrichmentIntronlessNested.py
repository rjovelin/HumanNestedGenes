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


# create lists to store the number of intronless nested genes and nested genes
# with introns that are on the same strand and opposite strand with their hosts
# [human, chimp, gorilla, orang-outan, macaque]
IntronlessSame, IntronlessOpposite, IntronSame, IntronOpposite = [], [], [], []


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
    NoIntronSame, NoIntronOpposite, WithIntronSame, WithIntronOpposite = 0, 0, 0, 0
    # loop over host-nested transcript pairs
    for i in range(len(SpHostNestedPairs)):
        # get the orientation of each transcript [+,+]
        OrientationPair = GenePairOrientation(SpHostNestedPairs[i], SpTranscriptCoordinates)
        # determine if the nested transcript has introns
        if SpHostNestedPairs[i][1] in SpIntronCoord:
            IntronLess = False
        else:
            IntronLess = True
        # check presence of intron and orientation of host and nested transcripts
        if IntronLess == True and len(set(OrientationPair)) == 1:
            # intronless, same orientation
            NoIntronSame += 1
        elif IntronLess == True and len(set(OrientationPair)) == 2:
            # intronless, differente orientation
            NoIntronOpposite += 1
        elif IntronLess == False and len(set(OrientationPair)) == 1:
            # with introns and same orientation
            WithIntronSame += 1
        elif IntronLess == False and len(set(OrientationPair)) == 2:
            # with introns and opposite orientation
            WithIntronOpposite += 1
    # populate lists
    IntronlessSame.append(NoIntronSame)
    IntronlessOpposite.append(NoIntronOpposite)
    IntronSame.append(WithIntronSame)
    IntronOpposite.append(WithIntronOpposite)      

# create a list to store the P values, in the same order as the lists with number of nested genes
PVal = []
# perform statistical tests to find out if there is an over-representation of intronless genes in same orientation
for i in range(len(SpeciesNames)):
    P = stats.fisher_exact([[IntronlessSame[i], IntronlessOpposite[i]], [IntronSame[i], IntronOpposite[i]]])[1]
    PVal.append(P)

# create table
# species, intronless-same, intronless-opposite, intron-same, intron-opposite, Pval
newfile = open('ContingencyTableIntronlessGenes.txt', 'w')
header = '\t'.join(['Species', 'Intronless-Same', 'Intronless-Opposite', 'Intron-Same', 'Intron-Opposite', 'P'])
newfile.write(header + '\n')

# loop over species:
for i in range(len(SpeciesNames)):
    line = '\t'.join([SpeciesNames[i], str(IntronlessSame[i]), str(IntronlessOpposite[i]), str(IntronSame[i]), str(IntronOpposite[i]), str(PVal[i])])
    newfile.write(line + '\n')
newfile.close()