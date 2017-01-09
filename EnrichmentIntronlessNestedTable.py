# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 11:57:29 2016

@author: RJovelin
"""


# use this script to test for over-representation of intronless nested genes in same sense orientation
# generate a contagency table with number of nested intronless-same sense, intronless opposite sense,
# with introns same sense, and with introns opposite sense

# usage python3 EnrichmentIntronlessNestedTable.py


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
with open('HumanNestedGenes.json') as human_json_data:
    Nested = json.load(human_json_data)

# get GFF file
GFF = 'Homo_sapiens.GRCh38.86.gff3'
 
# count internal (nested) intronless genes, genes with introns that are on the
# same strand and opposite strand relative to their hosts
NoIntronSame, NoIntronOpposite, WithIntronSame, WithIntronOpposite = 0, 0, 0, 0

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

# loop over host-nested transcript pairs
for i in range(len(HostNestedPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedPairs[i][0] in MapTranscriptGene
    assert HostNestedPairs[i][0] in TranscriptCoordinates    
    assert HostNestedPairs[i][1] in MapTranscriptGene
    assert HostNestedPairs[i][1] in TranscriptCoordinates    
    # get the orientation of each transcript [+,+]
    OrientationPair = GenePairOrientation(HostNestedPairs[i], TranscriptCoordinates)
    # determine if the nested transcript has introns
    if HostNestedPairs[i][1] in IntronCoord:
        IntronLess = False
    else:
        IntronLess = True
    # check presence of intron and orientation of host and nested transcripts
    if IntronLess == True and len(set(OrientationPair)) == 1:
        # intronless, same orientation
        NoIntronSame += 1
        assert set(OrientationPair) == {'-', '-'} or set(OrientationPair) == {'+', '+'}
    elif IntronLess == True and len(set(OrientationPair)) == 2:
        # intronless, differente orientation
        NoIntronOpposite += 1
        assert set(OrientationPair) == {'-', '+'} 
    elif IntronLess == False and len(set(OrientationPair)) == 1:
        # with introns and same orientation
        WithIntronSame += 1
        assert set(OrientationPair) == {'-', '-'} or set(OrientationPair) == {'+', '+'}
    elif IntronLess == False and len(set(OrientationPair)) == 2:
        # with introns and opposite orientation
        WithIntronOpposite += 1
        assert set(OrientationPair) == {'-', '+'}


# test that the proportions of intron-less and with-intron genes are the same on both strands
P = stats.fisher_exact([[NoIntronSame, NoIntronOpposite], [WithIntronSame, WithIntronOpposite]])[1]

# open file for writing results
newfile = open('ContingencyTableIntronlessGenes.txt', 'w')



# create contingency table table
newfile.write('Table 1. Number of internal (nested genes) with and without introns and their orientation\n')
newfile.write('\t'.join(['', 'Intronless', 'WithIntron', 'Ratio Intronless/Total', 'P']) + '\n')
newfile.write('\t'.join(['Same', str(NoIntronSame), str(WithIntronSame), str(round(NoIntronSame / (NoIntronSame + WithIntronSame), 4)), str(P)]) + '\n')
newfile.write('\t'.join(['Opposite', str(NoIntronOpposite), str(WithIntronOpposite), str(round(NoIntronOpposite / (NoIntronOpposite + WithIntronOpposite), 4)), str(P)]) + '\n')


newfile.close()


## make a set of un-nested genes
#UnNestedGenes = MakeNonOverlappingGeneSet(HostGenes, GeneCoord)



