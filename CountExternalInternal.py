# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 17:39:11 2017

@author: RJovelin
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


# map genes to their longest transcript {gene: longest_transcript}
GeneLongestTranscript = LongestTranscript(TranscriptCoordinates, MapGeneTranscript)
# match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
Matches = MatchHostTranscriptWithNestedTranscript(Nested, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# list all host, nested transcript pairs [[host, nested]]
HostNestedPairs = GetHostNestedPairs(Matches)

# make pairs of host-nested genes with intronless and intron-containing internal genes
ExtWithIntrons, ExtNoIntrons = set(), set()
IntWithIntrons, IntNoIntrons = set(), set()
PairsWithIntrons, PairsNoIntrons = [], []
# loop over host-nested transcript pairs
for i in range(len(HostNestedPairs)):
    # check that both transcripts have coordinates and have corresponding gene names    
    assert HostNestedPairs[i][0] in MapTranscriptGene and HostNestedPairs[i][0] in TranscriptCoordinates    
    assert HostNestedPairs[i][1] in MapTranscriptGene and HostNestedPairs[i][1] in TranscriptCoordinates    
    if HostNestedPairs[i][1] in IntronCoord:
        # internal gene has introns, add gene pairs to list 
        PairsWithIntrons.append([MapTranscriptGene[HostNestedPairs[i][0]], MapTranscriptGene[HostNestedPairs[i][1]]])
        ExtWithIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
        IntWithIntrons.add(MapTranscriptGene[HostNestedPairs[i][1]])
    else:
        # internal gene is intronless, add gene pair to list
        PairsNoIntrons.append([MapTranscriptGene[HostNestedPairs[i][0]], MapTranscriptGene[HostNestedPairs[i][1]]])
        ExtNoIntrons.add(MapTranscriptGene[HostNestedPairs[i][0]])
        IntNoIntrons.add(MapTranscriptGene[HostNestedPairs[i][1]])
        
print('N pairs intronless', len(PairsNoIntrons))
print('N pairs with introns', len(PairsWithIntrons))
print('N external with no intron gene', len(ExtNoIntrons))
print('N intronless internal', len(IntNoIntrons))
print('N external with intron-containing genes', len(ExtWithIntrons))
print('N intron-containing internal', len(IntWithIntrons))

