# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:53:18 2016

@author: RJovelin
"""

# use this script to classify overlapping genes into groups

# usage ClassifyOverlappingGenes.py [options]
# [human/chimp/gorilla]: species to consider

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


# get the species to consider from the command line
CurrSpecies = sys.argv[1]
assert CurrSpecies in ['human', 'chimp', 'gorilla'], 'species name is not valid'

# load dictionary of overlapping genes
if CurrSpecies == 'human':
    json_data = open('HumanOverlappingGenes.json')
elif CurrSpecies == 'chimp':
    json_data = open('ChimpOverlappingGenes.json')
elif CurrSpecies == 'gorilla':
    json_data = open('GorillaOverlappingGenes.json')
# load dictionary
OverlappingGenes = json.load(json_data)
json_data.close()

# make pairs of overlapping genes
OverlappingPairs = GetHostNestedPairs(OverlappingGenes)

# get GFF file
if CurrSpecies == 'human':
    GFF = 'Homo_sapiens.GRCh38.86.gff3'
elif CurrSpecies == 'chimp':
    GFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
elif CurrSpecies == 'gorilla':
    GFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
    
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
print('got gene coordinates on each chromosome')
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
print('mapped each gene to its mRNA transcripts', len(MapGeneTranscript))
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
print('removed genes lacking a mRNA transcript')
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
print('got gene coordinates', len(GeneCoord))
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
print('ordered genes on chromsomes')
# Map Transcript names to gene names {transcript: gene}
MapTranscriptGene = TranscriptToGene(GFF)
print('mapped transcripts to their parent gene', len(MapTranscriptGene))
# get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
ExonCoord = GeneExonCoord(GFF)
print('got exon coordinates', len(ExonCoord))
ExonCoord = CleanGeneFeatureCoord(ExonCoord, MapTranscriptGene)
print('cleaned up exon coordinates of non-mRNA transcripts', len(ExonCoord))
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
IntronCoord = GeneIntronCoord(ExonCoord)
print('got intron coordinates', len(IntronCoord))
IntronCoord = CleanGeneFeatureCoord(IntronCoord, MapTranscriptGene)
print('cleaned up intron coordinates of non-mRNA transcripts', len(IntronCoord))
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
CombinedIntronCoord = CombineAllGeneRegions(IntronCoord, MapTranscriptGene)
print('combined introns for each gene', len(CombinedIntronCoord))
# get the CDS coordinates of all transcripts {transcript: [[CDS_exon, CDS_end]]}
CDSCoord = GeneCDSCoord(GFF)
print('got CDS coordinates', len(CDSCoord))
CDSCoord = CleanGeneFeatureCoord(CDSCoord, MapTranscriptGene)
print('cleaned up CDS coordinates of non-mRNA transcripts', len(CDSCoord))


# create lists for groups of overlapping genes
Nested, SharingCDS, Interleaved, Piggyback, Convergent, Divergent = [], [], [], [], [], []




# loop over overlapping gene pairs
for pair in OverlappingPairs:
    for i in MapGeneTranscript[pair[0]]:
        for j in MapGeneTranscript[pair[1]]:
            for m in CDSCoord[i]:
                for k in CDSCoord[j]:
                    coord1 = set(range(m[0], m[1]))
                    coord2 = set(range(k[0], k[1]))
                    if len(coord1.intersection(coord2)) != 0:
                        SharingCDS.append([pair[0], pair[1]])



print(len(SharingCDS))
Shared = []
for i in SharingCDS:
    if i not in Shared:
        Shared.append(i)
print(len(Shared))
            
same = []
diff = []
            
for pair in Shared:
    assert GeneCoord[pair[0]][-1] in ['+', '-']
    if GeneCoord[pair[0]][-1] != GeneCoord[pair[1]][-1]:
        diff.append([pair[0], pair[1]])
    else:
        same.append([pair[0], pair[1]])            
            
            
            
            
print('same', len(same), same[0])
print('diff', len(diff), diff[0])









## identify itnronic nested genes {host_gene: [intronic_nested_gene]}
#HostGenes = FindIntronicNestedGenePairs(ContainedGenes, CombinedIntronCoord, GeneCoord)
#print('identified intronic nested genes', len(HostGenes))

    
