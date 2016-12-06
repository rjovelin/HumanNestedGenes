# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 16:25:38 2016

@author: RJovelin
"""

# use this script to identify overlapping genes, nested and itnronic nested genes
# in human, chimp and gorilla
# save dictionaries of overlapping and hosy:nested genes as json files

# usage python3 FindNestedGenes.py


import os
import sys
import numpy as np
import random
import json
import math
from HsaNestedGenes import *


# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
GgoGFF = 'Gorilla_gorilla.gorGor3.1.86.gff3'
    
# make a list of primate GFF files
GFFs = [HsaGFF, PtrGFF, GgoGFF]

# loop over GFF files, find nested and intronic=nested genes in each species 
for i in range(len(GFFs)):
    print(GFFs[i][:GFFs[i].index('.')])
    # get the coordinates of human genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    SpGeneChromoCoord = ChromoGenesCoord(GFFs[i])
    print('got gene coordinates on each chromosome')
    # map each gene to its mRNA transcripts
    SpMapGeneTranscript = GeneToTranscripts(GFFs[i])
    print('mapped each gene to its mRNA transcripts', len(SpMapGeneTranscript))
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    SpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(SpGeneChromoCoord, SpMapGeneTranscript)
    print('removed genes lacking a mRNA transcript')
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    SpGeneCoord = FromChromoCoordToGeneCoord(SpGeneChromoCoord)
    print('got gene coordinates', len(SpGeneCoord))
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    SpOrderedGenes = OrderGenesAlongChromo(SpGeneChromoCoord)
    print('ordered genes on chromsomes')
    # Find overlapping genes {gene1: [gene2, gene3]}
    SpOverlappingGenes = FindOverlappingGenePairs(SpGeneChromoCoord, SpOrderedGenes)
    print('found overlapping genes', len(SpOverlappingGenes))
    # Find genes fully contained in another gene {containing: [contained1, contained2]}
    SpContainedGenes = FindContainedGenePairs(SpGeneCoord, SpOverlappingGenes)
    print('found genes contained in other genes', len(SpContainedGenes))
    # Map Transcript names to gene names {transcript: gene}
    SpMapTranscriptGene = TranscriptToGene(GFFs[i])
    print('mapped transcripts to their parent gene', len(SpMapTranscriptGene))
    # get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
    SpExonCoord = GeneExonCoord(GFFs[i])
    print('got exon coordinates', len(SpExonCoord))
    SpExonCoord = CleanGeneFeatureCoord(SpExonCoord, SpMapTranscriptGene)
    print('cleaned up exon coordinates of non-mRNA transcripts', len(SpExonCoord))
    # get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
    SpIntronCoord = GeneIntronCoord(SpExonCoord)
    print('got intron coordinates', len(SpIntronCoord))
    SpIntronCoord = CleanGeneFeatureCoord(SpIntronCoord, SpMapTranscriptGene)
    print('cleaned up intron coordinates of non-mRNA transcripts', len(SpIntronCoord))
    # Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
    SpCombinedIntronCoord = CombineAllGeneRegions(SpIntronCoord, SpMapTranscriptGene)
    print('combined introns for each gene', len(SpCombinedIntronCoord))
    # identify itnronic nested genes {host_gene: [intronic_nested_gene]}
    SpHostGenes = FindIntronicNestedGenePairs(SpContainedGenes, SpCombinedIntronCoord, SpGeneCoord)
    print('identified intronic nested genes', len(SpHostGenes))

    
    if i == 0:
        # save contained genes as json file
        newfile = open('HumanContainedGenes.json', 'w')
        json.dump(SpContainedGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        # save intronic nested genes as json file
        newfile = open('HumanHostNestedGenes.json', 'w')
        json.dump(SpHostGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        # save overlapping genes as json file
        newfile = open('HumanOverlappingGenes.json', 'w')
        json.dump(SpOverlappingGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
    elif i == 1:
        # save contained genes as json file
        newfile = open('ChimpContainedGenes.json', 'w')
        json.dump(SpContainedGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        # save intronic nested genes as json file
        newfile = open('ChimpHostNestedGenes.json', 'w')
        json.dump(SpHostGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        # save overlapping genes as json file
        newfile = open('ChimpOverlappingGenes.json', 'w')
        json.dump(SpOverlappingGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
    elif i == 2:
        newfile = open('GorillaContainedGenes.json', 'w')
        json.dump(SpContainedGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        # save intronic nested genes as json file
        newfile = open('GorillaHostNestedGenes.json', 'w')
        json.dump(SpHostGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        # save overlapping genes as json file
        newfile = open('GorillaOverlappingGenes.json', 'w')
        json.dump(SpOverlappingGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        

