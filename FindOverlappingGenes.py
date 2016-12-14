# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 16:25:38 2016

@author: RJovelin
"""

# use this script to identify overlapping genes in human, chimp and gorilla
# save dictionaries of overlapping genes as json files

# usage python3 FindOverlappingGenes.py


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
    # save dictionaries as json file
    if i == 0:
        # save overlapping genes as json file
        newfile = open('HumanOverlappingGenes.json', 'w')
        json.dump(SpOverlappingGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
    elif i == 1:
        # save overlapping genes as json file
        newfile = open('ChimpOverlappingGenes.json', 'w')
        json.dump(SpOverlappingGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
    elif i == 2:
        # save overlapping genes as json file
        newfile = open('GorillaOverlappingGenes.json', 'w')
        json.dump(SpOverlappingGenes, newfile, sort_keys = True, indent = 4)
        newfile.close()
        

