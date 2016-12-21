# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 13:48:04 2016

@author: RJovelin
"""


# use this script to plot expression divergence between external their
# un-nested orthologs and between internal and their un-nested orthologs

# usage python3 PlotExpDivergYoungOld.py 


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

# load dictionaries with overlapping genes
with open('HumanOverlappingGenes.json') as human_json_data:
    HumanOverlappingGenes = json.load(human_json_data)
# get human GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
# get the coordinates of human genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
HumanGeneChromoCoord = ChromoGenesCoord(HsaGFF)
# map each gene to its mRNA transcripts
HumanMapGeneTranscript = GeneToTranscripts(HsaGFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
HumanGeneChromoCoord = FilterOutGenesWithoutValidTranscript(HumanGeneChromoCoord, HumanMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
HumanGeneCoord = FromChromoCoordToGeneCoord(HumanGeneChromoCoord)

Counts = {}
Pairs = GetHostNestedPairs(HumanOverlappingGenes)
for i in Pairs:
    chromo = HumanGeneCoord[i[0]][0]
    if chromo in Counts:
        Counts[chromo] += 1
    else:
        Counts[chromo] = 1
LG = list(Counts.keys())
LG.sort()
for i in LG:
    print(i, Counts[i], len(HumanGeneChromoCoord[i]), (Counts[i]/ len(HumanGeneChromoCoord[i]))*100)
    
    
Genes = {}    
# open file to read
infile = open(HsaGFF, 'r')
for line in infile:
    line = line.rstrip()
    if 'gene' in line and not line.startswith('#'):
        line = line.split('\t')
        if line[2] == 'gene':
            # get biotype
            biotype = line[8][line[8].index('biotype=')+8: line[8].index(';', line[8].index('biotype=')+8)]
            # record only protein coding genes
            if biotype == 'protein_coding':
                # get the gene name
                gene = line[8][line[8].index('ID=gene:')+8: line[8].index(';')]
                # get chromo, start, end positions 0-based, and orientation
                chromo, start, end, sense = line[0], int(line[3]) -1, int(line[4]), line[6]            
                # check if chromo is recorded
                if gene not in Genes:
                    Genes[gene] = [[chromo, start, end, sense]]
                else:
                    Genes[gene].append([chromo, start, end, sense])
# close file after reading
infile.close()
total = 0
for gene in Genes:
    if len(Genes[gene]) > 1:
        total += 1
print(total)        