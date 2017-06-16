# -*- coding: utf-8 -*-
"""
Created on Thu May 11 13:06:34 2017

@author: RJovelin
"""

# use this script to make a table with counts of overlapping pairs and genes in each species

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


# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Macaque', 'Marmoset',
           'Hedgehog', 'Shrew', 'Cat', 'Dog', 'Mouse', 'Cow', 'Horse', 'Sloth',
           'Armadillo', 'Opossum', 'Platypus']

# make a list of GFF files
GFF = ['Homo_sapiens.GRCh38.88.gff3', 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
       'Gorilla_gorilla.gorGor3.1.88.gff3', 'Pongo_abelii.PPYG2.88.gff3',
       'Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
       'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Sorex_araneus.COMMON_SHREW1.88.gff3',
       'Felis_catus.Felis_catus_6.2.88.gff3', 'Canis_familiaris.CanFam3.1.88.gff3',
       'Mus_musculus.GRCm38.88.gff3', 'Bos_taurus.UMD3.1.88.gff3', 
       'Equus_caballus.EquCab2.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
       'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Monodelphis_domestica.BROADO5.88.gff3',
       'Ornithorhynchus_anatinus.OANA5.88.gff3']


# create a dict with gene counts for all species
GeneCounts = {}
# make a list with types of gene overlaps
GeneTypes = ['Overlapping', 'Nested', 'PiggyBack', 'Convergent', 'Divergent']

# loop over species
for i in range(len(Species)):
    GeneCounts[Species[i]] = []
    # count the number of valid protein coding genes
    # get the coordinates of genes on each chromo {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # add total number of genes to list  
    total = len(set(GeneCoord.keys()))
    GeneCounts[Species[i]].append(total)    
    # loop over type of overlapping genes
    for j in range(len(GeneTypes)):
        # load dictionary
        json_data = open(Species[i] + GeneTypes[j] + 'Genes.json')
        overlapping = json.load(json_data)
        json_data.close()
        # count the number of gene pairs 
        Pairs = len(GetHostNestedPairs(overlapping))
        # count the number of overlapping genes
        Genes = len(MakeFullPartialOverlapGeneSet(overlapping))
        # compute percent of overlapping genes
        percent = round(Genes / total * 100, 2)
        # populate dict
        for item in [Pairs, Genes, percent]:
            GeneCounts[Species[i]].append(item)

# save counts to file
newfile = open('OverlappingGeneCounts.txt', 'w')
# build headers
Header1 = ['Species', 'AllGenes', '\t', 'All', '\t', '\t', 'Nst', '\t', '\t', 'Pbk', '\t', '\t', 'Con', '\t', '\t', 'Div']
Header2 = ['\t', '\t']
for i in range(5):
    Header2.expend(['Pairs', 'Genes', '(%)'])
# write headers to file, tab seperated
newfile.write('\t'.join(Header1) + '\n')
newfile.write('\t'.join(Header2) + '\n')
# write counts to file for each species, keep the phylogenetic order
for i in range(len(Species)):
    line = [Species[i]] + list(map(lambda x: str(x), GeneCounts[Species[i]]))
    newfile.write('\t'.join(line) + '\n')
newfile.close()

