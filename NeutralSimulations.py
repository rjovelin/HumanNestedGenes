# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 11:38:22 2018

@author: RJovelin
"""

#  use this script to perform simulations and get expected number of overlap 
# and non-overlap genes under a model of randomization or gene extention-reduction

# usage python NeutralSimulations.py model species iterations

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

# get species name and number of iterations from command
model, species, iterations = sys.argv[1], sys.argv[2], int(sys.argv[3])


# make a dict of species names: [GFF file, assembly_ref] pairs
Species = {'Armadillo': ['Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Dasnov3.0'],
 'Cat': ['Felis_catus.Felis_catus_6.2.88.gff3', 'Felis_catus_6.2'],
 'Chimp': ['Pan_troglodytes.CHIMP2.1.4.88.gff3', 'CHIMP2.1.4'],
 'Cow': ['Bos_taurus.UMD3.1.88.gff3', 'UMD3.1'],
 'Dog': ['Canis_familiaris.CanFam3.1.88.gff3', 'CanFam3.1'],
 'Gorilla': ['Gorilla_gorilla.gorGor3.1.88.gff3', 'gorGor3.1'],
 'Hedgehog': ['Erinaceus_europaeus.HEDGEHOG.88.gff3', 'HEDGEHOG'],
 'Horse': ['Equus_caballus.EquCab2.88.gff3', 'EquCab2'],
 'Human': ['Homo_sapiens.GRCh38.88.gff3', 'GRCh38'],
 'Macaque': ['Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Mmul_8.0.1'],
 'Marmoset': ['Callithrix_jacchus.C_jacchus3.2.1.88.gff3', 'C_jacchus3.2.1'],
 'Mouse': ['Mus_musculus.GRCm38.88.gff3', 'GRCm38'],
 'Opossum': ['Monodelphis_domestica.BROADO5.88.gff3', 'BROADO5'],
 'Orangutan': ['Pongo_abelii.PPYG2.88.gff3', 'PPYG2'],
 'Platypus': ['Ornithorhynchus_anatinus.OANA5.88.gff3', 'OANA5'],
 'Shrew': ['Sorex_araneus.COMMON_SHREW1.88.gff3', 'COMMON_SHREW1'],
 'Sloth': ['Choloepus_hoffmanni.choHof1.88.gff3', 'choHof1']}

# get GFF and assembly_ref
GFF = Species[species][0]
AssemblyRef = Species[species][1]

# get the coordinates of genes on each chromo {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
# get chromosome length
ChromoL = GetChromosomeLength(GFF, GeneChromoCoord, AssemblyRef)
# do some QC
# make a set of chromo
a, b = set(ChromoL.keys()), set(GeneChromoCoord.keys())
assert a == b

# create a dict to store the count of overlapping and non-overlapping genes for each simulation
Simulations = {}

# perform N simulations for given species    
for i in range(iterations):
    # check model
    if model == 'randomization':
        # randomize gene positions along each chromo
        RandomCoord = RandomizeGenePosition(GeneChromoCoord, ChromoL)
    elif model == 'extension':
        # randomize gene coordinates by allowing gene extension and gene reduction
        RandomCoord = RandomizeGeneLength(GeneChromoCoord)    
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    RandomGeneCoord = FromChromoCoordToGeneCoord(RandomCoord)
    # order randomized genes on each chromo
    OrderedRandom = OrderGenesAlongChromo(RandomCoord)
    # find overlapping genes from randomized genes
    OverlapRandom = FindOverlappingGenePairs(RandomCoord, OrderedRandom)
    # make a set of overlapping genes
    ovlp = MakeFullPartialOverlapGeneSet(OverlapRandom)
    # make a set of non-overlapping gene
    nonovlp = MakeNonOverlappingGeneSet(OverlapRandom, RandomGeneCoord)    
    # record number of overlapping and non-overlapping genes
    Simulations[species + '_' + model + '_' + str(i)] = [len(ovlp), len(nonovlp)]
    

# save results of simulation as json file
newfile = open('Simulations_' + species.title() + '_' + model.title() + '.json', 'w')  
# save dictionary to file        
json.dump(Simulations, newfile, sort_keys = True, indent = 4)
newfile.close()
