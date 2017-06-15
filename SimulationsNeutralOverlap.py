# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 18:14:00 2017

@author: RJovelin
"""

#  use this script to perform simulations and get expected number of overlap and non-overlap genes 

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


# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan',
           'Macaque', 'Marmoset', 'Hedgehog', 'Shrew',
           'Cat', 'Dog', 'Mouse', 'Cow', 'Horse', 'Sloth',
           'Armadillo', 'Opossum', 'Platypus']

GFF = ['Homo_sapiens.GRCh38.88.gff3', 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
       'Gorilla_gorilla.gorGor3.1.88.gff3', 'Pongo_abelii.PPYG2.88.gff3',
       'Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
       'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Sorex_araneus.COMMON_SHREW1.88.gff3',
       'Felis_catus.Felis_catus_6.2.88.gff3', 'Canis_familiaris.CanFam3.1.88.gff3',
       'Mus_musculus.GRCm38.88.gff3', 'Bos_taurus.UMD3.1.88.gff3', 
       'Equus_caballus.EquCab2.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
       'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Monodelphis_domestica.BROADO5.88.gff3',
       'Ornithorhynchus_anatinus.OANA5.88.gff3']

# make a list with gene coordinates in each species {gene:[chromosome, start, end, sense]}
SpGeneCoordinates = []
# make a list with gene corodinates on each chromo in each species {chromo: {gene:[chromosome, start, end, sense]}}
SpChromoGeneCoord = []
# make a list with chromosome length for each species
SpChromoLength = []

# populate lists
for i in range(len(GFF)):
    # get the coordinates of genes on each chromo {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    SpChromoGeneCoord.append(GeneChromoCoord)    
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    SpGeneCoordinates.append(GeneCoord)    
    # odrder genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    # get chromosome length
    ChromoL = GetChromosomeLength(GFF[0])
    SpChromoLength.append(ChromoL)
    
# create a dict to store the count of overlapping and non-overlapping genes for each simulation
Simulations = {}
for species in Species:
    Simulations[species] = []     

# perform 10 simulations for each species    
for i in range(10):
    # loop over species
    for j in range(len(Species)):
        print(i, Species[j])
        # get the gene coordinates on each chromo {chromo: {gene:[chromosome, start, end, sense]}}
        GeneChromoCoord = SpChromoGeneCoord[j]
        # get the length of each chromosome
        ChromoLength = SpChromoLength[j]
        # randomize gene positions along each chromo
        RandomCoord = RandomizeGenePosition(GeneChromoCoord, ChromoLength)
        print('done with randomization')
        # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
        RandomGeneCoord = FromChromoCoordToGeneCoord(RandomCoord)
        # order randomized genes on each chromo
        OrderedRandom = OrderGenesAlongChromo(RandomCoord)
        print('ordered genes')
        # find overlapping genes from randomized genes
        OverlapRandom = FindOverlappingGenePairs(RandomCoord, OrderedRandom)
        print('found overlapping genes')   
        # make a set of overlapping genes
        ovlp = MakeFullPartialOverlapGeneSet(OverlapRandom)
        # make a set of non-overlapping gene
        nonovlp = MakeNonOverlappingGeneSet(OverlapRandom, RandomGeneCoord)    
        print(i, Species[j], len(ovlp), len(nonovlp))
        # record number of overlapping and non-overlapping genes
        Simulations[Species[j]].append(len(ovlp), len(nonovlp))

# save results of simulation as json file
newfile = open('SimulationsOverlap.json', 'w')  
# save dictionary to file        
json.dump(Simulations, newfile, sort_keys = True, indent = 4)
newfile.close()