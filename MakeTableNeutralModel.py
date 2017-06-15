# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 13:51:07 2017

@author: RJovelin
"""

#  use this script to compare observed number of overlapping genes to a null neutral model

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
# make a list with ordered genes in each chromo
SpOrderedGenes = []
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
    SpOrderedGenes.append(OrderedGenes)
    # get chromosome length
    ChromoL = GetChromosomeLength(GFF[0])
    SpChromoLength.append(ChromoL)
    
    
# count the number of observed overlapping and non-overlapping genes in each species {sp: [overlap, non-overlap]}
ObservedOverlaGenes = {}
# loop over species
for i in range(len(Species)):
    # load dict of overlapping genes
    json_data = open(Species[i] + 'Overlapping' + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    # count the number of overlapping genes, store in list
    ObservedOverlaGenes[Species[i]] = [len(MakeFullPartialOverlapGeneSet(overlapping))]
    # count the number of non-overlapping genes, store in list
    ObservedOverlaGenes[Species[i]].append(len(MakeNonOverlappingGeneSet(overlapping, SpGeneCoordinates[i])))    

for i in range(len(Species)):
    print(Species[i], ObservedOverlaGenes[Species[i]])
    
    
    
for i in range(3):
    # get the coordinates of human genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[0])
    # get the length of each chromosome
    ChromoLength = GetChromosomeLength(GFF[0])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[0])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # order genes along chromosomes    
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    # randomize gene positions along each chromo
    RandomCoord = RandomizeGenePosition(ChromoGeneCoord, ChromoLength)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    RandomGeneCoord = FromChromoCoordToGeneCoord(RandomCoord)    
    # order randomized genes on each chromo
    OrderedRandom = OrderGenesAlongChromo(RandomCoord)
    # find overlapping genes from randomized genes
    OverlapRandom = FindOverlappingGenePairs(RandomCoord, OrderedRandom)
    # make a set of overlapping genes
    ovlp = MakeFullPartialOverlapGeneSet(OverlapRandom)
    # make a set of non-overlapping gene
    nonovlp = MakeNonOverlappingGeneSet(ovlp, RandomGeneCoord)    
    print(i, len(ovlp), len(nonovlp))

  


