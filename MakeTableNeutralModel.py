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


# load dictionary of counts of overlapping and non-overlapping genes from simulations
infile = open('SimulationsOverlap.json')  
Simulations = json.load(infile)
infile.close()

# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan',
           'Macaque', 'Marmoset', 'Hedgehog', 'Shrew',
           'Cat', 'Dog', 'Mouse', 'Cow', 'Horse', 'Sloth',
           'Armadillo', 'Opossum', 'Platypus']
# make a list a GFF files
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
for i in range(len(GFF)):
    # get the coordinates of genes on each chromo {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    SpGeneCoordinates.append(GeneCoord)    
    
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

# make a dict with mean and standard error from simulations
Mean = {}, SEM = {}
for species in Simulations:
    NO = [Simulations[species][0] for species in Simulations]
    NN = [Simulations[species][1] for species in Simulations]
    Mean[species] = [np.mean(NO), np.mean(NN)]    
    SEM[species] = [np.std(NO) / math.sqrt(len(NO)), np.std(NN) / math.sqrt(len(NN))]     
        
# compute P values of differences between observations and expectations
PVals = {}
for species in ObservedOverlaGenes:
    P = stats.fisher_exact([ObservedOverlaGenes[species], Mean[species]])[1]     
    PVals[species] = P
    
for species in ObservedOverlaGenes:
    print(species, ObservedOverlaGenes[species][0], ObservedOverlaGenes[species][1], Mean[species][0], Mean[species][1], PVals[species])
    
    
    
    
    
    
    