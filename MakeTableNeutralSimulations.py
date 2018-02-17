# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 13:00:01 2018

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


# make a list of simulation runs for all species for randomization model
ShufflingJson = [i for i in os.listdir() if 'Simulations_' in i and '_Randomization.json']
# make a list of simulation runs for all species for gene extension model
GeneExtensionJson = [i for i in os.listdir() if 'Simulations_' in i and '_Extension.json']
JsonFiles = [ShufflingJson, GeneExtensionJson]

# make a dict for each model {species: [json simulations runs]}
Shuffling, GeneExtension = {}, {}

# get the gene counts for each run for each species and model
for i in range(len(JsonFiles)):
    for filename in JsonFiles[i]:
        # get the species name
        species = filename[filename.index('_')+1: filename.find('_', filename.index('_')+1)]
        # load json dict
        infile = open(filename)
        runs = json.load(infile)
        infile.close()
        # get counts of ovlp and non-ovl genes for each run for that species
        if i == 0:
            Shuffling[species] = [runs[k] for k in runs]
        elif i == 1:
            GeneExtension[species] = [runs[k] for k in runs]


# make a dict of {species names: GFF file} pairs
GFF = {'Armadillo': 'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 
       'Cat': 'Felis_catus.Felis_catus_6.2.88.gff3',
       'Chimp': 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
       'Cow': 'Bos_taurus.UMD3.1.88.gff3',
       'Dog': 'Canis_familiaris.CanFam3.1.88.gff3',
       'Gorilla': 'Gorilla_gorilla.gorGor3.1.88.gff3',
       'Hedgehog': 'Erinaceus_europaeus.HEDGEHOG.88.gff3',
       'Horse': 'Equus_caballus.EquCab2.88.gff3',
       'Human': 'Homo_sapiens.GRCh38.88.gff3',
       'Macaque': 'Macaca_mulatta.Mmul_8.0.1.88.gff3',
       'Marmoset': 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
       'Mouse': 'Mus_musculus.GRCm38.88.gff3',
       'Opossum': 'Monodelphis_domestica.BROADO5.88.gff3',
       'Orangutan': 'Pongo_abelii.PPYG2.88.gff3',
       'Platypus': 'Ornithorhynchus_anatinus.OANA5.88.gff3',
       'Shrew': 'Sorex_araneus.COMMON_SHREW1.88.gff3',
       'Sloth': 'Choloepus_hoffmanni.choHof1.88.gff3'}

# make an ordered list of species names
SpeciesNames = ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Macaque',
                'Marmoset', 'Mouse', 'Cat', 'Dog', 'Cow', 'Horse', 'Hedgehog',
                'Shrew', 'Sloth', 'Armadillo', 'Opossum', 'Platypus']


# make a parallel list of gene coordinates in each species {gene:[chromosome, start, end, sense]}
SpGeneCoordinates = []
for species in SpeciesNames:
    # get the coordinates of genes on each chromo {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[species])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[species])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    SpGeneCoordinates.append(GeneCoord)    
    
# count the number of observed overlapping and non-overlapping genes in each species {sp: [overlap, non-overlap]}
ObservedOvlGenes = {}
# loop over species
for i in range(len(SpeciesNames)):
    # load dict of overlapping genes
    json_data = open(SpeciesNames[i] + 'Overlapping' + 'Genes.json')
    overlapping = json.load(json_data)
    json_data.close()
    # count the number of overlapping genes, store in list
    ObservedOvlGenes[SpeciesNames[i]] = [len(MakeFullPartialOverlapGeneSet(overlapping))]
    # count the number of non-overlapping genes, store in list
    ObservedOvlGenes[SpeciesNames[i]].append(len(MakeNonOverlappingGeneSet(overlapping, SpGeneCoordinates[i])))    


# compare the number of observed and simulated overlapping genes 
ResultsShuffle, ResultsGeneExtention = {}, {}

# loop over species
for i in range(2):
    for species in SpeciesNames:
        # count the mean number of overlapping and non-overlapping genes from all runs
        if i == 0:
            OvlNum = [item[0] for item in Shuffling[species]]
            NonOvlNum = [item[1] for item in Shuffling[species]]
        elif i == 1:
            OvlNum = [item[0] for item in GeneExtension[species]]
            NonOvlNum = [item[0] for item in GeneExtension[species]]
        SEMOvlNum = np.std(OvlNum) / math.sqrt(len(OvlNum))
        SEMNonOvlNum = np.std(NonOvlNum) / math.sqrt(len(NonOvlNum))
        # compute P values of differences between observations and expectations
        P = stats.fisher([ObservedOvlGenes[species], [round(np.mean(OvlNum)), round(np.mean(NonOvlNum))]])[1]
        if i == 0:
            ResultsShuffle[species] = [np.mean(OvlNum), SEMOvlNum, np.mean(NonOvlNum), SEMNonOvlNum, P]
        elif i == 1:
            ResultsGeneExtention[species] = [np.mean(OvlNum), SEMOvlNum, np.mean(NonOvlNum), SEMNonOvlNum, P]

# write results to file
newfile = open('NeutralSimulations.txt', 'w')
Header = ['Species', 'Obs_Ovl', 'Obs_NonOvl', 'Shuffle_Ovl (SEM)', 'Shuffle_NonOvl (SEM)', 'P', 'Extension_Ovl (SEM), Extension_NonOvl (SEM)', 'P']    
newfile.write('\t'.join(Header) + '\n')
for species in SpeciesNames:
    # make the line
    line = [species, ObservedOvlGenes[species][0], ObservedOvlGenes[species][1],
            '{0} ({1})'.format(ResultsShuffle[species][0], ResultsShuffle[species][1]),
            '{0} ({1})'.format(ResultsShuffle[species][2], ResultsShuffle[species][3]),
            ResultsShuffle[species][4],
             + '{0} ({1})'.format(ResultsGeneExtention[species][0], ResultsGeneExtention[species][1]),
             + '{0} ({1})'.format(ResultsGeneExtention[species][2], ResultsGeneExtention[species][3]),
             ResultsGeneExtention[species][4]]
    newfile.write('\t'.join(list(map(lambda x: str(x), line))) + '\n')
newfile.close()


    
    
    