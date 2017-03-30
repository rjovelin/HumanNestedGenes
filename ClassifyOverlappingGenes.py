# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:53:18 2016

@author: RJovelin
"""

# use this script to classify overlapping genes into 4 groups

# usage ClassifyOverlappingGenes.py 


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



# make a list of GFF files
GFF_Files = ['Bos_taurus.UMD3.1.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
             'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
             'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
             'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
             'Gorilla_gorilla.gorGor3.1.88.gff3', 'Homo_sapiens.GRCh38.88.gff3',
             'Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Monodelphis_domestica.BROADO5.88.gff3',
             'Mus_musculus.GRCm38.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3', 
             'Pan_troglodytes.CHIMP2.1.4.88.gff3', 'Pongo_abelii.PPYG2.88.gff3', 
             'Sorex_araneus.COMMON_SHREW1.88.gff3']

# make a parallel list of Species names
Species = ['Cow', 'Marmoset', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
           'Cat', 'Gorilla', 'Human', 'Macaque', 'Opossum', 'Mouse', 'Platypus',
           'Chimp', 'Orangutan', 'Shrew']

# make a parallel list of json files with overlapping genes
JsonFiles = [i + 'OverlappingGenes.json' for i in Species]


# loop over json files
for i in range(len(JsonFiles)):
    json_data = open(JsonFiles[i])
    # load dictionary
    OverlappingGenes = json.load(json_data)
    json_data.close()
    # make pairs of overlapping genes
    OverlappingPairs = GetHostNestedPairs(OverlappingGenes)
    # get GFF file
    GFF = GFF_Files[i]
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF)
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF)
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    # Map Transcript names to gene names {transcript: gene}
    MapTranscriptGene = TranscriptToGene(GFF)
    # get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
    ExonCoord = GeneExonCoord(GFF)
    ExonCoord = CleanGeneFeatureCoord(ExonCoord, MapTranscriptGene)
    # get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
    IntronCoord = GeneIntronCoord(ExonCoord)
    IntronCoord = CleanGeneFeatureCoord(IntronCoord, MapTranscriptGene)
    # Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
    CombinedIntronCoord = CombineAllGeneRegions(IntronCoord, MapTranscriptGene)


    # define 4 categories of overlapping genes:
    # 1) nested gene pairs: one gene is fully contained within the intron of another gene
    #    both genes can be on the same strand or on opposite strands
    # 2) piggyback gene pairs: both gene have the same orientation
    #    overlap can be partial or complete
    # 3) convergent gene pairs: both genes have different orientation
    #    first gene on chromosome is +, second gene on chromo is -
    # 4) divergent gene pairs: both genes have different orientation
    #    first gene on chromosome is -, second gene on chromo is +

    # create lists for groups of overlapping genes
    Nested, Piggyback, Convergent, Divergent = [], [], [], []

    # find intronic nested genes first
    ContainedGenes = FindContainedGenePairs(GeneCoord, OverlappingGenes)
    # identify itnronic nested genes {host_gene: [intronic_nested_gene]}
    HostGenes = FindIntronicNestedGenePairs(ContainedGenes, CombinedIntronCoord, GeneCoord)
    # make a list of Host, nested gene pairs
    Nested = GetHostNestedPairs(HostGenes)

    # make a list of set of genes to remove nested relationships
    NestedSets = []
    for pair in Nested:
        NestedSets.append(set(pair))

    # find piggyback genes, convergent and divergent genes
    # loop over overlapping gene pairs
    # first gene in the pair comes first on chromosome
    # but sometimes both members of the pair have same starting position
    for pair in OverlappingPairs:
        # check that pair is not nested
        if set(pair) not in NestedSets:
            # check orientation of both genes
            if GeneCoord[pair[0]][-1] == GeneCoord[pair[1]][-1]:
                # same orientation, piggyback gene pairs
                Piggyback.append(pair)
            elif GeneCoord[pair[0]][-1] != GeneCoord[pair[1]][-1]:
                # check orientation of first and second gene in the pair
                assert GeneCoord[pair[0]][1] <= GeneCoord[pair[1]][1]
                if GeneCoord[pair[0]][-1] == '+':
                    assert GeneCoord[pair[1]][-1] == '-'
                    # convergent gene pair
                    Convergent.append(pair)
                elif GeneCoord[pair[0]][-1] == '-':
                    assert GeneCoord[pair[1]][-1] == '+'
                    # divergent gene pairs
                    Divergent.append(pair)
            
    print(Species[i] + ' nested:{0}, piggyback:{1}, convergent:{2}, divergent:{3}'.format(len(Nested), len(Piggyback), len(Convergent), len(Divergent)))
    assert len(OverlappingPairs) == len(Nested) + len(Piggyback) + len(Convergent) + len(Divergent)
    
    # save overlapping relationships to json files

    # save nested genes as json file
    newfile = open(Species[i] + 'NestedGenes.json', 'w')
    json.dump(HostGenes, newfile, sort_keys = True, indent = 4)
    newfile.close()

    # convert Piggyback gene pairs into a dict and save as a json file
    SaveOverlappingPairsToFile(Piggyback, Species[i] + 'PiggyBackGenes')

    # convert Convergent gene pairs into a dict and save as a json file
    SaveOverlappingPairsToFile(Convergent, Species[i] + 'ConvergentGenes')

    # convert Divergent gene pairs into a dict and save as a json file
    SaveOverlappingPairsToFile(Divergent, Species[i] + 'DivergentGenes')

 