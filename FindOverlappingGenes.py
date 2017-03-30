# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 16:25:38 2016

@author: RJovelin
"""

# use this script to identify overlapping genes in mammals
# save dictionaries of overlapping genes as json files

# usage python3 FindOverlappingGenes.py 

import os
import sys
import numpy as np
import random
import json
import math
import sys
from HsaNestedGenes import *


# make a list of GFF files
GFF = ['Bos_taurus.UMD3.1.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
       'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
       'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
       'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
       'Gorilla_gorilla.gorGor3.1.88.gff3', 'Homo_sapiens.GRCh38.88.gff3',
       'Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Monodelphis_domestica.BROADO5.88.gff3',
       'Mus_musculus.GRCm38.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3', 
       'Pan_troglodytes.CHIMP2.1.4.88.gff3', 'Pongo_abelii.PPYG2.88.gff3', 
       'Sorex_araneus.COMMON_SHREW1.88.gff3']

# make a parallel list of species names
Species = ['Cow', 'Marmoset', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
           'Cat', 'Gorilla', 'Human', 'Macaque', 'Opossum', 'Mouse', 'Platypus',
           'Chimp', 'Orangutan', 'Shrew']
           
# loop over GFF files
for i in range(len(GFF)):
    print(GFFs[i][:GFFs[i].index('.')], Species[i])
    # get the coordinates of human genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    SpGeneChromoCoord = ChromoGenesCoord(GFFs[i])
    # map each gene to its mRNA transcripts
    SpMapGeneTranscript = GeneToTranscripts(GFFs[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    SpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(SpGeneChromoCoord, SpMapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    SpGeneCoord = FromChromoCoordToGeneCoord(SpGeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    SpOrderedGenes = OrderGenesAlongChromo(SpGeneChromoCoord)
    # Find overlapping genes {gene1: [gene2, gene3]}
    SpOverlappingGenes = FindOverlappingGenePairs(SpGeneChromoCoord, SpOrderedGenes)
    
    # save dictionaries as json file
    print(GFFs[i][:GFFs[i].index('.')], Species[i], len(SpOverlappingGenes))    
    newfile = open(Species + 'OverlappingGenes.json', 'w')  
    # save dictionary to file        
    json.dump(SpOverlappingGenes, newfile, sort_keys = True, indent = 4)
    newfile.close()