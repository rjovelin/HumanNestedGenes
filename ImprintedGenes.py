# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:29:11 2017

@author: RJovelin
"""

# use this script to test enrichement of imprinted genes among nested genes



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




# get option from command
SisterSp = sys.argv[1]
Analysis = sys.argv[2]
assert SisterSp in ['mouse', 'chimp']
assert Analysis in ['pairs', 'orthos']

# load dictionaries of overlapping genes
if SisterSp == 'chimp':
    # sister species is chimp and outgroup is gorilla
    jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
                 'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json',
                 'GorillaOverlappingGenes.json', 'GorillaNestedGenes.json']
    # get GFF file
    GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3', 'Gorilla_gorilla.gorGor3.1.86.gff3']
elif SisterSp == 'mouse':
    # sister species is mouse and outgroup is dog
    jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
                 'MouseOverlappingGenes.json', 'MouseNestedGenes.json',
                 'DogOverlappingGenes.json', 'DogNestedGenes.json']
    # get GFF file
    GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Mus_musculus.GRCm38.86.gff3', 'Canis_familiaris.CanFam3.1.87.gff3']
    

# make a list of dictionaries
AllOverlap = []
# loop over files
for i in range(len(jsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(jsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    AllOverlap.append(overlapping)

# make a list of gene coordinates       
AllCoordinates, AllOrdered = [], []
# loop over GFF files
for i in range(len(GFF)):
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    AllCoordinates.append(GeneCoord)
    AllOrdered.append(OrderedGenes)

# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# make pairs of overlapping genes
HumanPairs = AllPairs[:2]
SisterPairs = AllPairs[2:4]
OutGroupPairs = AllPairs[4:]

# make list with sets of non-overlapping genes
NonOverlappingSets = []
for i in range(3):
    j = i * 2
    # make a set of non-overlapping gene
    nonoverlap = MakeNonOverlappingGeneSet(AllOverlap[j], AllCoordinates[i])
    NonOverlappingSets.append(nonoverlap)    

# make sets of host and nested nested genes
NestedSets = []
for i in range(1, len(AllOverlap), 2):
    nestedset = MakeFullPartialOverlapGeneSet(AllOverlap[i])
    NestedSets.append(nestedset)

# make sets of overlapping genes
OverlapSets = []
for i in range(0, len(AllOverlap), 2):
    overlap = MakeFullPartialOverlapGeneSet(AllOverlap[i])
    OverlapSets.append(overlap)

# get 1:1 orthologs between human and sister-species
# get 1:1 orthologs between human, sister-species and outgroup {human:[sistersp,outgroup]}
if SisterSp == 'chimp':
    OrthoPairs = MatchOrthologPairs('HumanChimpOrthologs.txt')
    OrthoTrios = MatchOrthologTrios('HumanChimpGorillaOrthologs.txt')
elif SisterSp == 'mouse':
    OrthoPairs = MatchOrthologPairs('HumanMouseOrthologs.txt')
    OrthoTrios = MatchOrthologTrios('HumanMouseDogOrthologs.txt')

# reverse dict with human and sister-species orthologs
SisterOrthos = {}
for gene in OrthoPairs:
    SisterOrthos[OrthoPairs[gene]] = gene
# make a dict of ortho trios with sister-species genes as key
SisterOrthoTrios = {}
for gene in OrthoTrios:
    sistergene, outgroupgene = OrthoTrios[gene][0], OrthoTrios[gene][1]
    SisterOrthoTrios[sistergene] = [gene, outgroupgene]






# use this function to parse the files with imprinted gene information
def ParseImprinted(ImprintedGenesFile):
    '''
    (file) -> dict
    Take the file with information on imprinted genes and return a dictionary 
    of imprinted genes (gene name) and gene status (paternal/maternal) pairs    
    '''
    
    # make a dictionary with gene name and imprinting status
    imprinted = {}
    infile = open(ImprintedGenesFile)
    infile.readline()
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            # get gene status
            if line[3] in ['Imprinted', 'Predicted']:
                # get expressed allele
                allele = line[-1]
                # check for typo in the file
                if allele == 'Paterna':
                    allele = 'Paternal'
                # add gene name: allele pair to dict
                if line[0] in imprinted:
                    assert imprinted[line[0]] == allele
                else:
                    imprinted[line[0]] = allele
                # check if gene has aliases
                line[1] = line[1].replace(' ', '')
                if line[1] != '':
                    if ',' not in line[1]:
                        if line[1] in imprinted:
                            assert imprinted[line[1]] == allele
                        else:
                            imprinted[line[1]] = allele
                    else:
                        # get all alias names
                        line[1] = line[1].split(',')
                        for i in line[1]:
                            if i in imprinted:
                                assert imprinted[i] == allele
                            else:
                                imprinted[i] = allele
    infile.close()
    return imprinted
    

# use this function to map gene names with gene ID
def MapNametoID(GFF):
    '''
    (file) -> dict    
    Take the GFF annotation file and return a dictionary with gene ID: gene name pairs
    '''
    # create a dict {gene_ID: Name}    
    Names = {}
    infile = open(GFF, 'r')
    for line in infile:
        line = line.rstrip()
        if 'gene' in line and 'protein_coding' in line:
            line = line.rstrip().split('\t')
            # check that gene is protein-coding            
            assert line[2] == 'gene'
            biotype = line[8][line[8].index('biotype=')+8: line[8].index(';', line[8].index('biotype=')+8)]
            assert biotype == 'protein_coding'
            # get the gene ID
            gene = line[8][line[8].index('ID=gene:')+8: line[8].index(';')]
            # get gene name
            name = line[8][line[8].index('Name=') + len('Name='): line[8].index(';', line[8].index(';')+1)]
            assert gene not in Names and name not in Names.values()            
            Names[gene] = name
    infile.close()
    return Names    