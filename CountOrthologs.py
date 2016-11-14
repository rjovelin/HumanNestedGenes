# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 17:19:12 2016

@author: RJovelin
"""




import os
import sys
import numpy as np
import random
import json
import math
from HsaNestedGenes import *




# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
    
    
        
# get the coordinates of human genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
HumanGeneChromoCoord = ChromoGenesCoord(HsaGFF)
print('got gene coordinates on each chromosome')

# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
HumanGeneCoord = FromChromoCoordToGeneCoord(HumanGeneChromoCoord)
print('got gene coordinates', len(HumanGeneCoord))

# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
HumanOrderedGenes = OrderGenesAlongChromo(HumanGeneChromoCoord)
print('ordered genes on chromsomes')

# Find overlapping genes {gene1: [gene2, gene3]}
HumanOverlappingGenes = FindOverlappingGenePairs(HumanGeneChromoCoord, HumanOrderedGenes)
print('found overlapping genes', len(HumanOverlappingGenes))

# Find genes fully contained in another gene {containing: [contained1, contained2]}
HumanContainedGenes = FindContainedGenePairs(HumanGeneCoord, HumanOverlappingGenes)
print('found genes contained in other genes', len(HumanContainedGenes))

# Map Transcript names to gene names {transcript: gene}
HumanMapTranscriptGene = TranscriptToGene(HsaGFF)
print('mapped transcripts to their parent gene', len(HumanMapTranscriptGene))

# get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
HumanExonCoord = GeneExonCoord(HsaGFF)
print('got exon coordinates', len(HumanExonCoord))
HumanExonCoord = CleanGeneFeatureCoord(HumanExonCoord, HumanMapTranscriptGene)
print('cleaned up exon coordinates of non-mRNA transcripts', len(HumanExonCoord))

# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
HumanIntronCoord = GeneIntronCoord(HumanExonCoord)
print('got intron coordinates', len(HumanIntronCoord))
HumanIntronCoord = CleanGeneFeatureCoord(HumanIntronCoord, HumanMapTranscriptGene)
print('cleaned up intron coordinates of non-mRNA transcripts', len(HumanIntronCoord))

# get the transcript coordinates
HumanTranscriptCoord = TranscriptsCoord(HsaGFF)
print('got transcript coordinates', len(HumanTranscriptCoord))
   
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
HumanCombinedIntronCoord = CombineAllGeneRegions(HumanIntronCoord, HumanMapTranscriptGene)

# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
HumanHostGenes = FindIntronicNestedGenePairs(HumanContainedGenes, HumanCombinedIntronCoord, HumanGeneCoord)
print('identified intronic nested genes', len(HumanHostGenes))

# make a list of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
print('host-gene pairs in human', len(HumanHostNestedPairs))


# get the 1:1 orthologs between human and mouse
HumanMouseOrthos = ParseOrthologFile('HumanMouseOrthologs.txt')
print('mapped human genes to their orthologs in mouse', len(HumanMouseOrthos))
HumanDogOrthos = ParseOrthologFile('HumanDogOrthologs.txt')
HumanChickenOrthos = ParseOrthologFile('HumanChickenOrthologs.txt')
HumanFishOrthos = ParseOrthologFile('HumanZebrafishOrthologs.txt')

orthomouse, orthochicken, orthodog, orthofish = 0, 0, 0, 0
for pair in HumanHostNestedPairs:
    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos:
        orthomouse += 1
    if pair[0] in HumanDogOrthos and pair[1] in HumanDogOrthos:
        orthodog += 1
    if pair[0] in HumanChickenOrthos and pair[1] in HumanChickenOrthos:
        orthochicken += 1
    if pair[0] in HumanFishOrthos and pair[1] in HumanFishOrthos:
        orthofish += 1
print('ortho mouse', orthomouse)
print('ortho dog', orthodog)
print('ortho chicken', orthochicken)
print('ortho fish', orthofish)

