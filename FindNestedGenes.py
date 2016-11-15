# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 16:25:38 2016

@author: RJovelin
"""

# use this script to identify nested genes in human, chimp, mouse and dog
# save dictionary of genes contained in a host genes as a json file
# save dictionary of intronic nested genes in a host gene as a json file

# usage python3 FindNestedGenes.py


import os
import sys
import numpy as np
import random
import json
import math
from HsaNestedGenes import *


# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
MmuGFF = 'Mus_musculus.GRCm38.86.gff3'    
CfaGFF = 'Canis_familiaris.CanFam3.1.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
Mca = 'Macaca_mulatta.Mmul_8.0.1.86.gff3'    
    
# find contained and intronic genes in human   
    
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
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
HumanCombinedIntronCoord = CombineAllGeneRegions(HumanIntronCoord, HumanMapTranscriptGene)
print('combined introns for each gene', len(HumanCombinedIntronCoord))
# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
HumanHostGenes = FindIntronicNestedGenePairs(HumanContainedGenes, HumanCombinedIntronCoord, HumanGeneCoord)
print('identified intronic nested genes', len(HumanHostGenes))

# get gene expression in human
# filter nested genes based on expression

# save contained genes as json file
newfile = open('HumanContainedGenes.json', 'w')
json.dump(HumanContainedGenes, newfile, sort_keys = True, indent = 4)
newfile.close()
# save intronic nested genes as json file
newfile = open('HumanHostNestedGenes.json', 'w')
json.dump(HumanHostGenes, newfile, sort_keys = True, indent = 4)
newfile.close()


# find contained and intronic nested genes in mouse

# get the coordinates of mouse genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
MouseGeneChromoCoord = ChromoGenesCoord(MmuGFF)
print('got gene coordinates on each chromosome')
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
MouseGeneCoord = FromChromoCoordToGeneCoord(MouseGeneChromoCoord)
print('got gene coordinates', len(MouseGeneCoord))
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
MouseOrderedGenes = OrderGenesAlongChromo(MouseGeneChromoCoord)
print('ordered genes on chromsomes')
# Find overlapping genes {gene1: [gene2, gene3]}
MouseOverlappingGenes = FindOverlappingGenePairs(MouseGeneChromoCoord, MouseOrderedGenes)
print('found overlapping genes', len(MouseOverlappingGenes))
# Find genes fully contained in another gene {containing: [contained1, contained2]}
MouseContainedGenes = FindContainedGenePairs(MouseGeneCoord, MouseOverlappingGenes)
print('found genes contained in other genes', len(MouseContainedGenes))
# Map Transcript names to gene names {transcript: gene}
MouseMapTranscriptGene = TranscriptToGene(MmuGFF)
print('mapped transcripts to their parent gene', len(MouseMapTranscriptGene))
# get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
MouseExonCoord = GeneExonCoord(MmuGFF)
print('got exon coordinates', len(MouseExonCoord))
MouseExonCoord = CleanGeneFeatureCoord(MouseExonCoord, MouseMapTranscriptGene)
print('cleaned up exon coordinates of non-mRNA transcripts', len(MouseExonCoord))
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
MouseIntronCoord = GeneIntronCoord(MouseExonCoord)
print('got intron coordinates', len(MouseIntronCoord))
MouseIntronCoord = CleanGeneFeatureCoord(MouseIntronCoord, MouseMapTranscriptGene)
print('cleaned up intron coordinates of non-mRNA transcripts', len(MouseIntronCoord))
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
MouseCombinedIntronCoord = CombineAllGeneRegions(MouseIntronCoord, MouseMapTranscriptGene)
print('combined introns for each gene', len(MouseCombinedIntronCoord))
# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
MouseHostGenes = FindIntronicNestedGenePairs(MouseContainedGenes, MouseCombinedIntronCoord, MouseGeneCoord)
print('identified intronic nested genes', len(MouseHostGenes))


# get gene expression in mouse
# filter nested genes based on expression

# save contained genes as json file
newfile = open('MouseContainedGenes.json', 'w')
json.dump(MouseContainedGenes, newfile, sort_keys = True, indent = 4)
newfile.close()
# save intronic nested genes as json file
newfile = open('MouseHostNestedGenes.json', 'w')
json.dump(MouseHostGenes, newfile, sort_keys = True, indent = 4)
newfile.close()


# find contained and intronic nested genes in dog

# get the coordinates of human genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
DogGeneChromoCoord = ChromoGenesCoord(CfaGFF)
print('got gene coordinates on each chromosome')
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
DogGeneCoord = FromChromoCoordToGeneCoord(DogGeneChromoCoord)
print('got gene coordinates', len(DogGeneCoord))
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
DogOrderedGenes = OrderGenesAlongChromo(DogGeneChromoCoord)
print('ordered genes on chromsomes')
# Find overlapping genes {gene1: [gene2, gene3]}
DogOverlappingGenes = FindOverlappingGenePairs(DogGeneChromoCoord, DogOrderedGenes)
print('found overlapping genes', len(DogOverlappingGenes))
# Find genes fully contained in another gene {containing: [contained1, contained2]}
DogContainedGenes = FindContainedGenePairs(DogGeneCoord, DogOverlappingGenes)
print('found genes contained in other genes', len(DogContainedGenes))
# Map Transcript names to gene names {transcript: gene}
DogMapTranscriptGene = TranscriptToGene(CfaGFF)
print('mapped transcripts to their parent gene', len(DogMapTranscriptGene))
# get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
DogExonCoord = GeneExonCoord(CfaGFF)
print('got exon coordinates', len(DogExonCoord))
DogExonCoord = CleanGeneFeatureCoord(DogExonCoord, DogMapTranscriptGene)
print('cleaned up exon coordinates of non-mRNA transcripts', len(DogExonCoord))
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
DogIntronCoord = GeneIntronCoord(DogExonCoord)
print('got intron coordinates', len(DogIntronCoord))
DogIntronCoord = CleanGeneFeatureCoord(DogIntronCoord, DogMapTranscriptGene)
print('cleaned up intron coordinates of non-mRNA transcripts', len(DogIntronCoord))
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
DogCombinedIntronCoord = CombineAllGeneRegions(DogIntronCoord, DogMapTranscriptGene)
print('combined introns for each gene', len(DogCombinedIntronCoord))
# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
DogHostGenes = FindIntronicNestedGenePairs(DogContainedGenes, DogCombinedIntronCoord, DogGeneCoord)
print('identified intronic nested genes', len(DogHostGenes))

# get gene expression in human
# filter nested genes based on expression

# save contained genes as json file
newfile = open('DogContainedGenes.json', 'w')
json.dump(DogContainedGenes, newfile, sort_keys = True, indent = 4)
newfile.close()
# save intronic nested genes as json file
newfile = open('DogHostNestedGenes.json', 'w')
json.dump(DogHostGenes, newfile, sort_keys = True, indent = 4)
newfile.close()


# find contained and intronic genes in chimp   
    
# get the coordinates of human genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
ChimpGeneChromoCoord = ChromoGenesCoord(PtrGFF)
print('got gene coordinates on each chromosome')
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
ChimpGeneCoord = FromChromoCoordToGeneCoord(ChimpGeneChromoCoord)
print('got gene coordinates', len(ChimpGeneCoord))
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
ChimpOrderedGenes = OrderGenesAlongChromo(ChimpGeneChromoCoord)
print('ordered genes on chromsomes')
# Find overlapping genes {gene1: [gene2, gene3]}
ChimpOverlappingGenes = FindOverlappingGenePairs(ChimpGeneChromoCoord, ChimpOrderedGenes)
print('found overlapping genes', len(ChimpOverlappingGenes))
# Find genes fully contained in another gene {containing: [contained1, contained2]}
ChimpContainedGenes = FindContainedGenePairs(ChimpGeneCoord, ChimpOverlappingGenes)
print('found genes contained in other genes', len(ChimpContainedGenes))
# Map Transcript names to gene names {transcript: gene}
ChimpMapTranscriptGene = TranscriptToGene(PtrGFF)
print('mapped transcripts to their parent gene', len(ChimpMapTranscriptGene))
# get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
ChimpExonCoord = GeneExonCoord(PtrGFF)
print('got exon coordinates', len(ChimpExonCoord))
ChimpExonCoord = CleanGeneFeatureCoord(ChimpExonCoord, ChimpMapTranscriptGene)
print('cleaned up exon coordinates of non-mRNA transcripts', len(ChimpExonCoord))
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
ChimpIntronCoord = GeneIntronCoord(ChimpExonCoord)
print('got intron coordinates', len(ChimpIntronCoord))
ChimpIntronCoord = CleanGeneFeatureCoord(ChimpIntronCoord, ChimpMapTranscriptGene)
print('cleaned up intron coordinates of non-mRNA transcripts', len(ChimpIntronCoord))
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
ChimpCombinedIntronCoord = CombineAllGeneRegions(ChimpIntronCoord, ChimpMapTranscriptGene)
print('combined introns for each gene', len(ChimpCombinedIntronCoord))
# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
ChimpHostGenes = FindIntronicNestedGenePairs(ChimpContainedGenes, ChimpCombinedIntronCoord, ChimpGeneCoord)
print('identified intronic nested genes', len(ChimpHostGenes))

# get gene expression in human
# filter nested genes based on expression

# save contained genes as json file
newfile = open('ChimpContainedGenes.json', 'w')
json.dump(ChimpContainedGenes, newfile, sort_keys = True, indent = 4)
newfile.close()
# save intronic nested genes as json file
newfile = open('ChimpHostNestedGenes.json', 'w')
json.dump(ChimpHostGenes, newfile, sort_keys = True, indent = 4)
newfile.close()



# find contained and intronic genes in macaque 
    
# get the coordinates of human genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
MacaqueGeneChromoCoord = ChromoGenesCoord(McaGFF)
print('got gene coordinates on each chromosome')
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
MacaqueGeneCoord = FromChromoCoordToGeneCoord(MacaqueGeneChromoCoord)
print('got gene coordinates', len(MacaqueGeneCoord))
# Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
MacaqueOrderedGenes = OrderGenesAlongChromo(MacaqueGeneChromoCoord)
print('ordered genes on chromsomes')
# Find overlapping genes {gene1: [gene2, gene3]}
MacaqueOverlappingGenes = FindOverlappingGenePairs(MacaqueGeneChromoCoord, MacaqueOrderedGenes)
print('found overlapping genes', len(MacaqueOverlappingGenes))
# Find genes fully contained in another gene {containing: [contained1, contained2]}
MacaqueContainedGenes = FindContainedGenePairs(MacaqueGeneCoord, MacaqueOverlappingGenes)
print('found genes contained in other genes', len(MacaqueContainedGenes))
# Map Transcript names to gene names {transcript: gene}
MacaqueMapTranscriptGene = TranscriptToGene(McaGFF)
print('mapped transcripts to their parent gene', len(MacaqueMapTranscriptGene))
# get the exon coordinates of all transcript {transcript: [[exon_start, exon_end]]}
MacaqueExonCoord = GeneExonCoord(McaGFF)
print('got exon coordinates', len(MacaqueExonCoord))
MacaqueExonCoord = CleanGeneFeatureCoord(MacaqueExonCoord, MacaqueMapTranscriptGene)
print('cleaned up exon coordinates of non-mRNA transcripts', len(MacaqueExonCoord))
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
MacaqueIntronCoord = GeneIntronCoord(MacaqueExonCoord)
print('got intron coordinates', len(MacaqueIntronCoord))
MacaqueIntronCoord = CleanGeneFeatureCoord(MacaqueIntronCoord, MacaqueMapTranscriptGene)
print('cleaned up intron coordinates of non-mRNA transcripts', len(MacaqueIntronCoord))
# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
MacaqueCombinedIntronCoord = CombineAllGeneRegions(MacaqueIntronCoord, MacaqueMapTranscriptGene)
print('combined introns for each gene', len(MacaqueCombinedIntronCoord))
# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
MacaqueHostGenes = FindIntronicNestedGenePairs(MacaqueContainedGenes, MacaqueCombinedIntronCoord, MacaqueGeneCoord)
print('identified intronic nested genes', len(MacaqueHostGenes))

# get gene expression in human
# filter nested genes based on expression

# save contained genes as json file
newfile = open('MacaqueContainedGenes.json', 'w')
json.dump(MacaqueContainedGenes, newfile, sort_keys = True, indent = 4)
newfile.close()
# save intronic nested genes as json file
newfile = open('MacaqueHostNestedGenes.json', 'w')
json.dump(MacaqueHostGenes, newfile, sort_keys = True, indent = 4)
newfile.close()