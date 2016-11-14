# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 16:25:38 2016

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
MmuGFF = 'Mus_musculus.GRCm38.86.gff3'    
    
        
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

# Map genes with all their transcripts {gene: [transcripts]}
HumanMapGeneTranscript = GeneToTranscripts(HsaGFF)
print('mapped genes to all their transcripts', len(HumanMapGeneTranscript))

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

# get the CDS coordinates of all transcripts {transcript: [[CDS_start, CDS_end]]}
HumanCDSCoord = GeneCDSCoord(HsaGFF)
print('got CDS coordinates', len(HumanCDSCoord))
HumanCDSCoord = CleanGeneFeatureCoord(HumanCDSCoord, HumanMapTranscriptGene)
print('cleaned up CDS coordinates of non-mRNA transcripts', len(HumanCDSCoord))
 
# get the transcript coordinates
HumanTranscriptCoord = TranscriptsCoord(HsaGFF)
print('got transcript coordinates', len(HumanTranscriptCoord))
   
# get the longest mRNA transcript of each gene
HumanLongestTranscripts = LongestTranscript(HumanTranscriptCoord, HumanMapGeneTranscript)
print('got longest transcripts', len(HumanLongestTranscripts))    



# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
HumanCombinedIntronCoord = CombineAllGeneRegions(HumanIntronCoord, HumanMapTranscriptGene)

# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
HumanHostGenes = FindIntronicNestedGenePairs(HumanContainedGenes, HumanCombinedIntronCoord, HumanGeneCoord)
print('identified intronic nested genes', len(HumanHostGenes))

# make a list of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
print('host-gene pairs in human', len(HumanHostNestedPairs))



same, opposite = 0, 0
for pair in  HumanHostNestedPairs:
    orientation = GenePairOrientation(pair, HumanGeneCoord)
    assert '-' in orientation or '+' in orientation
    if len(set(orientation)) == 1:
        same += 1
    elif len(set(orientation)) == 2:
        opposite += 1
print('same', same)        
print('opposite', opposite)













## Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
#CelIntronicCoord = CombineAllGeneRegions(CelIntronCoordChromo, CelMapTranscriptGene)
#print('combined intron coordinates per gene')
## Find nested genes {host: [intronic nested genes]}
#CelHostGenes = FindIntronicNestedGenePairs(CelContainedGenes, CelIntronicCoord, CelGeneCoord)
#print('identified nested genes', len(CelHostGenes))
## get elegans expression {wormbase ID gene: [expression at 10 stages]}
#CelExpressionStage = ExpressionDevelopemtStages('WBPaper00041190.ce.mr.csv')
#print('expression', len(CelExpressionStage))
## remove host and nested genes without expression
#CelHostGenes = FilterHostNestedExpression(CelHostGenes, CelExpressionStage)
#print('removed host, nested genes without expression', len(CelHostGenes))



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

# Map genes with all their transcripts {gene: [transcripts]}
MouseMapGeneTranscript = GeneToTranscripts(MmuGFF)
print('mapped genes to all their transcripts', len(MouseMapGeneTranscript))

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

# get the CDS coordinates of all transcripts {transcript: [[CDS_start, CDS_end]]}
MouseCDSCoord = GeneCDSCoord(MmuGFF)
print('got CDS coordinates', len(MmuCDSCoord))
MouseCDSCoord = CleanGeneFeatureCoord(MouseCDSCoord, MouseMapTranscriptGene)
print('cleaned up CDS coordinates of non-mRNA transcripts', len(MouseCDSCoord))
 
# get the transcript coordinates
MouseTranscriptCoord = TranscriptsCoord(MmuGFF)
print('got transcript coordinates', len(MouseTranscriptCoord))
   
# get the longest mRNA transcript of each gene
MouseLongestTranscripts = LongestTranscript(MouseTranscriptCoord, MouseMapGeneTranscript)
print('got longest transcripts', len(MouseLongestTranscripts))    



# Combine all intron from all transcripts for a given gene {gene: [(region_start, region_end), ...]}
MouseCombinedIntronCoord = CombineAllGeneRegions(MouseIntronCoord, MouseMapTranscriptGene)

# identify itnronic nested genes {host_gene: [intronic_nested_gene]}
MouseHostGenes = FindIntronicNestedGenePairs(MouseContainedGenes, MouseCombinedIntronCoord, MouseGeneCoord)
print('identified intronic nested genes', len(MouseHostGenes))

# make a list of host-nested gene pairs
MouseHostNestedPairs = GetHostNestedPairs(MouseHostGenes)
print('host-gene pairs in mouse', len(MouseHostNestedPairs))

same, opposite = 0, 0
for pair in  MouseHostNestedPairs:
    orientation = GenePairOrientation(pair, MouseGeneCoord)
    assert '-' in orientation or '+' in orientation
    if len(set(orientation)) == 1:
        same += 1
    elif len(set(orientation)) == 2:
        opposite += 1
print('same', same)        
print('opposite', opposite)



  
