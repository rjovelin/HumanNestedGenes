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




#
#
#
#
## get the coordinates of introns # {chromo: {transcript:[(intron_start, intron_end), ...]}}
#CelIntronCoordChromo = CDSIntronCoord(CelGFF, 'intron')
#print('extracted intronic coordinates')
#print('mapped transcripts to genes')
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
  
