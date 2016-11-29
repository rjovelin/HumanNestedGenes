# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:30:01 2016

@author: RJovelin
"""


# import modules
import os
import sys
from HsaNestedGenes import *


# use this script to 










# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
    
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
HumanGeneChromoCoord = ChromoGenesCoord(HsaGFF)
# map each gene to its mRNA transcripts
HumanMapGeneTranscript = GeneToTranscripts(HsaGFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
HumanGeneChromoCoord = FilterOutGenesWithoutValidTranscript(HumanGeneChromoCoord, HumanMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
HumanGeneCoord = FromChromoCoordToGeneCoord(HumanGeneChromoCoord)
# Map Transcript names to gene names {transcript: gene}
HumanMapTranscriptGene = TranscriptToGene(HsaGFF)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
HumanTranscriptCoordinates = TranscriptsCoord(HsaGFF)
# map genes to their longest transcript {gene: longest_transcript}
HumanGeneLongestTranscript = LongestTranscript(HumanTranscriptCoordinates, HumanMapGeneTranscript)
# get all the CDS sequences in human
HumanCDS = ConvertCDSToFasta('Homo_sapiens.GRCh38.cds.all.fa')
print('human CDS', len(HumanCDS))

# match the sequence of the longest transcript to its parent gene
HumanCodingSeq = {}
for gene in HumanGeneLongestTranscript:
    # get transcript name
    transcript = HumanGeneLongestTranscript[gene]
    # get the corresponding sequence
    sequence = HumanCDS[transcript]
    HumanCodingSeq[gene] = sequence

print(len(HumanCodingSeq), len(HumanGeneCoord))



# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
ChimpGeneChromoCoord = ChromoGenesCoord(PtrGFF)
# map each gene to its mRNA transcripts
ChimpMapGeneTranscript = GeneToTranscripts(PtrGFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
ChimpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(ChimpGeneChromoCoord, ChimpMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
ChimpGeneCoord = FromChromoCoordToGeneCoord(ChimpGeneChromoCoord)
# Map Transcript names to gene names {transcript: gene}
ChimpMapTranscriptGene = TranscriptToGene(PtrGFF)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
ChimpTranscriptCoordinates = TranscriptsCoord(PtrGFF)
# map genes to their longest transcript {gene: longest_transcript}
ChimpGeneLongestTranscript = LongestTranscript(ChimpTranscriptCoordinates, ChimpMapGeneTranscript)
# get all the CDS sequences in human
ChimpCDS = ConvertCDSToFasta('Pan_troglodytes.CHIMP2.1.4.cds.all.fa')
print('chimp CDS', len(ChimpCDS))






