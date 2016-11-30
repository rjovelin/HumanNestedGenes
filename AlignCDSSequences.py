# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:30:01 2016

@author: RJovelin
"""


# import modules
import os
import sys
from HsaNestedGenes import *


# use this script to generate codon-based alignments of human-chimp orthologous coding sequences

# usage python3 AlignCDSSequences.py


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
# remove genes with sequences not multiple of 3, lacking start codons, with internal stop codons
HumanCodingSeq = FilerOutCDSSequences(HumanCodingSeq)
# remove terminal stop codons
HumanCodingSeq = RemoveTerminalStop(HumanCodingSeq)
 


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
# match the sequence of the longest transcript to its parent gene
ChimpCodingSeq = {}
for gene in ChimpGeneLongestTranscript:
    # get transcript name
    transcript = ChimpGeneLongestTranscript[gene]
    # get the corresponding sequence
    sequence = ChimpCDS[transcript]
    ChimpCodingSeq[gene] = sequence
# remove genes with sequences not multiple of 3, lacking start codons, with internal stop codons
ChimpCodingSeq = FilerOutCDSSequences(ChimpCodingSeq)
# remove terminal stop codons
ChimpCodingSeq = RemoveTerminalStop(ChimpCodingSeq)
 

# get the 1:1 orthologs between human-chimp-gorilla
HsaPtrGgoOrthos = ParseOrthologFile('HumanChimpGorillaOrthologs.txt')
# create a dict with human-chimp orthologs
HumanChimpOrthos = {}
for gene in HsaPtrGgoOrthos:
    HumanChimpOrthos[gene] = HsaPtrGgoOrthos[gene][0]


# create a directory named pairs
os.mkdir('./pairs/')


# save orthologous sequences into separate files
for gene in HumanChimpOrthos:
    # check that both human and chimp genes have valid coding sequences
    if gene in HumanCodingSeq and HumanChimpOrthos[gene] in ChimpCodingSeq:
        # save human: chimp orthologous sequences in separate fasta files
        # use human and chimp gene names in filename 
        newfile = open('./pairs/' + gene + '_' + HumanChimpOrthos[gene] + '.tfa', 'w')
        # write human gene and its sequence
        newile.write('>' + gene + '\n')
        newfile.write(HumanCodingSeq[gene] + '\n')
        newfile.write('>' + HumanChimpOrthos[gene] + '\n')
        newfile.write(ChimpCodingSeq[HumanChimpOrthos[gene]] + '\n')
        newfile.close()
        
        
# run to-coffee command lines to align protein sequences and get codon-based DNA alignments

# make a list of fasta files with un-aligned orthologous sequences
files = os.listdir('./pairs')

# fasta files have a .tfa extension
for filename in files:
    # get a outputfile name for the protein translation
    protein = filename[:-3] + '_p.tfa'
    # translate dna sequence into protein sequence
    dnaToPepCommand = "t_coffee -other_pg seq_reformat -in " + './pairs/' + filename + " -action +translate -output fasta_seq > " + protein
    print(dnaToPepcommand)
    os.system(dnaToPepcommand)
    # when sequences have different length, t-coffee padds the short sequence with 'o'
    # replace 'o' with '-'
    os.system("sed -i 's/o/-/g '" + protein)

    # align protein sequences
    alignPepCommand = "t_coffee " + protein
    print(alignPepCommand)
    os.system(alignPepCommand)

    # get name for aligned proteins
    protein_ali = protein[:-3] + '.aln'
    # get name for aligned DNA
    DNA_ali = filename[:-3] +'.aln'
    backtransCommand = 't_coffee -other_pg seq_reformat -in ' + './pairs/' + filename + ' -in2 ' + './pairs/' + protein_ali + ' -action +thread_dna_on_prot_aln' +\
                       ' -output clustalw > ' + DNA_ali
    print(backtransCommand)
    os.system(backtransCommand)

    # convert clustal format to fasta format
    DNA_ali_fasta = filename[:-3] + '_aln.tfa'
    convertDNAalignToFastaCommand = 't_coffee -other_pg seq_reformat -in' + './pairs'/ + DNA_ali + '-output fasta_aln > ' + DNA_ali_fasta
    print(convertDNAalignToFastaCommand)
    os.system(convertDNAalignToFastaCommand)