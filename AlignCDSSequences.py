# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 15:30:01 2016

@author: RJovelin
"""


# import modules
import os
import sys
from HsaNestedGenes import *


# use this script to generate codon-based alignments of human-chimp or human-mouse orthologous coding sequences

# usage python3 AlignCDSSequences.py [options]
# -[mouse/chimp]: align human seqs with orthologs in mouse or chimp

SisterSpecies = sys.argv[1]
assert SisterSpecies in ['mouse', 'chimp']

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.88.gff3'
if SisterSpecies == 'chimp':
    SpeciesGFF = 'Pan_troglodytes.CHIMP2.1.4.88.gff3'
elif SisterSpecies == 'mouse':
    SpeciesGFF = 'Mus_musculus.GRCm38.88.gff3'


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
SpeciesGeneChromoCoord = ChromoGenesCoord(SpeciesGFF)
# map each gene to its mRNA transcripts
SpeciesMapGeneTranscript = GeneToTranscripts(SpeciesGFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
SpeciesGeneChromoCoord = FilterOutGenesWithoutValidTranscript(SpeciesGeneChromoCoord, SpeciesMapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
SpeciesGeneCoord = FromChromoCoordToGeneCoord(SpeciesGeneChromoCoord)
# Map Transcript names to gene names {transcript: gene}
SpeciesMapTranscriptGene = TranscriptToGene(SpeciesGFF)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
SpeciesTranscriptCoordinates = TranscriptsCoord(SpeciesGFF)
# map genes to their longest transcript {gene: longest_transcript}
SpeciesGeneLongestTranscript = LongestTranscript(SpeciesTranscriptCoordinates, SpeciesMapGeneTranscript)
# get all the CDS sequences in sister species
if SisterSpecies == 'chimp':
    SpeciesCDS = ConvertCDSToFasta('Pan_troglodytes.CHIMP2.1.4.cds.all.fa')
elif SisterSpecies == 'mouse':
    SpeciesCDS = ConvertCDSToFasta('Mus_musculus.GRCm38.cds.all.fa')
print(SisterSpecies + ' CDS', len(SpeciesCDS))
# match the sequence of the longest transcript to its parent gene
SpeciesCodingSeq = {}
for gene in SpeciesGeneLongestTranscript:
    # get transcript name
    transcript = SpeciesGeneLongestTranscript[gene]
    # get the corresponding sequence
    sequence = SpeciesCDS[transcript]
    SpeciesCodingSeq[gene] = sequence
# remove genes with sequences not multiple of 3, lacking start codons, with internal stop codons
SpeciesCodingSeq = FilerOutCDSSequences(SpeciesCodingSeq)
# remove terminal stop codons
SpeciesCodingSeq = RemoveTerminalStop(SpeciesCodingSeq)
 

############ continue here









# get the 1:1 orthologs and destination folder
if SisterSpecies == 'chimp':
    Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
    Folder = './HumanChimpPairs/'
elif SisterSpecies == 'mouse':
    Orthos = MatchOrthologPairs('HumanMouseOrthologs.txt')
    Folder = './HumanMousePairs/'

# create a directory named pairs
os.mkdir(Folder)

# save orthologous sequences into separate files
for gene in Orthos:
    # check that both human and ortholog genes have valid coding sequences
    if gene in HumanCodingSeq and Orthos[gene] in SpeciesCodingSeq:
        # save human: ortholog sequences in separate fasta files
        # use human and ortholog gene names in filename 
        newfile = open(Folder + gene + '_' + Orthos[gene] + '.tfa', 'w')
        # write human gene and its sequence
        newfile.write('>' + gene + '\n')
        newfile.write(HumanCodingSeq[gene] + '\n')
        newfile.write('>' + Orthos[gene] + '\n')
        newfile.write(SpeciesCodingSeq[Orthos[gene]] + '\n')
        newfile.close()
        
        
# run to-coffee command lines to align protein sequences and get codon-based DNA alignments

# make a list of fasta files with un-aligned orthologous sequences
files = os.listdir(Folder)

# fasta files have a .tfa extension
for filename in files:
    # get a outputfile name for the protein translation
    protein = filename[:-4] + '_p.tfa'
    # translate dna sequence into protein sequence
    #t_coffee -other_pg seq_reformat -in ORTHOMCL7896.tfa -action +translate -output fasta_seq > ORTHOMCL7896_p.tfa
    dnaToPepCommand = "t_coffee -other_pg seq_reformat -in " + Folder + filename + " -action +translate -output fasta_seq > " + protein
    print(dnaToPepCommand)
    os.system(dnaToPepCommand)
    # when sequences have different length, t-coffee padds the short sequence with 'o'
    # replace 'o' with '-'
    os.system("sed -i 's/o/-/g '" + protein)

    # align protein sequences
    #t_coffee ORTHOMCL7896_p.tfa
    alignPepCommand = "t_coffee " + protein
    print(alignPepCommand)
    os.system(alignPepCommand)

    # back-translate protein alignment to DNA in MSA format
    #t_coffee -other_pg seq_reformat -in ORTHOMCL7896.tfa -in2 ORTHOMCL7896_p.aln -action +thread_dna_on_prot_aln -output clustalw >ORTHOMCL7896.aln
    # get name for aligned proteins
    protein_ali = protein[:-4] + '.aln'
    # get name for aligned DNA
    DNA_ali = filename[:-4] +'.aln'
    backtransCommand = 't_coffee -other_pg seq_reformat -in ' + Folder + filename + ' -in2 ' + protein_ali + ' -action +thread_dna_on_prot_aln' +\
                       ' -output clustalw > ' + DNA_ali
    print(backtransCommand)
    os.system(backtransCommand)

    # convert clustal format to fasta format
    #t_coffee -other_pg seq_reformat -in ORTHOMCL7896.aln -output fasta_aln > ORTHOMCL7896_aln.tfa
    DNA_ali_fasta = filename[:-4] + '_aln.tfa'
    convertDNAalignToFastaCommand = 't_coffee -other_pg seq_reformat -in ' + DNA_ali + ' -output fasta_aln > ' + DNA_ali_fasta
    print(convertDNAalignToFastaCommand)
    os.system(convertDNAalignToFastaCommand)


##translate multi-fasta DNA to AA
#t_coffee -other_pg seq_reformat -in ORTHOMCL7896.tfa -action +translate -output fasta_seq > ORTHOMCL7896_p.tfa
##align AA sequences
#t_coffee ORTHOMCL7896_p.tfa
##convert AA alignment to fasta format
#t_coffee -other_pg seq_reformat -in ORTHOMCL7896_p.aln -output fasta_aln > ORTHOMCL7896_p-aln.tfa
##back-translate AA alignment to DNA in MSA format
#t_coffee -other_pg seq_reformat -in ORTHOMCL7896.tfa -in2 ORTHOMCL7896_p.aln -action +thread_dna_on_prot_aln -output clustalw >ORTHOMCL7896.aln
##convert DNA alignment to fasta format
#t_coffee -other_pg seq_reformat -in ORTHOMCL7896.aln -output fasta_aln > ORTHOMCL7896_aln.tfa

print('done aligning sequences')

