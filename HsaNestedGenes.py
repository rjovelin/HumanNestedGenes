# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 21:04:34 2016

@author: Richard
"""

import numpy as np
import math
import os
import random
import json

# use this function to record the gene coordinates on each separate chromosome    
def ChromoGenesCoord(gff_file):
    '''
    (file) -> dict
    Return a dictionnary with chromosome as key and a dictionnary of protein-coding
    gene coordinates as value
    '''
      
    # open file to read
    infile = open(gff_file, 'r')
    # make a dictionnary with chromosome as key and a dictionnary as value
    # {chromo: {gene:[chromosome, start, end, sense]}}
    Linkage = {}
    for line in infile:
        line = line.rstrip()
        if 'gene' in line and not line.startswith('#'):
            line = line.split('\t')
            if line[2] == 'gene':
                # get biotype
                biotype = line[8][line[8].index('biotype=')+8: line[8].index(';', line[8].index('biotype=')+8)]
                # record only protein coding genes
                if biotype == 'protein_coding':
                    # get the gene name
                    gene = line[8][line[8].index('ID=gene:')+8: line[8].index(';')]
                    # get chromo, start, end positions 0-based, and orientation
                    chromo, start, end, sense = line[0], int(line[3]) -1, int(line[4]), line[6]            
                    # check if chromo is recorded
                    if chromo not in Linkage:
                        # create a dictionnary
                        Linkage[chromo] = {}
                    # populate inner dict with gene : coords
                    Linkage[chromo][gene] = [chromo, start, end, sense]
    # close file after reading
    infile.close()
    return Linkage


# use this function to remove genes that lack a mRNA transcript
def FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript):
    '''
    (dict, dict) -> dict    
    Take the dictionary with gene coordinates on each chromosome, the dictionary
    of gene: transcript pairs and return a modified dictionary of gene coordinates
    of each chromo for genes that have mRNA transcripts (ie. remove gene that only have
    abberant transcripts, NMD processed transcripts, etc)
    '''
    
    # GeneChromoCoord is in the form # {chromo: {gene:[chromosome, start, end, sense]}}
    # MapGeneTranscript is in the form {gene : [transcripts]}
    
    # remove genes lacking mRNA transcripts before ordering genes on chromos
    for chromo in GeneChromoCoord:
        to_remove = [gene for gene in GeneChromoCoord[chromo] if gene not in MapGeneTranscript]
        if len(to_remove) != 0:
            for gene in to_remove:
                del GeneChromoCoord[chromo][gene]
    return GeneChromoCoord
    

# use this function to get the coordinates of all protein coding genes
def FromChromoCoordToGeneCoord(GeneCoord):
    '''
    (dict) -> dict
    Take the dictionary with gene coordinates on each chromo and return a dictionary
    with gene and coordinates pairs    
    '''
    # GeneCoord is in the form {chromo: {gene:[chromosome, start, end, sense]}}
    # create a dictionary with gene as key {gene:[chromosome, start, end, sense]}
    Genes = {}
    for chromo in GeneCoord:
        for gene in GeneCoord[chromo]:
            Genes[gene] = list(GeneCoord[chromo][gene])
    return Genes


# use this function to record the coordinates of each transcripts
def TranscriptsCoord(gff_file):
    '''
    (file) -> dict
    Return a dictionary with the transcript name and the its coordinates
    '''
    
    # open file to read
    infile = open(gff_file, 'r')
    # make a dictionnary with chromosome as key and a dictionnary as value
    # {transcript:[chromosome, start, end, sense]}
    Transcripts = {}
    for line in infile:
        line = line.rstrip()
        if 'mRNA' in line:
            line = line.split('\t')
            if line[2] == 'mRNA':
                # get the gene name
                transcript = line[8][line[8].index('transcript:')+11 : line[8].index(';')]
                # get chromo, start, end positions 0-based, and orientation
                chromo, start, end, sense = line[0], int(line[3]) -1, int(line[4]), line[6]            
                # populate the dictionnary 
                Transcripts[transcript] = [chromo, start, end, sense]
    # close file after reading
    infile.close()
    return Transcripts
   
   
   
# use this function to create a dict of transcript : gene pairs
def TranscriptToGene(gff_file):
    '''
    (file) -> dict
    Returns a dictionnary with transcript : gene pairs from the gff annotation file
    '''
 
    # create a set of protein-coding genes
    ProteinCoding = set()
    # open file for reading
    infile = open(gff_file, 'r')
    for line in infile:
        line = line.rstrip()
        if 'gene' in line and not line.startswith('#'):
            line = line.split('\t')
            if line[2] == 'gene':
                # get biotype
                biotype = line[8][line[8].index('biotype=')+8: line[8].index(';', line[8].index('biotype=')+8)]
                # record only protein coding genes
                if biotype == 'protein_coding':
                    # get the gene name
                    gene = line[8][line[8].index('ID=gene:')+8: line[8].index(';')]
                    ProteinCoding.add(gene)
    infile.close()
    
    #create a dictionnary of transcript : gene pairs
    transcripts_genes = {}
    # open file for reading
    infile = open(gff_file, 'r')
    for line in infile:
        if 'mRNA' in line:
            line = line.rstrip().split('\t')
            if line[2] == 'mRNA':
                transcript = line[8][line[8].index('transcript:')+11 : line[8].index(';')]
                gene = line[8][line[8].index('Parent=gene:')+12 : line[8].index(';', line[8].index('Parent'))]
                # keep mRNAs of protein coding genes 
                if gene in ProteinCoding:
                    transcripts_genes[transcript] = gene
    infile.close()
    return transcripts_genes


# use this function to create a dict of gene : list of transcripts pairs
def GeneToTranscripts(gff_file):
    '''
    (file) -> dict
    Returns a dictionnary with gene as key and a list of transcripts as value
    '''

    # get the dictionnary of transcripts : gene names pairs
    transcripts_genes = TranscriptToGene(gff_file)
    # create a reverse dictionnary of {gene : [transcripts]} pairs
    genes = {}
    for transcript in transcripts_genes:
        gene_name = transcripts_genes[transcript]
        if gene_name in genes:
            genes[gene_name].append(transcript)
        else:
            genes[gene_name] = [transcript]
    return genes



# use this function to map genes with their longest mRNA transcripts
def LongestTranscript(TranscriptCoordinates, MapGeneTranscript):
    '''
    (dict, dict) -> dict
    Take a dictionary with transcrip coordinates, a dictionary matching all
    transcripts to each gene and return a dictionary of gene: longest transcript pairs
    '''
    # TranscriptCoordinates is in the form {transcript:[chromosome, start, end, sense]}
    # MapGeneTranscript is in the orm {gene: [transcripts]}
    # create a dict {gene: longest_transcript}
    Genes = {}
    for gene in MapGeneTranscript:
        # create list to store transcript and their length
        L = []
        for transcript in MapGeneTranscript[gene]:
            L.append([TranscriptCoordinates[transcript][2] - TranscriptCoordinates[transcript][1], transcript])
        # sort list, get the transcript with maximum length
        L.sort()
        longest = L[-1][1]
        Genes[gene] = longest
    return Genes


# use this function to get the CDS coordinates of all transcripts 
def GeneCDSCoord(gff_file):
    '''
    (file, str) -> dict
    Return a dictionnary with transcript : CDS positions list pairs    
    '''
    
    # open file for reading
    infile = open(gff_file, 'r')
    # make a dictionnary with chromosome as key and a dictionnary as value
    # {transcript:[(region_start, region_end), (region_start, region_end)]}
    CDSPositions = {}
    for line in infile:
        line = line.rstrip()
        if 'CDS' in line:
            line = line.split('\t')
            if line[2] == 'CDS':
                # get start and end positions 0-based
                start, end = int(line[3]) -1, int(line[4])
                # check that exon correspond to one transcript (in some GFF format,
                # exons for more than 1 transcript are collapsed in a single line)
                assert line[8].count('transcript') == 1, 'exon matches multiple transcripts'
                # get the parent transcript
                transcript = line[8][line[8].index('Parent=transcript:')+len('Parent=transcript:'):line[8].index(';protein')]
                # populate dict                
                if transcript not in CDSPositions:
                    CDSPositions[transcript] = [[start, end]]
                else:
                    CDSPositions[transcript].append([start, end])
    # sort exon coordinates
    for transcript in CDSPositions:
        CDSPositions[transcript].sort()
        
    for transcript in CDSPositions:
        for i in range(0, len(CDSPositions[transcript])-1):
            for j in range(i+1, len(CDSPositions[transcript])):
                assert CDSPositions[transcript][j][0] > CDSPositions[transcript][i][1]
    # close file after reading
    infile.close()
    return CDSPositions


# use this function to get the exon coordinates of all transcripts 
def GeneExonCoord(gff_file):
    '''
    (file) -> dict
    Return a dictionnary with transcript as key and exon positions list pairs
    as value. Note that exons define intron positions but are not necessarily
    entirely coding    
    '''
    
    # open file for reading
    infile = open(gff_file, 'r')
    # make a dictionnary of exon coorddinates for each transcript
    # {transcript:[(region_start, region_end), (region_start, region_end)]}
    ExonPositions = {}
    for line in infile:
        line = line.rstrip()
        if 'exon' in line:
            line = line.split('\t')
            if line[2] == 'exon':
                # get start and end positions 0-based
                start, end = int(line[3]) -1, int(line[4])
                # check that exon correspond to one transcript (in some GFF format,
                # exons for more than 1 transcript are collapsed in a single line)
                assert line[8].count('transcript') == 1, 'exon matches multiple transcripts'
                # get the parent transcript
                transcript = line[8][line[8].index('Parent=transcript:')+len('Parent=transcript:'):line[8].index(';')]
                # populate dict                
                if transcript not in ExonPositions:
                    ExonPositions[transcript] = [[start, end]]
                else:
                    ExonPositions[transcript].append([start, end])
    # sort exon coordinates
    for transcript in ExonPositions:
        ExonPositions[transcript].sort()
        
    for transcript in ExonPositions:
        for i in range(0, len(ExonPositions[transcript])-1):
            for j in range(i+1, len(ExonPositions[transcript])):
                assert ExonPositions[transcript][j][0] > ExonPositions[transcript][i][1]
    # close file after reading
    infile.close()
    return ExonPositions


# use this function to get the positions of introns for each transcript 
def GeneIntronCoord(ExonCoord):
    '''
    (dict) -> dict
    Take a dictionary with exon positions for each transcript and return a dictionary
    with intron positions for each transcript
    Precondition: exons which are not necessarily coding are delimiting introns, 
    and are already sorted according to their chromosomal positions
    '''
    # create dictionnary to store positions {transcript: [(s1, end1), (s2, end2)]}
    IntronCoord = {}
    
    # loop over transcript
    for transcript in ExonCoord:
        # check that transcripts has multiple exons (and thus at least 1 intron)
        if len(ExonCoord[transcript]) >= 2:
            # loop over the exon coordinates (already sorted)
            for i in range(len(ExonCoord[transcript]) -1):
                # intron start is end of first exon and intron end is start of next exon
                intronstart = ExonCoord[transcript][i][-1]
                intronend = ExonCoord[transcript][i+1][0]
                assert intronend > intronstart, 'end position should be geater than start position'
                introncoord = [intronstart, intronend]
                # populate dict
                if transcript in IntronCoord:
                    IntronCoord[transcript].append(introncoord)
                else:
                    IntronCoord[transcript] = [introncoord]
    return IntronCoord   


# use this function to remove non-valid transcripts from gene feature coordinates
def CleanGeneFeatureCoord(GeneFeatureCoord, TranscriptToGene):
    '''
    (dict, dict) -> dict
    Take the dictionary of transcript: gene feature coordinates (intron, exon or CDS) 
    and the dictionary with transcript: gene pairs and return a modified dictionary
    of feature coordinates with non-mRNA transcripts removed
    Precondition: Transcript: gene pairs correspond to mRNA transcripts and their
    parent gene. Gene features includes coordinates of abberant transcript, 
    non-coding RNAs, etc     
    '''
    
    # create a list of transcripts to remove
    to_remove = [i for i in GeneFeatureCoord if i not in TranscriptToGene]
    # remove any non-mRNA transcript not paired with their parent gene
    for i in to_remove:
        del GeneFeatureCoord[i]
    return GeneFeatureCoord


# use this function to order genes along chromosomes
def OrderGenesAlongChromo(GeneCoord):
    '''
    (dict) -> dict
    Take a dictionary of gene coordinates per chromosome and return a dictionary
    with chromo as key and a list of ordered gene names    
    '''
    # GeneCoord is in the form {'I': {gene:[chromosome, start, end, sense]}}
    # make a dictionary {chromo: [[start, gene1], [start, gene2]}
    StartPos = {}
    for chromo in GeneCoord:
        # initialize list
        StartPos[chromo] = []
        for gene in GeneCoord[chromo]:
            StartPos[chromo].append([GeneCoord[chromo][gene][1], gene])
    # sort list on start positions
    for chromo in StartPos:
        StartPos[chromo].sort()
    # create dict of ordered gene names {chromo :[ordred gene names]}     
    OrderedGenes = {}
    for chromo in StartPos:
        OrderedGenes[chromo] = []
        for i in StartPos[chromo]:
            OrderedGenes[chromo].append(i[1])
    return OrderedGenes
  

# use this function to compute the distance between 2 genes
def ComputeDistanceBetweenGenes(gene1, gene2, GeneCoord):
    '''
    (str, str, dict) -> int
    Take 2 genes and the dictionary of gene coordinates and return the distance
    betwen the 2 genes
    Precondition: the 2 genes are located on the same chromosome
    '''
    
    # check that genes are on the same chromosome
    assert GeneCoord[gene1][0] == GeneCoord[gene2][0]    
    # get start positions
    S1, S2 = GeneCoord[gene1][1], GeneCoord[gene2][1]
    # get end positions
    E1, E2 = GeneCoord[gene1][2], GeneCoord[gene2][2]
    if S1 <= S2:
        D = S2 - E1
    elif S2 < S1:
        D = S1 - E2
    # return 0 for overlapping genes
    if D < 0:
        D = 0
    return D

  
# use this function to find overlapping genes
def FindOverlappingGenePairs(GeneCoord, OrderedGenes):
    '''
    (dict, dict) -> dict
    Take a dictionary of gene coordinates per chromosome and a dictionary of
    ordered genes along each chromosome and return a dictionary of host protein
    coding genes pairs that overlap    
    '''
    
    # GeneCoord is in the form {chromo: {gene:[chromo, start, end, sense]}}
    # OrderedGenes is a dict in the form {chromo: [gene1, gene2, gene3...]}    
        
    # create a dict to store all overlapping gene pairs 
    # key is always the fisrt gene along chromo {gene1: [gene2, gene3]}
    Overlap = {}
    # loop over each chromosome
    for chromo in GeneCoord:
        # loop over ordered genes on chromo:
        for i in range(len(OrderedGenes[chromo])-1):
            for j in range(i+1, len(OrderedGenes[chromo])):
                # compare each gene until the next gene doesn't overlap
                gene1, gene2 = OrderedGenes[chromo][i], OrderedGenes[chromo][j]
                # check for overlap between the 2 genes
                coordinates1 = set(range(GeneCoord[chromo][gene1][1], GeneCoord[chromo][gene1][2]))
                coordinates2 = set(range(GeneCoord[chromo][gene2][1], GeneCoord[chromo][gene2][2]))
                if len(coordinates1.intersection(coordinates2)) != 0:
                    # genes overlap, record gene pairs
                    # check if gene1 already in dict
                    if gene1 in Overlap:
                        Overlap[gene1].append(gene2)
                    else:
                        Overlap[gene1] = [gene2]
                else:
                    # exit loop to skip to the next gene, no need to evaulate all gene pairs
                    break
    return Overlap


# use this function to identify gene pairs in which a gene is fully contained in another gene
def FindContainedGenePairs(GeneCoord, Overlap):
    '''
    (dict, dict) -> dict
    Take the dictionary with gene coordinate, the dictionary with overlapping
    gene pairs and return a dictionary with gene pairs in which one gene is
    fully contained in another gene
    '''
 
    # GeneCoord is in the form {gene:[chromo, start, end, sense]}
    # Overlap is in the form {gene1: [gene2, gene3]}
        
    # create a dict with gene containing other genes {containing: [contained1, contained2]}
    ContainingGenes = {}
    # loop over each overlapping gene
    for gene1 in Overlap:
        # check if each gene in the overlapping gene is fully located within the over gene 
        for gene2 in Overlap[gene1]:
            # check if one of the 2 genes is fully contained in the other gene
            coord1 = set(range(GeneCoord[gene1][1], GeneCoord[gene1][2])) 
            coord2 = set(range(GeneCoord[gene2][1], GeneCoord[gene2][2]))
            FullyContained = False            
            if coord1.issubset(coord2):
                # gene1 is contained within gene2
                contained, containing = gene1, gene2
                # update boolean, found a gene within the second gene
                FullyContained = True
            elif coord2.issubset(coord1):
                # gene2 is contained within gene1
                contained, containing = gene2, gene1
                # update boolean, found a gene within the second gene
                FullyContained = True
            # check if one of the 2 genes is fully contained in the other gene
            if FullyContained == True:
                # populate dict of containing/contained gene pairs
                if containing in ContainingGenes:
                    ContainingGenes[containing].append(contained)
                else:
                    ContainingGenes[containing] = [contained]
    return ContainingGenes



# use this function to collapse all introns/CDS of all transcripts for a given gene
def CombineAllGeneRegions(GeneRegionCoord, TranscriptToGene):
    '''
    (dict, dict) -> dict
    Take the dictionary of CDS, exon or intron coordinates for each transcript
    the dictionary of transcript: gene pairs and return a dictionary with gene
    name and the coordinates of all gene region (intron, exon or CDS) combined
    for all transcripts of that gene
    '''
    # GeneRegionCoord in the form {transcript:[(region_start, region_end), (region_start, region_end)]}
    # TranscriptGeneMatches in the form {transcript: gene}

    # create a dict {gene: [(region_start, region_end), (region_start, region_end)]}
    AllGeneRegions = {}
    for TS in GeneRegionCoord:
        # get gene name
        gene = TranscriptToGene[TS]
        # populate dict with gene name : region coord
        # add coordinates for that transcript, transform in tuple
        if gene not in AllGeneRegions:
            AllGeneRegions[gene] = []
        for item in GeneRegionCoord[TS]:
            AllGeneRegions[gene].append(tuple(item))
    # collapse all identical coordinates
    for gene in AllGeneRegions:
        # transform coordinates into set to remove duplicate regions
        AllGeneRegions[gene] = set(AllGeneRegions[gene])
        # transform back into list to order coordinates
        AllGeneRegions[gene] = list(AllGeneRegions[gene])
        # transform all coordinates into lists
        for i in range(len(AllGeneRegions[gene])):
            AllGeneRegions[gene][i] = list(AllGeneRegions[gene][i])
        # order coordinates by start positions
        AllGeneRegions[gene].sort()
    return AllGeneRegions
   


# use this function to filter out contained genes that do share exons or introns with their host genes
def FindContainedGenesSharingExonIntron(ContainedGenes, IntronicCoord, ExonicCoord, GeneCoord):
    '''
    (dict, dict, dict, dict) -> dict
    Take the dictionary of genes with genes fully located within them,
    the dictionary of intron coordinates for each gene, the dictionary of exon 
    coordinates for each gene, the genomic coordinates of each gene and return
    a dictionary of host genes and nested genes with overlapping exons and/or introns    
    '''
    
    # ContainedGenes in the form {containing gene: [list of contained genes]}
    # IntronicCoord is in the form {gene: [list of intron coord]}
    # ExonicCoord is in the form {gene: [list of exon coord]}    
    # GeneCoordinates is the form {gene:[chromosome, start, end, sense]}   
   
    # create a dict {host: {contained genes}    
    HostGenes = {}
    for gene in ContainedGenes:
        # check if the contained genes share exon and/or introns with the host gene
        for contained in ContainedGenes[gene]:
            assert GeneCoord[gene][0] == GeneCoord[contained][0], 'host and nested have different chromo'
            # set up boolean to be updated if the contained share exon/intron
            ShareFeature = False
            # check if the contained gene has introns
            if contained in IntronicCoord:
                # loop over of introns of the contained gene
                for containedintron in IntronicCoord[contained]:
                    # get coordinates of the introns of the contained gene
                    containedintroncoord = set(range(containedintron[0], containedintron[1]))
                    # check if host gene has introns
                    if gene in IntronicCoord:
                        # loop over intron coordinates of the containing gene
                        for intron in IntronicCoord[gene]:
                            introncoord = set(range(intron[0], intron[1]))
                            # check if introns overlap
                            if len(containedintroncoord.intersection(introncoord)) != 0:
                                # update boolean
                                ShareFeature = True
                    # loop over exon coordinates of the host gene
                    for exon in ExonicCoord[gene]:
                        exoncoord = set(range(exon[0], exon[1]))
                        # check if intron and exon overlap
                        if len(containedintroncoord.intersection(exoncoord)) != 0:
                            # update boolean
                            ShareFeature = True
            # loop over exons of the contained gene
            for containedexon in ExonicCoord[contained]:
                # get coordinates of the exons of the contained gene
                containedexoncoord = set(range(containedexon[0], containedexon[1]))
                # check if host gene has introns
                if gene in IntronicCoord:
                    # loop over intron coordinates of the containing gene
                    for intron in IntronicCoord[gene]:
                        introncoord = set(range(intron[0], intron[1]))
                        # check if exon and intron overlap
                        if len(containedexoncoord.intersection(introncoord)) != 0:
                            # update boolean
                            ShareFeature = True
                # loop over exon coordinates of the host gene
                for exon in ExonicCoord[gene]:
                    exoncoord = set(range(exon[0], exon[1]))
                    # check if exons overlap
                    if len(containedexoncoord.intersection(exoncoord)) != 0:
                        # update boolean
                        ShareFeature = True       
            # check if contained genes overlaps with feature
            if ShareFeature == True:
                # record host: nested gene pairs
                if gene not in HostGenes:
                    HostGenes[gene] = set()
                HostGenes[gene].add(contained)
    # convert sets to lists
    for gene in HostGenes:
        HostGenes[gene] = list(HostGenes[gene])
    return HostGenes
    

# use this function to identify nested genes contained in the intron of their host genes
def FindIntronicNestedGenePairs(ContainedGenes, IntronicCoord, GeneCoord):
    '''
    (dict, dict, dict) -> dict
    Take the dictionary of genes with genes fully located within them,
    the dictionary of intron coordinates for each gene, the genomic coordinates
    of each gene and return a dictionary of host genes with nested genes fully
    contained within introns    
    '''
    
    # ContainedGenes in the form {containing gene: [list of contained genes]}
    # IntronicCoord is in the form {gene: [list of intron coord]}
    # GeneCoordinates is the form {gene:[chromosome, start, end, sense]}   
   
    # create a dict {host: {intronic nested genes}    
    HostGenes = {}
    for gene in ContainedGenes:
        # check if gene has introns
        if gene in IntronicCoord:
            # check if the contained genes are located within introns of the containing gene
            # get coordinates of each contained gene
            for contained in ContainedGenes[gene]:
                containedcoord = set(range(GeneCoord[contained][1], GeneCoord[contained][2]))
                # loop over intron coordinates of the containing gene
                for intron in IntronicCoord[gene]:
                    introncoord = set(range(intron[0], intron[1]))
                    # check if contained gene resides within intron
                    if containedcoord.issubset(introncoord):
                        # record host: nested gene pairs
                        if gene not in HostGenes:
                            HostGenes[gene] = set()
                        HostGenes[gene].add(contained)
    # convert sets to lists
    for gene in HostGenes:
        HostGenes[gene] = list(HostGenes[gene])
    return HostGenes


# use this function to generate pairs of host and nested genes
def GetHostNestedPairs(HostNested):
    '''
    (dict) -> list
    Take the dictionary woth host and nested genes and return a list of lists
    each containing a pair of host and nested gene
    '''
    
    # create a list of lists with [[host, nested]]
    HostNestedPairs = []
    for host in HostNested:
        for nested in HostNested[host]:
            pair = [host, nested]
            HostNestedPairs.append(pair)
    return HostNestedPairs


# use this script to remove gene pairs from the higher hierarchical level
def RemoveGenePairsFromHigherLevel(HigherLevelPairs, LowerLevelPairs):
    '''
    (list, list) -> dict
    Take the list of gene pairs from a higher inclusive level (ie. overlapping, 
    contained) and remove the gene pairs from the lower lovel (contained,
    intronic nested) to generate pairs of genes that do not overlapp and that 
    correspond to a unique gene organization
    Precondition: first gene in pair is the first gene on chromo    
    '''
    
    # create a list of gene pairs to remove
    to_remove = []
    # loop over gene pairs from lower hierarchical level
    for i in range(len(LowerLevelPairs)):
        # loop over gene pairs from higher hierarchical level 
        for j in range(len(HigherLevelPairs)):
            # compare the 2 gene pairs, ignore gene order
            if set(LowerLevelPairs[i]) == set(HigherLevelPairs[j]):
                # remove pair from higher hierarchical level 
                to_remove.append(HigherLevelPairs[j])
    # remove gene pairs from higher hierarchical level    
    for pair in to_remove:
        HigherLevelPairs.remove(pair) 
    return HigherLevelPairs
    
      
# use this function to get the orientation of a pair of genes
def GenePairOrientation(GenePair, GeneCoord):
    '''
    (list, dict) -> list
    Take a pair of gene, a dictionary with the genes coordinates and return a 
    list with the strand orientation of the 2 genes
    '''
    return [GeneCoord[GenePair[0]][-1], GeneCoord[GenePair[1]][-1]] 
  
  
# use this function to get the proportions of same and opposite strand pairs  
def GetSameOppositeStrandProportions(GenePairs, GeneCoord):
    '''
    (list, dict) -> (float, float)
    Take the list of gene pairs, the dictionary of gene coordinates and return
    the proportions of gene pairs with same strand and with opposite strand
    Precondition: genes without valid mRNAs have been filtered out 
    '''
  
    # count the number of same strand and opposite strand pairs
    same, opposite = 0, 0
    for pair in GenePairs:
        orientation = set(GenePairOrientation(pair, GeneCoord))
        if len(orientation) == 1:
            same += 1
        else:
             assert len(orientation) == 2, 'there should be only 2 different signs'
             opposite += 1
    # return proportions
    return (same / (same + opposite), opposite / (same + opposite))
  
  
  
## use this function to match human genes with orrthologs in a single other species
#def MatchOrthologPairs(OrthoFile):
#    '''
#    (file) -> dict
#    Take a file with orthology assignment between 2 species and return a dictionary
#    of 1:1 orthologs between the 2 species. 
#    '''
#    # the file colums foloww the format:
#    # Sp1GeneID Sp1TranscriptID Sp2GeneID Sp2GeneName HomologyType OrthologyC
#    
#    # create a dict of orthologs
#    Orthos = {}
#    infile = open(OrthoFile)
#    infile.readline()
#    for line in infile:
#        # consider only 1:1 orthologs
#        if 'ortholog_one2one' in line:
#            line = line.rstrip().split('\t')
#            # get gene IDs of the 2 species
#            gene1, gene2 = line[0], line[2]
#            # check that genes are ensembl gene IDs
#            for i in [gene1, gene2]:
#                assert 'ENS' in i, 'gene id is not valid'
#                assert 'ortholog' in line[4], 'ortholog should be in homology type'
#            # orthologous gene names appear multiple times in file because of multiple transcripts
#            if gene1 in Orthos:
#                assert gene2 == Orthos[gene1], 'gene is already matched to a 1:1 ortholog'
#            Orthos[gene1] = gene2
#    infile.close()                      
#    # check that all orthologs are 1;1 orthologs
#    values = list(Orthos.values())
#    for i in values:
#        assert values.count(i) == 1
#    return Orthos
 
# use this function to match human genes with orrthologs in a single other species
def MatchOrthologs(OrthoFile):
    '''
    (file) -> dict
    Take a file with orthology assignment between 2 species and return a dictionary
    orthologs between the 2 species. 
    '''
    # the file colums foloww the format:
    # Sp1GeneID Sp1TranscriptID Sp2GeneID Sp2GeneName HomologyType OrthologyC
    
    # create a dict of orthologs
    Orthos = {}
    infile = open(OrthoFile)
    infile.readline()
    for line in infile:
        # consider all orthology types
        if 'ortholog' in line:
            line = line.rstrip().split('\t')
            # get gene IDs of the 2 species
            # note that format of ortho file with Platypus isdifferent           
            if 'Platypus'in OrthoFile:
                gene1, gene2 = line[0], line[3]
            else:
                gene1, gene2 = line[0], line[2]
            # check that genes are ensembl gene IDs
            for i in [gene1, gene2]:
                assert 'ENS' in i, 'gene id is not valid'
                assert 'ortholog' in line[4], 'ortholog should be in homology type'
            # orthologous gene names appear multiple times in file because of multiple transcripts
            if gene1 in Orthos:
                Orthos[gene1].add(gene2)
            else:
                Orthos[gene1] = set()
                Orthos[gene1].add(gene2)
    infile.close()                      
    # convert set to lists
    for gene in Orthos:
        Orthos[gene] = list(Orthos[gene])
    return Orthos

# Map humangenes to their orthologs in 2 other species  
def MatchOrthologTrios(OrthoFile):
    '''
    (file) -> dict
    Take a file with orthology assignment between 3 species and return a dictionary
    of orthologs between the 3 species. Consider only 1:1 orthologs, use gene ID
    of 1st species as key
    '''
    
    # create a dict of orthologs
    Orthos, Orthologs = {}, {}
    
    infile = open(OrthoFile)
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            # check that there is a ortholog in both species
            if 'ortholog' in line:
                if line.count('ortholog') == 2:
                    line = line.split('\t')
                    assert len(line) == 10
                    # get gene IDs of the 3 species
                    gene1, gene2, gene3 = line[0], line[2], line[6]
                    # check that genes are ensembl gene IDs
                    for i in [gene1, gene2, gene3]:
                        assert 'ENS' in i, 'gene id is not valid'
                    assert 'ortholog' in line[4] and 'ortholog' in line[8], 'ortholog should be in homology type'
                    # record 1:1 orthologs
                    if line[4] == 'ortholog_one2one' and line[8] == 'ortholog_one2one':
                        if gene1 not in Orthos:
                            Orthos[gene1] = [set(), set()]
                            Orthos[gene1][0].add(gene2)
                            Orthos[gene1][1].add(gene3)
                        else:
                            Orthos[gene1][0].add(gene2)
                            Orthos[gene1][1].add(gene3)
    infile.close()                      
 
    # check that all orthologs are 1;1 orthologs
    # make a dict {gene1: [gene2, gene3]}
    for gene in Orthos:
        assert len(Orthos[gene][0]) == len(Orthos[gene][1]) == 1, 'there is more than 1 ortholog'
        Orthos[gene][0] = list(Orthos[gene][0])
        Orthos[gene][1] = list(Orthos[gene][1])        
        Orthologs[gene] = [Orthos[gene][0][0], Orthos[gene][1][0]]
    print(len(Orthos), len(Orthologs))
    #return Orthologs  
    return Orthologs



# use this function to parse the primate expression data file
def ParsePrimateExpressionData(ExpressionFile, species):
    '''
    (file) -> dict
    Take the file with expression (RPKM) for each gene and individual and tissu
    in primates and return a dictionary for a given species with median expression
    across individuals for each of the tissues. (note that orangutan has no
    expression in testis)
    '''
    
    # create a lambda function to convert the string values into floats
    ExpVal = lambda x: float(x)
    
    infile = open(ExpressionFile)
    infile.readline()
    
    # create a dict with medium expression across individuals for the different organs
    # {gene_ID: [brain, cerebellum, heart, kidney, liver, testis]}    
    expression = {}
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get gene ID, extract the expression values for each organ
            if species == 'human':
                gene = line[0]
                brexp, cbexp, htexp = list(map(ExpVal, line[5:11])), list(map(ExpVal, line[11:13])), list(map(ExpVal, line[13:16]))
                kdexp, lvexp, tsexp = list(map(ExpVal, line[16: 19])), list(map(ExpVal, line[19: 21])), list(map(ExpVal, line[21: 23]))
            elif species == 'chimp':
                gene = line[1]
                brexp, cbexp, htexp = list(map(ExpVal, line[23: 29])), list(map(ExpVal, line[29: 31])), list(map(ExpVal, line[31: 33]))
                kdexp, lvexp, tsexp = list(map(ExpVal, line[33: 35])), list(map(ExpVal, line[35: 37])), list(map(ExpVal, line[37:38]))
            elif species == 'gorilla':
                gene = line[2]
                brexp, cbexp, htexp = list(map(ExpVal, line[50: 52])), list(map(ExpVal, line[52: 54])), list(map(ExpVal, line[54: 56]))
                kdexp, lvexp, tsexp = list(map(ExpVal, line[56: 58])), list(map(ExpVal, line[58: 60])), list(map(ExpVal, line[60:61]))
            elif species == 'orangoutan':
                gene = line[3]
                brexp, cbexp, htexp = list(map(ExpVal, line[61: 63])), list(map(ExpVal, line[63:64])), list(map(ExpVal, line[64: 66]))
                kdexp, lvexp = list(map(ExpVal, line[66: 68])), list(map(ExpVal, line[68:70])) 
            elif species == 'macaque':
                gene = line[4]
                brexp, cbexp, htexp = list(map(ExpVal, line[70: 73])), list(map(ExpVal, line[73: 75])), list(map(ExpVal, line[75: 77]))
                kdexp, lvexp, tsexp = list(map(ExpVal, line[77: 79])), list(map(ExpVal, line[79: 81])), list(map(ExpVal, line[81:]))
            # get the median expression level per tissue
            assert gene not in expression, 'gene is already recorded'
            # check species (orang outan has no expression in testis)
            if species == 'orangoutan':
                expression[gene] = [np.median(brexp), np.median(cbexp), np.median(htexp), np.median(kdexp), np.median(lvexp)]
            else:
                expression[gene] = [np.median(brexp), np.median(cbexp), np.median(htexp), np.median(kdexp), np.median(lvexp), np.median(tsexp)]
    infile.close()
    return expression


# use this function to parse the GTEX or Encode median expression file
def ParseExpressionFile(ExpressionFile):
    '''
    (file) -> dict
    Take the file with median gene expression in each tissue from GTEX in human
    of from Encode in mouse and return a dictionary with gene: list of expression
    in tissues key, value pairs
    '''
    # create a dict to store the expression (median tissue expression FPKM normalized by upper quartiles))
    Expression = {}
    infile = open(ExpressionFile)
    header = infile.readline()
    for line in infile:
        if line.startswith('ENS'):
            line = line.rstrip().split('\t')
            # get gene name
            gene = line[0]
            # get expression profile
            profile = line[1:]
            # convert strings to float
            profile = list(map(lambda x: float(x), profile))
            Expression[gene] = profile
    infile.close()
    return Expression


# use this function to remove genes without expression
def RemoveGenesLackingExpression(ExpressionProfile):
    '''
    (dict) -> dict
    Take the dictionary of expression profile for each gene and return a modified
    dictionary excluding genes that are not expressed in any tissue
    '''
    
    # create a list of genes to remove
    to_remove = [gene for gene in ExpressionProfile if sum(ExpressionProfile[gene]) == 0]    
    for gene in to_remove:
        del ExpressionProfile[gene]
    return ExpressionProfile

# use this function to transform absulte expression into relative expression
def TransformRelativeExpression(ExpressionProfile):
    '''
    (dict) -> dict
    Take the dictionary with absolute expression expression level in each tissue
    for each gene and return a dictionary with gene as key and list of relative
    expression as value
    '''
    # create a dict {gene: [list of relative expression]}
    RelativeExpression = {}
    for gene in ExpressionProfile:
        # initialize value with empty list
        RelativeExpression[gene] = []
        # transform the expression level at each stage
        for i in range(len(ExpressionProfile[gene])):
            # divide absolue level by sum of absolute expressions
            RelativeExpression[gene].append(ExpressionProfile[gene][i] / sum(ExpressionProfile[gene]))
    return RelativeExpression


# use this function to compute the expression breadth of expressed genes
def ExpressionBreadth(ExpressionProfile):
    '''
    (dict) -> dict
    Take the dictionary with expression profiles for each gene and return a 
    dictionary with the number of tissues in which the gene is expressed
    Precondition: genes without expression have been removed    
    '''
    
    # create a dictionary {gene: number of tissues with expression}
    breadth = {}
    # loop over genes with expression profiles
    for gene in ExpressionProfile:
        # count the number of tissues with non-0 expression
        Ntissues = [i for i in ExpressionProfile[gene] if i != 0]
        Ntissues = len(Ntissues)
        breadth[gene] = Ntissues
    return breadth


# use this function to compute expression specificity of a single gene
def ComputeTau(L):
    '''
    (list) -> float
    Take a list with expression level values of a single gene and return the
    expression specificity (tau) index of this gene. tau ranges from 0 to 1 and
    high values indicates tissue-specific expression     
    '''
    
    tau = 0
    # loop over list of expression values, add 1 - (value / maximum value)
    for i in range(len(L)):
        tau += (1 - (L[i] / max(L)))
    # divide tau by number of tissues -1
    return tau / (len(L) -1)

# use this function to compute expression specificity of each gene
def ExpressionSpecificity(ExpressionProfile):
    '''
    (dict) -> dict
    Take the dictionary with expression profiles for each gene and return a 
    dictionary with expression specificity
    Precondition: genes without expression have been removed
    '''
    
    # create a dict {gene: tau}
    specificity = {}
    # loop over genes with expression profiles
    for gene in ExpressionProfile:
        specificity[gene] = ComputeTau(ExpressionProfile[gene])
    return specificity

# use this function to compute the euclidian distance between expression vectors
def EuclidianDistance(L1, L2):
    '''
    (list, list) -> float
    Take 2 lists of expression levels (lists have the same length) and return 
    the euclidian distance between the two vectors
    '''
    # set up distance
    D = 0
    # loop over list 1
    for i in range(len(L1)):
        D = D + (L1[i] - L2[i])**2
    # take the square root
    D = math.sqrt(D)
    return D


# use this function to compute the expression distance (euclidian distance) among gene pairs
def ComputeExpressionDivergenceGenePairs(L, ExpressionProfile):
    '''
    (list, dict) -> list
    Take a list of inner lists with gene pairs, the dictionary of expression
    profiles for each gene and return a list of euclidian distances measuring
    the expression divergence between the 2 genes in each pair    
    '''
    # create a list to store the expression divergence
    Divergence = []
    # loop over the list of gene pairs
    for i in range(len(L)):
        # compute divergence between the 2 genes in the given pair
        gene1, gene2 = L[i][0], L[i][1]        
        D = EuclidianDistance(ExpressionProfile[gene1], ExpressionProfile[gene2])
        Divergence.append(D)
    return Divergence


# use this function to compute the expression distance (euclidian distance) among pairs of orthologs
def ComputeExpressionDivergenceOrthologs(L, ExpressionProfileSp1, ExpressionProfileSp2):
    '''
    (list, dict, dict) -> list
    Take a list of inner lists with pairs of orthologs and return a list of euclidian
    distances measuring the expression divergence between the 2 genes in each pair    
    '''
    # create a list to store the expression divergence
    Divergence = []
    # loop over the list of gene pairs
    for i in range(len(L)):
        # compute divergence between the 2 genes in the given pair
        gene1, gene2 = L[i][0], L[i][1]        
        D = EuclidianDistance(ExpressionProfileSp1[gene1], ExpressionProfileSp2[gene2])
        Divergence.append(D)
    return Divergence


# use this function to create lists of orthologs with both genes expressed
def ExpressedOrthologousPairs(Sp1Expression, Sp2Expression, Genes, Orthologs):
    '''
    (dict, dict, set, dict) -> list
    Take the dictionaries of expression profiles for species 1 and 2, the set
    of genes of interest in species 1, and the dictionary of orthologs and return
    a list of expressed orthologous pairs
    '''
    # create a list of gene pairs
    ExpressedOrthos = []
    # loop over gene set of interest
    for gene in Genes:
        # check that gene has ortholog
        if gene in Orthologs:
            # check that gene and its orthologs are expressed
            if gene in Sp1Expression and Orthologs[gene] in Sp2Expression:
                ExpressedOrthos.append([gene, Orthologs[gene]])
    return ExpressedOrthos


# use this function to generate sets of gene pairs separated by a given distance
def GenerateSetsGenePairsDistance(GeneCoord, OrderedGenes, ExpressionProfile):
    '''
    (dict, dict, dict) -> tuple
    Take the dictionary of gene coordinates, the dictionary of ordered genes
    along each chromosome, the dictionary of gene: expression pairs and return a
    tuple with lists of gene pairs separated by a given distance (in bp):
    proximal, intermediate and distant
    '''
        
    # make lists of gene pairs [[gene1, gene2], ....[gene n, gene n+1]]
    Proximal, Moderate, Intermediate, Distant = [], [], [], []
    # loop over chromosomes
    for chromo in OrderedGenes:
        # loop over the list of ordered genes
        for i in range(len(OrderedGenes[chromo])):
            # check that gene is not host or nested, has expression
            if OrderedGenes[chromo][i] in ExpressionProfile:
                # get the end position of gene 1
                EndGene1 = GeneCoord[OrderedGenes[chromo][i]][2]                
                # grab 2nd gene to form a pair                
                for j in range(i+1, len(OrderedGenes[chromo])):
                    # check that gene is not host or nested and has expression
                    if OrderedGenes[chromo][j] in ExpressionProfile:
                        # get the start position of gene 2
                        StartGene2 = GeneCoord[OrderedGenes[chromo][j]][1]
                        # check if distance is less that 500 bp
                        D = StartGene2 - EndGene1
                        if D >= 0 and D < 1000:
                            # add gene pair to Proximal
                            Proximal.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
                        elif D >= 1000 and D < 10000:
                            # add gene pair to Intermediate
                            Moderate.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
                        elif D >= 10000 and D < 50000:
                            # add gene pair to Intermediate
                            Intermediate.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
                        elif D >= 50000:
                            # add gene pair to Distant
                            Distant.append([OrderedGenes[chromo][i], OrderedGenes[chromo][j]])
    return Proximal, Moderate, Intermediate, Distant


# use this function to remove nested gene pairs without expression
def FilterGenePairsWithoutExpression(HostNestedPairs, ExpressionProfile, Stringency):
    '''
    (list, dict, str) -> list
    Take the list of gene pairs and the dictionary of expression profile for
    each gene, the level of stringency to remove gene pairs without expression
    and return a modified list of gene pairs in which pairs of gene
    lacking expression are removed
    Precondition: genes without expression in the expression profile dictionary
    have been removed    
    '''
    
    # filter gene pairs based on expression
    to_remove = []
    for pair in HostNestedPairs:
        if Stringency == 'loose':
            # remove gene pair if both genes are not expressed
            if pair[0] not in ExpressionProfile and pair[1] not in ExpressionProfile:
                to_remove.append(pair)
        elif Stringency == 'strict':
            # remove gene pair if any gene is not expressed
            if pair[0] not in ExpressionProfile or pair[1] not in ExpressionProfile:
                to_remove.append(pair)
    # remove pairs        
    for pair in to_remove:
        HostNestedPairs.remove(pair)
    return HostNestedPairs

# use this function to make a set of genes with full or partial overlap
def MakeFullPartialOverlapGeneSet(OverlapGenes):
    '''
    (dict) -> set
    Take the dictionary of genes that are fully or partially overlapping 
    (ie, overlapping or nested) and return a set of genes with coordinates that
    overlap partially or fully 
    '''
    # OverlapGenes is in the form {gene: [gene1, gene2]}
    WithOverlap = set()
    for i in OverlapGenes:
        WithOverlap.add(i)
        for j in OverlapGenes[i]:
            WithOverlap.add(j)
    return WithOverlap

# use this function to make a set of non-overlapping genes
def MakeNonOverlappingGeneSet(OverlapGenes, GeneCoord):
    '''
    (dict, dict) -> set     
    Take the dictionary of genes that are overlapping and the dictionary of
    gene coordinates and return a set of genes that are not overlapping
    '''
    # make a set of genes with overlapping coordinates
    WithOverlap = MakeFullPartialOverlapGeneSet(OverlapGenes)
    # make a set of all genes
    AllGenes = set(GeneCoord.keys())
    # return a set of genes that do not overlap
    return AllGenes - WithOverlap 
    

# use this function to match the longest transcript of the nested genes to transcripts of the host genes
def MatchHostTranscriptWithNestedTranscript(HostGenes, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoordinates):
    '''
    (dict, dict, dict, dict, dict) -> int
    Take the dictionary of host: nested genes, the dictionary of gene and
    transcript pairs, the dictionary of gene: longest transcript pairs,
    the dictionary of transcript coordinates, the dictionary of intron positions
    for each transcript and return a dictionary of host transcripts (in priority
    the longest) with the longest transcript of their nested genes
    '''

    # HostGenes is in the form {gene: [list of nested genes]}
    # MapGeneTranscript is in the form {gene: [transcripts]}
    # GeneLongestTranscript is the form {gene: longest_transcript}
    # TranscriptCoordinates is in the form {transcript: [chromo, start, end, orientation]}
    # IntronCoordinates if the form {transcript: [(intron_start, intron_end), ...]}  
    
    # create a dict {host_TS: [nested_TS]}
    HostNestedTranscripts = {}    
    # create a dict to store the host: nested gene pairs
    # for which the longest nested transcript is not in the longest host transcript    
    CheckAgain = {}
    
    # loop over nested genes
    # check if the longest transcript of the nested gene is in the longest
    # transcript of the host gene
    for gene in HostGenes:
        # get the gene's longest transcript
        LongestHost = GeneLongestTranscript[gene]
        # loop over the gene's nested genes
        for nested in HostGenes[gene]:
            # get the longest transcript of the nested gene
            LongestNested = GeneLongestTranscript[nested]
            # get the coordinates of the longest transcript of the nested gene
            nestedcoord = set(range(TranscriptCoordinates[LongestNested][1], TranscriptCoordinates[LongestNested][2]))
            # set up boolean, update when nested transcript is found in the host transcript            
            FoundNested = False
            # loop over the intron in the longest transcript of the host gene
            for intron in IntronCoordinates[LongestHost]:
                introncoord = set(range(intron[0], intron[1]))
                # check if transcript is within the intron
                if nestedcoord.issubset(introncoord):
                    # nested transcript is contained in the intron
                    # match nested transcript to host transcript, exit loop
                    if LongestHost not in HostNestedTranscripts:
                        HostNestedTranscripts[LongestHost] = []
                    HostNestedTranscripts[LongestHost].append(LongestNested)
                    FoundNested = True
                    break
            # record host gene: nested gene pair if not found            
            if FoundNested == False:
                if gene not in CheckAgain:
                    CheckAgain[gene] = []
                CheckAgain[gene].append(nested)
    
    # check if some nested genes have not been matched
    if len(CheckAgain) != 0:
        # create a set of nested genes that have already been matched
        AlreadyMatched = set()
        # loop over the host genes
        for gene in CheckAgain:
            # loop over all transcripts of the host gene
            for TS in MapGeneTranscript[gene]:
                # check if host transcript has introns
                if TS in IntronCoordinates:
                    # do not consider the longest transcript
                    if TS != GeneLongestTranscript[gene]:
                        # loop over the nested genes
                        for nested in CheckAgain[gene]:
                            # check that nested gene is not already recorded
                            if nested not in AlreadyMatched:
                                # always grab the longest transcript of the nested gene
                                LongestNested = GeneLongestTranscript[nested]
                                # get the coordinates of the longest transcript of the nested gene
                                nestedcoord = set(range(TranscriptCoordinates[LongestNested][1], TranscriptCoordinates[LongestNested][2]))
                                # set up boolean, update when nested transcript is found in the host transcript            
                                FoundNested = False
                                # loop over the introns of the host transcript
                                for intron in IntronCoordinates[TS]:
                                    introncoord = set(range(intron[0], intron[1]))
                                    # check if transcript is within the intron
                                    if nestedcoord.issubset(introncoord):
                                        # nested transcript is contained within the intron
                                        # match nested transcript to host transcript, exit loop
                                        if TS not in HostNestedTranscripts:
                                            HostNestedTranscripts[TS] = []
                                        HostNestedTranscripts[TS].append(LongestNested)
                                        FoundNested = True
                                        # nested gene is recorded
                                        AlreadyMatched.add(nested)
                                        break
        # check that all nested genes (and therefore all host transcripts) have been matched
        nestedgenes = set()
        for gene in CheckAgain:
            for i in CheckAgain[gene]:
                nestedgenes.add(i)
        assert nestedgenes == AlreadyMatched, 'not all nested genes have been matched'
    return HostNestedTranscripts


# use this function to collect the number of introns for host, nested and un-nested genes 
def CollectIntronNumbers(UnNestedGenes, HostNestedTranscriptMatches, GeneLongestTranscript, TranscriptCoordinates, IntronCoordinates):
    '''
    (set, dict, dict, dict, dict) -> (list, list, list)
    Take the dictionary of host: nested gene genes, the dictionary of host
    transcripts: longest nested transcripts pairs, the dictionary of gene coordinates,
    the dictionary of gene: longest transcript pairs, the dictionary of intron
    positions for each transcript and the dictionary of transcript coordinates 
    and return a tuple with lists containing the number of introns for
    the transcript of the host, nested and un-nested genes    
    '''
    # GeneLongestTranscript is in the form {gene: longest_transcript}
    # HostNestedTranscriptMatches if in the form {host_TS: [nested_TS]}
    # GeneCoord is in the form {gene: [chromo, start, end, orientation]}
    # GeneLongestTranscript is in the form {gene: longest_transcript}
    # IntronCoordinates is in the form {transcript: [(intron_start, intron_end), ...]}
    # TranscriptCoordinates is in the form {transcript:[chromosome, start, end, sense]}

    # initialize empty list to store the number of itnrons
    HostNum, NestedNum, OthersNum = [], [], []
    # loop over host transcripts 
    for hostTS in HostNestedTranscriptMatches:
        # check that this transcript has intron
        assert hostTS in IntronCoordinates, 'transcript of host gene should have introns'
        # count the number of introns        
        HostNum.append(len(IntronCoordinates[hostTS]))
        # loop over nested transcripts
        for nestedTS in HostNestedTranscriptMatches[hostTS]:
            # check that transcript has introns
            if nestedTS in IntronCoordinates:
                # count the number of introns
                NestedNum.append(len(IntronCoordinates[nestedTS]))
            else:
                # check that transcript has coordinates
                assert nestedTS in TranscriptCoordinates, 'nested TS should have coordinates'
                # transcript is intronless, add 0
                NestedNum.append(0)
    # loop over un-nested genes
    for gene in UnNestedGenes:
        # get the longest transcript
        TS = GeneLongestTranscript[gene]
        # check if gene has introns
        if TS in IntronCoordinates:
            # count the number of introns
            OthersNum.append(len(IntronCoordinates[TS]))
        else:
            # check that transcript has coordinates
            assert TS in TranscriptCoordinates, 'un-nested gene should have coordinates'
            # transcript is intronless, add 0
            OthersNum.append(0)
    return HostNum, NestedNum, OthersNum


# use this function to get the intron length of the host transcripts
def CollectHostGeneIntronLength(HostNestedTranscriptMatches, TranscriptCoordinates, IntronCoordinates):
    '''
    (dict, dict, dict) -> (list, list)
    Take the dictionary of host transcripts: nested transcripts, the dictionary
    of transcript coordinates, the dictionary of intron coordinates and return
    a tuple with lists of intron length for gene-containing introns (minus the
    length of all transcripts inside) and introns without genes.     
    '''
    # HostNestedTranscriptMatches if in the form {host_TS: [nested_TS]}
    # TranscriptCoordinates is in the form {transcript: [chromo, start, end, orientation]}
    # IntronCoordinates if the form {transcript: [(intron_start, intron_end), ...]}  
    
    # create lists to store intron length
    GeneContainingIntron, NoGeneIntron = [], []    
    # loop over host transcripts 
    for hostTS in HostNestedTranscriptMatches:
        # loop over the host transcript's introns
        for HostIntron in IntronCoordinates[hostTS]:
            IntronWithGene = False
            # get intron coordinates and intron length
            introncoord = set(range(HostIntron[0], HostIntron[1]))
            intronlength = set(range(HostIntron[0], HostIntron[1]))
            # loop over nested transcripts, remove length of transcript
            # length of any transcripts found in intron      
            for nestedTS in HostNestedTranscriptMatches[hostTS]:
                # get the coordinates of the nested transcript
                nestedcoord = set(range(TranscriptCoordinates[nestedTS][1], TranscriptCoordinates[nestedTS][2]))
                # check if nested transcript is contained in intron
                if nestedcoord.issubset(introncoord):
                    # update boolean if any transcript is found in intron
                    IntronWithGene = True
                    # remove the positions of the nested transcript from intron
                    intronlength -= nestedcoord 
            # check if intron contains a gene
            if IntronWithGene == True:
                GeneContainingIntron.append(len(intronlength))
            elif IntronWithGene == False:
                NoGeneIntron.append(len(intronlength))
    return GeneContainingIntron, NoGeneIntron


# use this function to get the intron length of nested transcripts
def CollectNestedGeneIntronLength(HostNestedTranscriptMatches, TranscriptCoordinates, IntronCoordinates):
    '''
    (dict, dict, dict, dict) -> list
    Take the dictionary of host transcripts: nested transcripts, the dictionary
    of transcript coordinates, the dictionary of intron coordinates and return
    a list of intron length for nested genes     
    '''
    # HostNestedTranscriptMatches if in the form {host_TS: [nested_TS]}
    # TranscriptCoordinates is in the form {transcript: [chromo, start, end, orientation]}
    # IntronCoordinates if the form {transcript: [(intron_start, intron_end), ...]}  
    
    # create a list of intron length
    IntronLength = []
    # loop over host transcripts 
    for hostTS in HostNestedTranscriptMatches:
        # loop over nested transcripts
        for nestedTS in HostNestedTranscriptMatches[hostTS]:
            if nestedTS in IntronCoordinates:
                # loop over the introns of the nested transcript
                for intron in IntronCoordinates[nestedTS]:
                    IntronLength.append(len(set(range(intron[0], intron[1]))))
            else:
                assert nestedTS in TranscriptCoordinates, 'transcript should have coordinates'
                IntronLength.append(0)
    return IntronLength


# use this function to get the intron length of un-nested transcripts
def CollectUnNestedGeneIntronLength(UnNestedGenes,  GeneLongestTranscript, TranscriptCoordinates, IntronCoordinates):
    '''
    (set, dict, dict, dict) -> list
    Take the set of un-nested genes, the dictionary of gene: longest transcript,
    the dictionary of transcript coordinates, the dictionary of intron coordinates
    and return a list of intron length for nested genes     
    '''
    # GeneLongestTranscript is in the form {gene: longest_transcript}
    # TranscriptCoordinates is in the form {transcript: [chromo, start, end, orientation]}
    # IntronCoordinates if the form {transcript: [(intron_start, intron_end), ...]}  
    
    # create a list of intron length
    IntronLength = []
    # loop over un-nested genes
    for gene in UnNestedGenes:
        # get the longest transcript
        TS = GeneLongestTranscript[gene]
        # check if gene has introns
        if TS in IntronCoordinates:
            # loop over introns of the transcript
            for intron in IntronCoordinates[TS]:
                IntronLength.append(len(set(range(intron[0], intron[1]))))
        else:
            # check that transcript has coordinates
            assert TS in TranscriptCoordinates, 'un-nested gene should have coordinates'
            # transcript is intronless, add 0
            IntronLength.append(0)
    return IntronLength


# use this function to generate un-nested genes to randomly draw
def GenerateAllUnNestedGenes(Overlap, OrderedGenes, ExpressionProfile):
    '''
    (set, dict, dict) -> dict
    Take the set of host and nested genes, the dictionary of ordered genes
    along each chromosome, the dictionary of gene: expression pairs and return
    a dictionary of pairs of number: un-nested gene on each chromsome
    '''
    
    # make a dictionary with chromsome as key and number: gene {chromo: {num: gene}}    
    ToDrawGenesFrom = {}
    # loop over chromosomes
    for chromo in OrderedGenes:
        # set up counter
        k = 0
        # add chromo as key and intialize inner dict
        ToDrawGenesFrom[chromo] = {}
        # loop over the list of ordered genes
        for i in range(len(OrderedGenes[chromo])):
            # check that gene does not overlap with any other gene and that gene is expressed
            if OrderedGenes[chromo][i] not in Overlap and OrderedGenes[chromo][i] in ExpressionProfile:
                # add gene pair and update counter
                ToDrawGenesFrom[chromo][k] = OrderedGenes[chromo][i]
                k += 1
    return ToDrawGenesFrom


# use this function to generate a pool of un-nested gene pairs randomly draw from
def GenerateMatchingPoolPairs(pair, ToDrawGenesFrom, GeneCoord, Distance):
    '''
    (list, dict, dict, int) -> list
    Take a gene pair gene, the dictionary of genes to draw from,
    the dictionary of gene coordinates and the matching distance and return a
    list of expressed and matching gene pairs (distance, orientation and
    chromosome) to randomly draw from
    '''
    
    PairPool = []
    gene1, gene2 = pair[0], pair[1]
    # get gene orientation
    orientation = set(GenePairOrientation(pair, GeneCoord))
    # get gene chromos
    chromo1, chromo2 = GeneCoord[gene1][0], GeneCoord[gene2][0]
    # compute distance between genes if genes are on the same chromosome
    if chromo1 == chromo2:
        D = ComputeDistanceBetweenGenes(gene1, gene2, GeneCoord)
    # check if genes are on the same chromosome
    if chromo1 != chromo2:
        # make lists of genes on each matching chromosome
        PossibleGenesChromo1 = [ToDrawGenesFrom[chromo1][i] for i in ToDrawGenesFrom[chromo1]]
        PossibleGenesChromo2 = [ToDrawGenesFrom[chromo2][i] for i in ToDrawGenesFrom[chromo2]]
        for i in range(len(PossibleGenesChromo1)):
            for j in range(len(PossibleGenesChromo2)):
                G1, G2 = PossibleGenesChromo1[i], PossibleGenesChromo2[j]
                # match gene orientation
                if {GeneCoord[G1][-1], GeneCoord[G2][-1]} == orientation:
                    PairPool.append([G1, G2]) 
    else:
        # make a list of genes on chromo
        PossibleGenes = [ToDrawGenesFrom[chromo1][i] for i in ToDrawGenesFrom[chromo1]]
        for i in range(0, len(PossibleGenes) -1):
            for j in range(i+1, len(PossibleGenes)):
                G1, G2 = PossibleGenes[i], PossibleGenes[j]
                # check that genes have matching orientation
                if {GeneCoord[G1][-1], GeneCoord[G2][-1]} == orientation:
                    # compute distance between genes
                    d = ComputeDistanceBetweenGenes(G1, G2, GeneCoord)
                    # check if distance is within limits
                    if D - Distance <= d <= D + Distance:
                        PairPool.append([G1, G2])
    return PairPool   
    
 
 
# use this function to determine if at least 1 ortholog is expressed
def IsOrthoExpressed(gene, Orthologs, SecondSpExpression):
    '''
    (str, dict) -> bool
    Take a gene of interest, the dictionary of gene: orthologs and the dictionary
    of expression profile of the orthologs in a second species and return a boolean
    True if at least one ortholog is expressed and False if no orthologs are expresed
    '''
    OrthoExpressed = False   
    for ortho in Orthologs[gene]:
        if ortho in SecondSpExpression:
            OrthoExpressed = True
    return OrthoExpressed
    

# use this function to generate a list of control genes with matching characteristics
def GenerateMatchingGenes(GeneList, GeneCoord, ToDrawGenesFrom, ExpressionSpecificity, Orthologs, SecondSpExpression):
    '''
    (list, dict, dict, dict, dict) -> list
    Take a list of genes of interest, the dictionary of gene coordinates,
    the dictionary of number: un-nested and expressed gene pairs, the dictionary of 
    gene expression specificity and the dictionary of expression profiles for orthologs
    in the sister-species and return a list of matching (chromosome and tissue specificity) un-nested genes
    Precondition: nested and non-expressed genes have been removed from the dict num: gene pairs
    '''
    
    # make a list of control genes with matching proporties    
    ControlGenes = []
    # loop of list of genes of interest
    for gene in GeneList:
        # get chromo
        chromo = GeneCoord[gene][0]
        # make a list of matching genes
        MatchingGenes = []
        for i in ToDrawGenesFrom[chromo]:
            # check that gene has ortholog
            if ToDrawGenesFrom[chromo][i] in Orthologs:
                # check if at least one ortholog is expressed
                if IsOrthoExpressed(ToDrawGenesFrom[chromo][i], Orthologs, SecondSpExpression) == True:
                    # check that gene has matching expression specificity
                    if ExpressionSpecificity[gene] - 0.01 <= ExpressionSpecificity[ToDrawGenesFrom[chromo][i]] <= ExpressionSpecificity[gene] + 0.01:
                        MatchingGenes.append(ToDrawGenesFrom[chromo][i])       
        if len(MatchingGenes) != 0:
            # pick a random gene in list of matching genes
            i = random.randint(0, len(MatchingGenes) -1)
            # add matching gene to list of control genes
            ControlGenes.append(MatchingGenes[i])
    return ControlGenes
 

# use this function to dtermine if a single nested pair is conserved in a second species
def IsNestedPairConserved(pair, Orthologs, OrthoNestedPairs):
    '''
    (list, dict, list) -> bool
    Take a list of host, nested genes in a focal species, the dictionary of
    gene orthologs between focal and second species, and the list of host, nested
    lists in the second species
    Precondition: genes in pair have both orthologs in second species
    '''
    # pair is a list of host, nested genes in focal species [host, nested]
    # Orthologs is a dictionary of orthologs {gene: [ortho1, ortho2..]}    
    # OrthoNestedPairs is a list of [host, nested] gene pairs in other species
    
    # generate a list with all combinations of orthologs of host:nested pairs
    orthopairs = []
    for ortho1 in Orthologs[pair[0]]:
        for ortho2 in Orthologs[pair[1]]:
            orthopairs.append([ortho1, ortho2])
    # set boolean to be updated if nesting is conserved
    ConservedNesting = False
    # check if at least one ortholog pair is nested in sister species
    for i in orthopairs:
        if i in OrthoNestedPairs:
            # nesting is conserved, update boolean and exit loop 
            ConservedNesting = True
            break
    return ConservedNesting


# use this function to sort young and ancestral nesting events
def InferYoungOldNestingEvents(NestedGenes, SpeciesNestedPairs, OrthologPairs): 
    '''
    (list, list, list) -> (list, list)    
    Take a list of host, nested gene pair lists in the focal species, a list
    of lists of host, nested pairs in sister species and outgroups, and 
    a list of dictionaries of orthologs and return a tuple with list of
    host: nested pairs that are infered to be old and young
    (before the divergence of the 2 species or after)
    '''   
 
    # require presence of both genes in human, sister-species and at least 1 outgroup
    # young nesting if nesting not present any other species 
    # old nesting if nesting present in sister species
    
    # NestedGenes is the list of gene pairs in the species of interest (eg. human)
    # siters_species = chimp or mouse if species of interest is human    
    # SpeciesNestedPairs include list of nested gene pairs [sister_species, outgroups]
    # OrthologPairs include dictionary of orthologs between human and other species [sister_species, outgroups]

    assert len(SpeciesNestedPairs) == len(OrthologPairs)

    # create lists of young and old nested gene pairs in the focal species
    Young, Old = [], []
    
    # loop over gene pairs in species of interest
    for pair in NestedGenes:
        # record indices (species) in which orthologs of both genes are present
        PresentInSpecies = []
        for i in range(len(OrthologPairs)):
            if pair[0] in OrthologPairs[i] and pair[1] in OrthologPairs[i]:
                PresentInSpecies.append(i)
        # check that orthologs are present >= 1 outgroup and in sister species
        if len(PresentInSpecies) >= 2 and PresentInSpecies[0] == 0:
            # check if pair is present in sister species
            ConservedNesting = IsNestedPairConserved(pair, OrthologPairs[0], SpeciesNestedPairs[0])
            if ConservedNesting == True:
                Old.append(pair)
            else:
                # check if pair is absent from all species, including outgroups
                # loop over species in which orthologs are present
                ConservedInOutGroups = False
                for i in PresentInSpecies:
                    NestedConserved = IsNestedPairConserved(pair, OrthologPairs[i], SpeciesNestedPairs[i])
                    if NestedConserved == True:
                        ConservedInOutGroups = True
                if ConservedInOutGroups == False:
                    Young.append(pair)
    return Old, Young


# use this function to translate a coding sequence into a protein
def TranslateCDS(CDS):
    '''
    (str) -> str
    Translate a coding sequence into a protein sequence according to the standard
    genetic code.
    Precondition: the cosing sequence needs to be in frame (ie, frame 1)
    >>> TranslateCDS('ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA')
    MAMAPRTEINSTRING*
    >>> cds_translate('ATGTACTAA')
    MY*
    '''
    # use standard genetic code
    genetic_code = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
                   'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
                   'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
                   'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
                   'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
                   'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                   'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                   'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                   'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
                   'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                   'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                   'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                   'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
                   'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                   'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                   'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
    # convert CDS sequence to upper cases to map codon to genetic code dict
    CDS = CDS.upper()
    # initialize protein sequence, and update with amino acids
    protein = ''
    for i in range(0, len(CDS), 3):
        codon = CDS[i:i+3]
        if codon not in genetic_code:
            protein += 'X'
        else:
            protein += genetic_code[codon]
    return protein


# use this function to reverse complement a DNA sequence
def ReverseComplement(dna):
    '''
    (str) -> (str)
    Return the reverse complementary sequence of string dna
    >>> ReverseComplement('atcg')
    'cgat'
    '''
    # use only valid nucleotides, all other bases are converted to Ns
    valid_bases = {'A', 'T', 'C', 'G'}
    # convert sequence to upper cases
    DNA = dna.upper()
    dna_comp = ''
    for i in DNA:
        if i == 'A':
            dna_comp += 'T'
        elif i == 'T':
            dna_comp += 'A'
        elif i == 'C':
            dna_comp += 'G'
        elif i == 'G':
            dna_comp += 'C'
        elif i not in valid_bases:
            dna_comp += 'N'
    reverse_comp_dna = ''
    for i in reversed(dna_comp):
        reverse_comp_dna += i
    if dna.islower():
        reverse_comp_dna = reverse_comp_dna.lower()
    return reverse_comp_dna


# use this function to convert the CDS file to a dictionary
def ConvertCDSToFasta(CDSFile):
    '''
    (file) -> dict
    Take a file with fasta sequences and return a dictionnary with
    gene name as key and single string sequence as value
    '''
    CDS = {}
    infile = open(CDSFile, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                name = line[1:line.index('.')]
                assert name not in CDS
                CDS[name] = ""
            else:
                CDS[name] += line.upper()
    infile.close
    return CDS

# use this function to convert a fasta file to a dictionary
def ConvertFasta(fasta):
    '''
    (file) -> dict
    Take a file with fasta sequences and return a dictionnary with gene: sequence
    '''
    # convert nematode genomes into a dictionnary
    genome = {}
    infile = open(fasta, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                genome[line[1:]] = ""
                seq_name = line[1:]
            else:
                genome[seq_name] += line
    infile.close
    return genome

    
# use this function to filter out genes with weird sequences
def FilerOutCDSSequences(CodingSeq):
    '''
    (dict) -> dict
    Take a dictionary with gene: CDS pairs and return modified dictionary
    in which sequences without start codon, not mulctiple of 3 and with extra 
    stop codons are removed
    '''
    
    # create a set of genes to remove
    to_remove = set()
    for gene in CodingSeq:
        protein = TranslateCDS(CodingSeq[gene])
        # remove genes with more than 1 stop codon
        if protein.count('*') > 1:
            to_remove.add(gene)
        # remove genes with internal stop codons
        if protein.count('*') == 1 and protein[-1] != '*':
            to_remove.add(gene)
        # remove genes without start codons
        if CodingSeq[gene][:3] != 'ATG':
            to_remove.add(gene)
        # remove sequences that are not multiple of 3
        if len(CodingSeq[gene]) % 3 != 0:
            to_remove.add(gene)
    for gene in to_remove:
        del CodingSeq[gene]
    return CodingSeq


# use this function to remove the terminal stop codon
def RemoveTerminalStop(CodingSeq):
    '''
    (dict) -> dict
    Take the dictionary of gene: CDS pairs and return a modified dictionary 
    in which the terminal stop codons have been removed from the coding sequences
    '''
    for gene in CodingSeq:
        if CodingSeq[gene][-3:].upper() in ['TAA', 'TAG', 'TGA']:
            # remove terminal stop
            CodingSeq[gene] = CodingSeq[gene][:-3]
    return CodingSeq
                   

# use this function to perform a permutation test
def PermutationResampling(Sample1, Sample2, Iteration, statistic = np.mean):
    '''
    (list, list, int, function) -> float    
    Take lists of values from 2 samples, the number of resampling iteration,
    and the statistic to compare (by default the mean of the 2 samples) and
    returns the p-value that statistic for sample1 is different from statistc
    for sample2
    '''
    # compute the mean difference between the 2 samples
    observed_diff = abs(statistic(Sample1) - statistic(Sample2))
    # get the sample size of sample 1 
    Nsample1 = len(Sample1)
    # combined the values from the 2 samples into a single array
    combined = np.concatenate([Sample1, Sample2])
    # create a list to store the mean differences from each iteration    
    diffs = []
    # compute the mean difference Iteration times
    for i in range(Iteration):
        # shuffle the values from the 2 samples
        xs = np.random.permutation(combined)
        # compute the mean difference between 2 samples of same size as samples 1 and 2
        diff = np.mean(xs[:Nsample1]) - np.mean(xs[Nsample1:])
        diffs.append(diff)
    # compute the 2 tails p-value 
    pval = (np.sum(diffs > observed_diff) + np.sum(diffs < -observed_diff)) / Iteration
    return pval
    
    
# use this function to compute the p-distance between 2 aligned sequences
def ProteinDistance(seq1, seq2):
    '''
    (str, str) -> float
    Return the p-distance between protein sequences seq1 and seq2: the fraction
    of mismatched amino acids between seq1 and seq2.
    Precondition: seq1 and seq2 must be aligned and have the same length
    '''

    # test prerequisite of equal length
    assert len(seq1) == len(seq2), 'sequences have different lengths'
    # accept - or ~ as gaps and make sure the sequences have same case
    SEQ1 = seq1.replace('~', '-').upper()
    SEQ2 = seq2.replace('~', '-').upper()
    # set up counters. count mismatches and count positions with gaps or non-valid AAs 
    D = 0
    undefined = 0
    # make a set of valid AAs
    valid_AA = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K',
                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}    
    
    # count differences between seq1 and seq2
    # gapped positions are excluded, undefined positions are excluded
    for i in range(len(SEQ1)):
        # count non-valid positions
        if SEQ1[i] not in valid_AA or SEQ2[i] not in valid_AA:
            undefined += 1
        # count mismatches
        elif SEQ1[i] in valid_AA and SEQ2[i] in valid_AA:
            if SEQ1[i] != SEQ2[i]:
                D += 1
    # distance = N differences / length of sequence (without undefined sites)    
    L = len(SEQ1) - undefined
    if L != 0:
        p = D/ L
        return round(p, 6)
    else:
        return 'NA'
        

# use this function to normalize expression counts in UQ-FPKM
def UpperQuartileFPKM(C, UQ, L):
    '''
    (int, numb, int) -> float
    Take the fragment counts, the upper quartile of fragment count, the gene length
    and return a FPKM estimate of gene expression normalized bu upper quartile
    '''
    
    FPKM_UP = (C * (10**9)) / (UQ * L)
    return FPKM_UP


        
# use this function to get upper-quartile normalized FPKM
def ExpressionNormalization(ModeEncodeExpressionFile):
    '''
    (file) -> dict
    Take a file with mouse expression from ModEncode and returns a dictionary
    with normalized upp-quartile FPKM expression extimate
    '''
    
    # create a dictionary {gene_ID: [count, length]}
    expression = {}
    infile = open(ModeEncodeExpressionFile)
    infile.readline()
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            gene, length, count  = line[0], float(line[2]), float(line[4]) 
            if 'ENS' in gene:
                assert gene not in expression
                expression[gene] = [count, length]
    infile.close()
    
    # create a list with counts, do not include 0 for computation of UQ
    Counts = [expression[gene][0] for gene in expression if expression[gene][0] != 0]
    # compute upper-quartile
    UQ = np.percentile(Counts, 75)
    # create a dict with normalized fpkm
    fpkm = {}
    for gene in expression:
        C, L = expression[gene][0], expression[gene][1]
        FPKM_UQ = UpperQuartileFPKM(C, UQ, L)
        fpkm[gene] = FPKM_UQ
    return fpkm
    

# use this function to compute the expression level of a collection of samples
def ComputeTissueExpression(Folder):
    '''
    (folder) -> list
    Return a list of dictionaries with expression level for a given tissue
    Precondition: all samples of a given tissue are organized in a same folder
    '''
    
    # create a list of files in folder
    files = [i for i in os.listdir(Folder) if '.tsv' in i]
    # create a list to store the dictionaries with normalized expression
    tissue_expression = []
    for i in range(len(files)):
        fpkm = ExpressionNormalization(Folder + '/' + files[i])
        tissue_expression.append(fpkm)
    return tissue_expression


# use this expression to merge expression samples
def MergeSamples(L):
    '''
    (list) -> dict    
    Take a list of dictionaries with gene expression for given samples
    and return a dictionary with gene: list of expression values pairs
    '''
    
    # create a dictionary to store the expression values
    expression = {}
    for i in range(len(L)):
        for gene in L[i]:
            if gene in expression:
                expression[gene].append(L[i][gene])
            else:
                expression[gene] = [L[i][gene]]
    return expression
    

# use this function to generate gene expression profile for the same tissues in human and mouse
def MatchHumanToMouseExpressionProfiles(HumanExpression):
    '''
    (dict) -> dict
    Take the dictionary of human gene expression profiles from GTEX and return
    a dictionary with expression profiles matching the tissues of gene expression
    profiles in mouse
    Precondition: tissues for mouse expression have already been filtered to 
    match human tissues
    '''
    # Human tissues from GTEX include:
    HumanTissues = ['Adipose Tissue', 'Adrenal Gland', 'Bladder', 'Blood',
                    'Blood Vessel', 'Brain', 'Breast', 'Cervix Uteri', 'Colon',
                    'Esophagus', 'Fallopian Tube', 'Heart', 'Kidney', 'Liver',
                    'Lung', 'Muscle', 'Nerve', 'Ovary', 'Pancreas', 'Pituitary',
                    'Prostate', 'Salivary Gland', 'Skin', 'Small Intestine', 'Spleen',
                    'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina']
    # mouse tissues from ModEncode include:
    MouseTissues = ['Adipose Tissue',  'Adrenal Gland', 'Bladder', 'Brain',
                    'Breast', 'Colon', 'Heart', 'Kidney', 'Liver', 'Lung',
                    'Ovary', 'Pancreas', 'Small Intestine', 'Spleen',
                    'Stomach', 'Testis']
    # make a list of indices for mouse tissues in human mouse tissues
    # to extract expression of human tissues in these tissues
    ToExtract = [HumanTissues.index(tissue) for tissue in MouseTissues]

    # create a new dictionary for human expression profiles
    Profile = {}
    # loop over human genes
    for gene in HumanExpression:
        # initialize gene: list pair
        Profile[gene] = []
        # loop over indices
        for i in ToExtract:
            # add expression level at index position
            Profile[gene].append(HumanExpression[gene][i])
    return Profile
    
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
        if 'gene' in line and not line.startswith('#'):
            line = line.rstrip().split('\t')
            # consider protein-coding genes            
            if line[2] == 'gene':
                # check that gene is protein-coding            
                biotype = line[8][line[8].index('biotype=')+8: line[8].index(';', line[8].index('biotype=')+8)]
                if biotype == 'protein_coding':
                    # get the gene ID
                    gene = line[8][line[8].index('ID=gene:')+8: line[8].index(';')]
                    # get gene name
                    name = line[8][line[8].index('Name=') + len('Name='): line[8].index(';', line[8].index(';')+1)]
                    Names[gene] = name
    infile.close()
    return Names    


# use this function to parse the Cosmic list of cancer genes
def ParseCosmicFile(Cosmic, GeneNames):
    '''
    (file, dict) -> set
    Take the file of cancer gene census from the Cosmic database, the dictionary
    of gene name: gene ID matches and return a set of gene for which mutations
    have been causally implicated in cancer
    '''
    
    CancerGenes = set()
    infile = open(Cosmic)
    infile.readline()
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            if line[0] in GeneNames:
                CancerGenes.add(GeneNames[line[0]])
    infile.close()
    return CancerGenes
 
# use this function to parse the file of Complex disease genes
def ParseComplexDisease(ComplexDisease, GeneNames):
    '''
    (file, dict) -> set    
    Take the file with associations between genes and complex diseases the 
    dictionary of gene name: gene ID matches and return a set of genes associated
    with diseases
    '''
    GAD = set()
    infile = open(ComplexDisease)
    infile.readline().rstrip().split('\t')
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            # get gene name
            gene = line[5]
            if '"' in gene:
                assert gene.count('"') == 2 and gene[0] == '"' and gene[-1] == '"'
                gene = gene[1:-1]
            # remove space in gene
            gene = gene.replace(' ', '')
            # check if multiple genes are listed
            if ',' in gene:
                assert ':' not in gene and ';' not in gene
                gene = gene.split(',')
            elif ';' in gene:
                assert ':' not in gene and ',' not in gene
                gene = gene.split(';')
            elif ':' in gene:
                assert ';' not in gene and ',' not in gene
                gene = gene.split(':')
            # populate set with gene names       
            if type(gene) == str:
                if gene in GeneNames:
                    GAD.add(GeneNames[gene])
            elif type(gene) == list:
                for item in gene:
                    if item in GeneNames:
                        GAD.add(GeneNames[item])
            # get alternative gene names
            alternative = line[1].split('|')
            while '"' in alternative:
                alternative.remove('"')
            while '' in alternative:
                alternative.remove('')
            for name in alternative:
                if name in GeneNames:
                    GAD.add(GeneNames[name])
    infile.close()
    return GAD

# use this function to parse the GWS file
def ParseGWASDisease(GWASFile, TraitsToRemove, GeneNames):
    '''
    (file, file, dict) -> set
    Take the file of GWAS associations between traits and genes, a file with
    non-disease traits to remove, a dict with gene name: gene ID matches and 
    return a set of genes associated with diseases
    '''
    # make a set of GWAS disease genes
    GWAS = set()
    # exclude non-disease traits
    ExcludeTraits = set()
    infile = open(TraitsToRemove)
    for line in infile:
        if line.rstrip() != '':
            line = line.strip()
            ExcludeTraits.add(line)
    infile.close()
    # open GWAS catalog
    infile = open(GWASFile, encoding='utf8')
    infile.readline()
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            # check that trait is disease only        
            if line[7] not in ExcludeTraits:
                genes = line[13].replace(' ', '')
                # check if multiple genes are listed
                if ', ' in genes:
                    genes = genes.split(',')
                if type(genes) == str:
                    if genes in GeneNames:
                        GWAS.add(GeneNames[genes])
                elif type(genes) == list:
                    for item in genes:
                        if item in GeneNames:
                            GWAS.add(GeneNames[item])       
    infile.close()
    return GWAS


# use this function to parse the file of tumor driver genes
def ParseTumorDrivers(DriverGenes):
    '''
    (file) -> set
    Prse the file of cancer drivers to return a set of driver genes
    '''
    # make a set of cancer driver genes
    Drivers = set()
    infile = open(DriverGenes)
    infile.readline()
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            Drivers.add(line[2][:line[2].index('.')])
    infile.close()
    return Drivers


# use this function to parse the OMIM datafile
def ParseOMIMDisease(TitleFile, OMIMFile, GeneNames):
    '''
    (file, file, dict) -> set
    Take the file with phenotypes, the file with gene: phenotype associations,
    a dictionary with gene name: gene ID matches and return a set of meadlean
    disease genes
    '''
    
    # make a set of mendelian disease genes
    OMIM = set()
    # get the mim IDs corresponding to associations between phenotypes and genes
    mimIDs = set()
    infile = open('mimTitles.txt')
    for line in infile:
        if (not line.startswith('#')) and line.rstrip() != '':
            line = line.rstrip().split('\t')
            if line[0] == 'Number Sign' or line[0] == 'Percent' or line[0] == 'Plus':
                mimIDs.add(line[1])
    infile.close()
    # get the set of phenotype associated genes
    infile = open('morbidmap.txt')
    for line in infile:
        if (not line.startswith('#')) and line.rstrip() != '':
            line = line.rstrip().split('\t')
            pheno = line[0].replace(' ', '')
            pheno = pheno.split(',')
            pheno = pheno[-1]
            pheno = pheno[:pheno.index('(')]
            if pheno in mimIDs:
                # extract the genes
                genes = line[-3].replace(' ', '')
                # check if multiple genes are listed
                if ',' in genes:
                    genes = genes.split(',')
                if type(genes) == str:
                    if genes in GeneNames:
                        OMIM.add(GeneNames[genes])
                elif type(genes) == list:
                    for item in genes:
                        if item in GeneNames:
                            OMIM.add(GeneNames[item])
    infile.close()
    return OMIM

# use this function to save a list of overlapping genes as a dict in a json file
def SaveOverlappingPairsToFile(OverlappingPairs, JsonFileName):
    '''
    (list, str) -> file
    Take a list of overlapping gene pairs and save overlapping relationships to a json file 
    '''

    # create a dictionary of overlapping genes
    OverlappingGenes = {}
    # loop over gene pairs in list
    for pair in OverlappingPairs:
        # use first gene in pair as key
        if pair[0] in OverlappingGenes:
            OverlappingGenes[pair[0]].append(pair[1])
        else:
            OverlappingGenes[pair[0]] = [pair[1]]
    newfile = open(JsonFileName + '.json', 'w')
    json.dump(OverlappingGenes, newfile, sort_keys = True, indent = 4)
    newfile.close()
    
    
# use this file to parse the GO annotation file
def ParseGOFile(GOFile):
    '''
    (file) -> dict
    Take the file with GO annotations and return a dictionary of gene name: set of GO IDs
    '''
    # create a dictionary to map gene names to GO terms
    Annotations = {}
    infile = open(GOFile)
    for line in infile:
        if not line.startswith('!'):
            if 'GO:' in line:
                line = line.rstrip().split('\t')
                # get gene and GO term ID
                gene  = line[2]
                GOid = [item for item in line if item.startswith('GO:')]
                # check if gene is already recorded
                if gene not in Annotations:
                    Annotations[gene] = set()
                for GO in GOid:
                    Annotations[gene].add(GO)
    infile.close()
    return Annotations             


# use this function to generate a set of valid GO IDs
def FilterGOTerms(GOClass):
    '''
    (file) -> set
    Take a file of GO IDs corresponding to a specific GIO class (biological processes,
    molecular function) and return a set of IDs for this class
    '''
    GOids = set()
    infile = open(GOClass)
    for line in infile:
        if line.startswith('GO:'):
            line = line.rstrip().split('\t')
            GOids.add(line[0])
    return GOids

# use this function to map Ensembl gene ID to GO annotations
def MapEnsemblGenesToGOTerms(Annotations, GeneNames):
    '''
    (dict, dict) -> dict
    Take the dictionary of gene_names:GO_terms pairs, the dictionary of
    gene_IDs:_gene_names pairs and return a dictionary of gene_IDs:GO_terms
    '''
    
    # create a dictionary to map ensembl IDs to GO terms {gene: {GO1, GO2...}}
    GeneOntology = {}
    # reverse  the diction of gene name {name: ID}    
    Names = {}
    for gene in GeneNames:
        Names[GeneNames[gene]] = gene
    # loop over gene names with GO annotations
    for name in Annotations:
        # map name to gene
        if name in Names:
            gene = Names[name]
            # get assign the GO terms to gene
            GeneOntology[gene] = Annotations[name]
    return GeneOntology

    
# use this function to compute jacaard similarity index
def JaccardIndex(A, B):
    '''
    (set, set) -> num
    Take 2 sets of values and return the Jaccard index, mesuring the similarity
    between the 2 sets
    '''
    # by default Jaccard index  = 1 when both sets are empty
    if len(A) == 0 and len(B) == 0:
        return 1
    else:
        return len(A.intersection(B)) / len(A.union(B))    
    

# use this function to add significance of pairwise comparisons to a graph
def AddSignificanceToBars(ax, SignificanceLevel, XLine1, XLine2, YLine, XText, YText):
    '''
    (ax, str, num, num, num, num, num) -> ax
    Take a matplotlib ax object, the significance level (as stars), the positions
    of the bracket and star and return the ax with annotated significance level
    '''
    ax.annotate("", xy=(XLine1, YLine), xycoords='data', xytext=(XLine2, YLine), textcoords='data',
                 arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 0.7))
    # add stars for significance
    ax.text(XText, YText, SignificanceLevel, horizontalalignment='center', verticalalignment='center',
            color = 'grey', fontname = 'Arial', size = 7)
    return ax   
    
    
# use this function to convert p values to star significance
def ConvertPToStars(PValues):
    '''
    (list) -> list
    Take a list of P-values and return a list with star strings showing significance    
    '''
    Significance = []
    for pvalue in PValues:
        if pvalue >= 0.05:
            Significance.append('')
        elif pvalue < 0.05 and pvalue >= 0.01:
            Significance.append('*')
        elif pvalue < 0.01 and pvalue >= 0.001:
            Significance.append('**')
        elif pvalue < 0.001:
            Significance.append('***')
    return Significance
    
    
# use this function to set values > cutoff equal to cutoff    
def CombineHighValues(L, cutoff):
    '''
    (list, int) -> list
    Take a list of overlap length and return a modified list with values higher
    than cutoff equal to cutoff
    '''
    for i in range(len(L)):
        if L[i] >= cutoff:
            L[i] = cutoff
    return L    
    

# use this function to count the number of disease and non-disease genes for 
# each disease type and each gene category
def CountDiseaseGenes(GeneList, DiseaseList):
    '''
    (list, list) -> list    
    Take a list of gene sets and a list of sets of disease genes and return
    a list with lists of disease and non-disease gene counts for each category
    for each disease class
    '''
    
    GeneCounts = []
    for DiseaseType in DiseaseList:
        counts = []
        for i in range(len(GeneList)):
            disease = len([gene for gene in GeneList[i] if gene in DiseaseType])
            nondisease = len([gene for gene in GeneList[i] if gene not in DiseaseType])
            counts.append([disease, nondisease])
        GeneCounts.append(counts)
    return GeneCounts
