# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 21:04:34 2016

@author: Richard
"""

import numpy as np
import math


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
    (file, str) -> dict
    Return a dictionnary with transcript as key and exon positions list pairs
    as value. Note that exons define intron positions but are not necessarily
    entirely coding    
    '''
    
    # open file for reading
    infile = open(gff_file, 'r')
    # make a dictionnary with chromosome as key and a dictionnary as value
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
  
# use this function to get the orientation of a pair of genes
def GenePairOrientation(GenePair, GeneCoord):
    '''
    (list, dict) -> list
    Take a pair of gene, a dictionary with the genes coordinates and return a 
    list with the strand orientation of the 2 genes
    '''
    return [GeneCoord[GenePair[0]][-1], GeneCoord[GenePair[1]][-1]] 
  
  
# Map humangenes to their orthologs in 2 other species  
def ParseOrthologFile(OrthoFile):
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
                    assert len(line) == 8
                    # get gene IDs of the 3 species
                    gene1, gene2, gene3 = line[0], line[2], line[5]
                    # check that genes are ensembl gene IDs
                    for i in [gene1, gene2, gene3]:
                        assert 'ENS' in i, 'gene id is not valid'
                    assert 'ortholog' in line[3] and 'ortholog' in line[6], 'ortholog should be in homology type'
                    # record 1:1 orthologs
                    if line[3] == 'ortholog_one2one' and line[6] == 'ortholog_one2one':
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
    across individuakls for each of the tissues
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
                kdexp, lvexp, tsexp = list(map(ExpVal, line[33: 35])), list(map(ExpVal, line[35: 37])), list(map(ExpVal, line[37]))
            elif species == 'macaque':
                gene = line[4]
                brexp, cbexp, htexp = list(map(ExpVal, line[70: 73])), list(map(ExpVal, line[73: 75])), list(map(ExpVal, line[75: 77]))
                kdexp, lvexp, tsexp = list(map(ExpVal, line[77: 79])), list(map(ExpVal, line[79: 81])), list(map(ExpVal, line[81:]))
            # get the median expression level per tissue
            assert gene not in expression, 'gene is already recorded'
            expression[gene] = [np.median(brexp), np.median(cbexp), np.median(htexp), np.median(kdexp), np.median(lvexp), np.median(tsexp)]
    infile.close()
    return expression


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



# use this function to make a set of host and nested genes
def MakeHostNestedGeneSet(HostGenes):
    '''
    (dict) -> set
    Take the dictionary of host: nested genes (or contained) and return a set
    of host and nested genes (or contained) 
    '''
    # HostGenes is in the form {gene: [transcripts]}
    hostnested = set()
    for i in HostGenes:
        hostnested.add(i)
        for j in HostGenes[i]:
            hostnested.add(j)
    return hostnested


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
def FilterGenePairsWithoutExpression(HostNestedPairs, ExpressionProfile):
    '''
    (list, dict) -> list
    Take the list of gene pairs and the dictionary of expression profile for
    each gene and return a modified list of gene pairs in whhich pairs of gene
    lacking expression are removed
    Precondition: genes without expression in the expression profile dictionary
    have been removed    
    '''
    
    # filter gene pairs based on expression
    to_remove = []
    for pair in HostNestedPairs:
        # remove gene pair if at least one of the gene is not expressed
        if pair[0] not in ExpressionProfile or pair[1] not in ExpressionProfile:
            to_remove.append(pair)
    for pair in to_remove:
        HostNestedPairs.remove(pair)
    return HostNestedPairs


