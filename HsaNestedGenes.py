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
  
  
  
# use this function to match human genes with orrthologs in a single other species
def MatchOrthologPairs(OrthoFile):
    '''
    (file) -> dict
    Take a file with orthology assignment between 2 species and return a dictionary
    of 1:1 orthologs between the 2 species. 
    '''
    # the file colums foloww the format:
    # Sp1GeneID Sp1TranscriptID Sp2GeneID Sp2GeneName HomologyType OrthologyC
    
    # create a dict of orthologs
    Orthos = {}
    infile = open(OrthoFile)
    infile.readline()
    for line in infile:
        # consider only 1:1 orthologs
        if 'ortholog_one2one' in line:
            line = line.rstrip().split('\t')
            # get gene IDs of the 2 species
            gene1, gene2 = line[0], line[2]
            # check that genes are ensembl gene IDs
            for i in [gene1, gene2]:
                assert 'ENS' in i, 'gene id is not valid'
                assert 'ortholog' in line[4], 'ortholog should be in homology type'
            assert gene1 not in Orthos, 'gene is already matched to a 1:1 ortholog'
            Orthos[gene1] = gene2
    infile.close()                      
    # check that all orthologs are 1;1 orthologs
    values = list(Orthos.values())
    for i in values:
        assert values.count(i) == 1
    return Orthos
    
  
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


# use this function to parse the GTEX median expression file
def ParseGTEXExpression(GTEXExpressionFile):
    '''
    (file) -> dict
    Take the file from GTEX with median expression in each tissue and return a
    dictionary with gene: list of expression in tissues key, value pairs
    '''
    # create a dict to store the expression (median tissue expression FPKM normalized by upper quartiles))
    Expression = {}
    infile = open(GTEXExpressionFile)
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


# use this function to generate un-nested pairs of genes to randomly draw
def GenerateUnNestedGenePairs(HostGenes, GeneCoord, OrderedGenes, ExpressionProfile):
    '''
    (dict, dict, dict, dict) -> dict
    Take the dictionary of host and nested genes, the dictionary of gene coordinates,
    the dictionary of ordered genes along each chromosome, the dictionary of
    gene: expression pairs and return a dictionary of pairs of number: gene pair on each chromsome
    that are distant of 500 bp at most
    '''
    
    # make a set of nested and host genes
    IncludingNestedGenes = set()
    for gene in HostGenes:
        IncludingNestedGenes.add(gene)
        for nested in HostGenes[gene]:
            IncludingNestedGenes.add(nested)
    
    # make a dictionary with chromsome as key and number: gene pairs {chromo: {num: [gene1, gene2]}}    
    ToDrawFrom = {}
    # loop over chromosomes
    for chromo in OrderedGenes:
        # set up counter
        k = 0
        # add chromo as key and intialize inner dict
        ToDrawFrom[chromo] = {}
        # loop over the list of ordered genes
        for i in range(len(OrderedGenes[chromo])):
            # check that gene is not host or nested, has expression
            if OrderedGenes[chromo][i] in ExpressionProfile and OrderedGenes[chromo][i] not in IncludingNestedGenes:
                # get the end position of gene 1
                EndGene1 = GeneCoord[OrderedGenes[chromo][i]][2]                
                # grab 2nd gene to form a pair                
                for j in range(i+1, len(OrderedGenes[chromo])):
                    # check that gene is not host or nested and has expression
                    if OrderedGenes[chromo][j] in ExpressionProfile and OrderedGenes[chromo][j] not in IncludingNestedGenes:
                        # get the start position of gene 2
                        StartGene2 = GeneCoord[OrderedGenes[chromo][j]][1]
                        # check if distance is less that 500 bp
                        D = StartGene2 - EndGene1
                        if D >= 0 and D <= 2000:
                            # add gene pair and update counter
                            ToDrawFrom[chromo][k] = [OrderedGenes[chromo][i], OrderedGenes[chromo][j]]
                            k += 1
    return ToDrawFrom



# use this function to generate un-nested genes to randomly draw
def GenerateAllUnNestedGenes(NoOverlap, OrderedGenes):
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
            # check that gene does not overlap with any other gene
            if OrderedGenes[chromo][i] in NoOverlap:
                # add gene pair and update counter
                ToDrawGenesFrom[chromo][k] = OrderedGenes[chromo][i]
                k += 1
    return ToDrawGenesFrom


# use this function to sort young and ancestral nesting events
def InferYoungOldNestingEvents(FirstSpOrthologs, SecondSpOverlapPairs, OutGroupOverlapPairs, FirstSpHostNestedPairs):
    '''
    (dict, list, list, list) -> (list, list)
    Take a dictionary with ortholog gene trios, the lists of overlapping gene
    pairs in the sister species and in the outgroup and the list of host:nested
    pairs in the focal species and return a tuple with list of host: nested pairs
    that are infered to be old and young (before the divergence of the 2 species or after)
    '''   
    # create lists of sets of gene pairs to remove the order between genes
    SecondOverlap, OutGroupOverlap = [], []
    for pair in SecondSpOverlapPairs:
        SecondOverlap.append(set(pair))
    for pair in OutGroupOverlapPairs:
        OutGroupOverlap.append(set(pair))
    # create lists of host-nested gene pairs that are old (present in second
    # species, or second species and outgroup) or young (not found in outgroup
    # and not found in the second species, ie species-specific)    
    old, young = [], []
    # loop over the gene pairs    
    for pair in FirstSpHostNestedPairs:
        # check that both genes have orthologs
        if pair[0] in FirstSpOrthologs and pair[1] in FirstSpOrthologs:
            # check if gene pair is overlapping in second species and outgroup
            # get ortholog of host in 2nd species
            hostorthosp2 = FirstSpOrthologs[pair[0]][0]
            # get ortholog of nested in 2nd species
            nestedorthosp2 = FirstSpOrthologs[pair[1]][0]
            # get ortholog of host in outgroup
            hostorthoout = FirstSpOrthologs[pair[0]][1]
            # get ortholog of nested gene in outgroup
            nestedorthoout = FirstSpOrthologs[pair[1]][1]
            # check if orthologs are overlapping in 2nd species
            if {hostorthosp2, nestedorthosp2} not in SecondOverlap and {hostorthoout, nestedorthoout} not in OutGroupOverlap:
                # nested is species-specific
                young.append(pair)
            elif {hostorthosp2, nestedorthosp2} in SecondOverlap and {hostorthoout, nestedorthoout} not in OutGroupOverlap:
                # gene pair is overlapping in 2nd species, suffcient to be called old nested pair
                old.append(pair)
            elif {hostorthosp2, nestedorthosp2} in SecondOverlap and {hostorthoout, nestedorthoout} in OutGroupOverlap:
                # gene pair is overlapping in 2nd species and outgroup, infer old nested pair
                old.append(pair)

    return old, young


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
                   


