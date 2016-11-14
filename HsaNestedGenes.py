# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 21:04:34 2016

@author: Richard
"""


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
def FindContainedGenePairs(GeneCoordinates, Overlap):
    '''
    (dict, dict) -> dict
    Take the dictionary with gene coordinate, the dictionary with overlapping
    gene pairs and return a dictionary with gene pairs in which one gene is
    fully contained in another gene
    '''
 
    # GeneCoordinates is in the form {gene:[chromo, start, end, sense]}
    # Overlap is in the form {gene1: [gene2, gene3]}
        
    # create a dict with gene containing other genes {containing: [contained1, contained2]}
    ContainingGenes = {}
    # loop over each overlapping gene
    for gene1 in Overlap:
        # check if each gene in the overlapping gene is fully located within the over gene 
        for gene2 in Overlap[gene1]:
            # check if one of the 2 genes is fully contained in the other gene
            coord1 = set(range(GeneCoordinates[gene1][1], GeneCoordinates[gene1][2])) 
            coord2 = set(range(GeneCoordinates[gene2][1], GeneCoordinates[gene2][2]))
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
    