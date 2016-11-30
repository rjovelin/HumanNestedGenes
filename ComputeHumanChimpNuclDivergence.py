# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 23:30:44 2016

@author: Richard
"""


# use this script to compute nucleotide divergence between human-chimp orthologs

# usage python3 ComputeHumanChimpNuclDivergence.py


import os
from HsaNestedGenes import *


# create a directory to 
os.mkdir('HumanChimpCodeml')

# convert alignment files to PAML format

# get a list with the file names
files = os.listdir('./HumanChimpPairs/')
for filename in files:
    # check that file is the t-coffee DNA alignment file
    if filename[-8:] == '_aln.tfa':
        # convert sequences to dictionnary
        alignment = ConvertFasta(filename)
        # get the name of the file without '_aln.tfa'
        name = filename[:-8]
        # make a list with gene names as they appear in the file
        seq_names = []
        infile = open('./HumanChimpPairs/' + filename)
        for line in infile:
            if line.startswith('>'):
                seq_names.append(line.rstrip()[1:])
        infile.close()
        # open newfile for writing
        newfile = open('HumanChimpCodeml/' + name + '.txt', 'w')
        # add the number of sequences and the length of the alignment in the first line
        newfile.write(str(len(alignment)) + ' ' + str(len(alignment[seq_names[0]])) + '\n')
        # add orthologous sequences in fasta format
        newfile.write('>' + seq_names[0] + '\n')
        newfile.write(alignment[seq_names[0]] + '\n')
        newfile.write('>' + seq_names[1] + '\n')
        newfile.write(alignment[seq_names[1]] + '\n')
        newfile.close()
print('done convertir alignments to PAML format')




######################## edit below 




# generate codeml control files
# get the list of alignment files
files = os.listdir('./')
for filename in files:
    if filename[-4:] == '.txt' and 'CRE' in filename and 'tre' not in filename:
        generate_codeml_control_file(filename, 'codeml.ctl', './')
print('done generating codeml control files')

# run codeml
# get the the list of control files
files = os.listdir('./')
for filename in files:
    if '.ctl' in filename and 'CRE' in filename:
        os.system('codeml ' + filename)
print('done running codeml')

# save divergent results to file
save_divergence_results_to_file('../Genome_Files/356_10172014.gff3', '../Genome_Files/534_10172014.gff3', 'CRM_CLA_protein_divergence.txt', './')
print('done computing sequence divergence')









            
# use this function to create a codeml control file for a given alignment file
def generate_codeml_control_file(alignment_file, control_file, destination): 
    '''
    (file, file, str) -> file
    Take a PAML alignment file and template codeml control file and generate
    a new control file, with parameters seqfile and outfile modified to 
    correspond to the alignment file
    '''

    # check that alignment_file is a text file
    if alignment_file[-4:] == '.txt':
        # get the file name
        name = alignment_file[:-4]
        
    # open control file for reading
    ctlfile = open(control_file, 'r')
        
    # open outputfile file for writing
    newfile = open(destination + name + '_codeml.ctl', 'w')
    
    # go through the control file, copy to new file with new input file and output file
    for line in ctlfile:
        if 'seqfile' not in line and 'outfile' not in line:
            newfile.write(line)
        elif 'seqfile' in line:
            line = line.rstrip().split()
            # replace the file name with the alignment file name
            line[2] = alignment_file
            for item in line[:-1]:
                newfile.write(item + ' ')
            newfile.write(line[-1] + '\n')
        elif 'outfile' in line:
            line = line.rstrip().split()
            # replace the outputfile name with alignment_file.out
            line[2] = name + '.out.txt'
            for item in line[:-1]:
                newfile.write(item + ' ')
            newfile.write(line[-1] + '\n')
    ctlfile.close()
    newfile.close()


# use this function to parse the codeml output files and combine results into a single table
def save_divergence_results_to_file(crem_gff, cla_gff, outputfile, source):
    '''
    (file, file, file) -> file
    Extract the dN and dS values from the codml output files, 
    and save a table to outputfile with columns including the gene names
    and transcript names of each sequence using the GFF annotation files
    '''
    
    # make a list of files in current directory
    files = os.listdir(source)

    # open outfile for writing
    newfile = open(source + outputfile, 'w')
    newfile.write('Cremanei_transcript' + '\t'  + 'Cremanei_gene' + '\t' + 'Clatens_transcript' + '\t' + 'Clatens_gene' + '\t' + 'dN' + '\t' + 'dS' + '\t' + 'dN/dS' + '\n')

    # get latens transcript : gene pairs
    cla_transcripts = transcript_to_gene(cla_gff)

    # get remanei transcript : gene pairs
    crem_transcripts = transcript_to_gene(crem_gff)
    
    # loop over the files
    for filename in files:
        if '.out.txt' in filename:
            infile = open(filename, 'r')
            crem_TS = filename[:filename.index('_CLA')]
            cla_TS = filename[filename.index('CLA'): filename.index('.out')]
            for line in infile:
                if 'tree length for dN' in line:
                    line = line.split()
                    dN = float(line[-1])
                elif 'tree length for dS' in line:
                    line = line.split()
                    dS = float(line[-1])
            crem_gene = crem_transcripts[crem_TS]
            cla_gene = cla_transcripts[cla_TS]
            if dS != 0:
                omega = dN / dS
            elif dS == 0:
                omega = 'NA'
            newfile.write(crem_TS + '\t' + crem_gene + '\t' + cla_TS + '\t' + cla_gene + '\t' + str(dN) + '\t' + str(dS) + '\t' + str(omega) + '\n')
            infile.close()

    newfile.close()
                   