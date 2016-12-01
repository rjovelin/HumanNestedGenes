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
        alignment = ConvertFasta('./HumanChimpPairs/' + filename)
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

# generate codeml control files
# get the list of alignment files
files = os.listdir('./HumanChimpCodeml/')
for filename in files:
    if filename[-4:] == '.txt' and 'ENS' in filename:
        # get the file name without extension including gene names
        name = filename[:-4]
        # open control file for reading
        ctlfile = open('codeml.ctl', 'r')
        # open outputfile file for writing
        newfile = open('./HumanChimpCodeml/' + name + '_codeml.ctl', 'w')
        # go through the control file, copy to new file with new input file and output file
        for line in ctlfile:
            if 'seqfile' not in line and 'outfile' not in line:
                newfile.write(line)
            elif 'seqfile' in line:
                line = line.rstrip().split()
                # replace the file name with the alignment file name
                line[2] = filename
                newfile.write(' '.join(line) + '\n')
            elif 'outfile' in line:
                line = line.rstrip().split()
                # replace the outputfile name with alignment_file.out
                line[2] = name + '.out.txt'
                newfile.write(' '.join(line) + '\n')
        ctlfile.close()
        newfile.close()
print('done generating codeml control files')

# copy necessary files to directory with input and control files
# copy files to working directory
os.system('cp 2NG.dN 2NG.dS 2NG.t 4fold.nuc lnf HumanChimpTree.tre.txt rst rst1 rub ./HumanChimpCodeml/')

# go to directory containing control files and input files
os.chdir('./HumanChimpCodeml/')

# run codeml
# get the the list of control files
ctlfiles = [i for i in os.listdir('./') if '_codeml.ctl' in i and 'ENS' in i]
for filename in ctlfiles:
    os.system('codeml ' + filename)
print('done running codeml')

# write summary of nucleotide divergence to file 

# make a list of codeml output files in current directory
outfiles = [i for i in os.listdir('./') if '.out.txt' in i]

# save summary file in parent directory
newfile = open('../HumanChimpSeqDiverg.txt', 'w')
# write header
newfile.write('\t'.join(['Human_Gene', 'Chimp_Gene', 'dN', 'dS', 'dN/dS']) + '\n')

# loop over outputfiles
for filename in outfiles:
    infile = open(filename, 'r')
    # get human gene name
    HumanGene = filename[:filename.index('_')]
    # get chimp gene name
    ChimpGene = filename[filename.index('_')+1: filename.index('.out')]
    # loop through file    
    for line in infile:
        if 'tree length for dN' in line:
            line = line.split()
            dN = float(line[-1])
        elif 'tree length for dS' in line:
            line = line.split()
            dS = float(line[-1])
    if dS != 0:
        omega = dN / dS
    elif dS == 0:
        omega = 'NA'
    # write results to file    
    newfile.write('\t'.join([HumanGene, ChimpGene, str(dN), str(dS), str(omega)]) + '\n')
    # close outputfile
    infile.close()
# close summary file
newfile.close()
print('done writing summary file')

# go back to parent directory
os.chdir('../')

