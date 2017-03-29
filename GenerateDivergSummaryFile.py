# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:43:46 2016

@author: RJovelin
"""

# use this script to write summary of nucleotide divergence to file 

import os
import sys

# usage GenerateDivergSummaryFile.py [options]
# -['mouse'/'chimp']: use human-mouse or human-chimp divergence

# get option from command
Species = sys.argv[1]
assert Species in ['chimp', 'mouse']

# get folder with divergence files
if Species == 'chimp':
    Folder = './HumanChimpCodeml/'
elif Species == 'mouse':
    Folder = './HumanMouseCodeml/'


# make a list of codeml output files in current directory
outfiles = [i for i in os.listdir(Folder) if '.out.txt' in i]

# save summary file in parent directory
if Species == 'chimp':
    SummaryFile = './HumanChimpSeqDiverg.txt'
    header = '\t'.join(['Human_Gene', 'Chimp_Gene', 'dN', 'dS', 'dN/dS'])
elif Species == 'mouse':
    SummaryFile = './HumanMouseSeqDiverg.txt'
    header = '\t'.join(['Human_Gene', 'Mouse_Gene', 'dN', 'dS', 'dN/dS'])

newfile = open(SummaryFile, 'w')
# write header
newfile.write(header + '\n')



#### continue here



# loop over outputfiles
for filename in outfiles:
    infile = open('./HumanChimpCodeml/' + filename, 'r')
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
    
    # filter out genes with large dN or dS or omega values
    if omega == 'NA':
        if dN <= 1 and dS <= 1:
            # write results to file    
            newfile.write('\t'.join([HumanGene, ChimpGene, str(dN), str(dS), str(omega)]) + '\n')
    else:
        if dN <= 1 and dS <= 1 and omega <= 10:
            newfile.write('\t'.join([HumanGene, ChimpGene, str(dN), str(dS), str(omega)]) + '\n')
            
    # close outputfile
    infile.close()
# close summary file
newfile.close()
print('done writing summary file')

