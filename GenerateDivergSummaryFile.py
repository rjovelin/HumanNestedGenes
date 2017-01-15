# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:43:46 2016

@author: RJovelin
"""

# use this script to write summary of nucleotide divergence to file 

import os

# make a list of codeml output files in current directory
outfiles = [i for i in os.listdir('./HumanChimpCodeml/') if '.out.txt' in i]

# save summary file in parent directory
newfile = open('./HumanChimpSeqDiverg.txt', 'w')
# write header
newfile.write('\t'.join(['Human_Gene', 'Chimp_Gene', 'dN', 'dS', 'dN/dS']) + '\n')

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
    if dN <= 1 and dS <= 1:
            newfile.write('\t'.join([HumanGene, ChimpGene, str(dN), str(dS), str(omega)]) + '\n')        
        
     
#    if omega == 'NA':
#        if dN <= 2 and dS <= 2:
#            # write results to file    
#            newfile.write('\t'.join([HumanGene, ChimpGene, str(dN), str(dS), str(omega)]) + '\n')
#    else:
#        if dN <= 2 and dS <= 2 and omega <= 2:
#            newfile.write('\t'.join([HumanGene, ChimpGene, str(dN), str(dS), str(omega)]) + '\n')
            
            
            
    # close outputfile
    infile.close()
# close summary file
newfile.close()
print('done writing summary file')

