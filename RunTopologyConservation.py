# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:40:43 2017

@author: RJovelin
"""

import os
import sys

species = sys.argv[1]


# open file for writing command
outputfile = 'CountOrthologsConservedTopology.sh' 
newfile = open(outputfile, 'w')
newfile.write('cd /RQusagers/rjovelin/richard/Richard/Nested_Genes/')
newfile.write('\n')
newfile.write('python3 PlotTopologyConservationDistance.py ' + species)
newfile.close()            

MyCommand = 'qsub -l walltime=10:00:00,nodes=1:ppn=4 ' + outputfile
os.system(MyCommand)
print(MyCommand)
        
