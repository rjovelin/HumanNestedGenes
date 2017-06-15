# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:40:43 2017

@author: RJovelin
"""

import os

# open file for writing command
outputfile = 'PerformNeutralSimulations.sh' 
newfile = open(outputfile, 'w')
newfile.write('cd /RQusagers/rjovelin/richard/Richard/Nested_Genes/')
newfile.write('\n')
newfile.write('python3 SimulationsNeutralOverlap.py')
newfile.close()            

MyCommand = 'qsub -l walltime=20:00:00,nodes=1:ppn=12,mem=46g ' + outputfile
os.system(MyCommand)
print(MyCommand)
        
