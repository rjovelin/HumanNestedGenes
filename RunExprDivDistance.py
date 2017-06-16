# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:24:39 2017

@author: RJovelin
"""

import os

# open file for writing command
outputfile = 'RunExprDivDistance.sh' 
newfile = open(outputfile, 'w')
newfile.write('cd /RQusagers/rjovelin/richard/Richard/Nested_Genes/')
newfile.write('\n')
newfile.write('python3 PlotExpressionDivergenceDistance.py')
newfile.close()            

MyCommand = 'qsub -l walltime=20:00:00,nodes=1:ppn=12,mem=46g ' + outputfile
os.system(MyCommand)
print(MyCommand)
        
