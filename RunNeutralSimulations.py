# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:33:35 2018
@author: RJovelin
"""

import os
import sys
import getopt

# use this script to generate shell scripts to run neutral simulations and launch jobs


def Usage():
    print("""
    usage: RunNeutralSimulations.py [-h|--Help|-w|--Walltime|-i|--Iterations]
    -h, --Help: help
    -w, --Walltime: job run time in hours (96)
    -i, --Iterations: number of simulation runs
    """)

try:
    opts, args = getopt.getopt(sys.argv[1:], 'hw:i:', ['Help', 'Walltime=', 'Iterations='])
except getopt.GetoptError:
    Usage()
    sys.exit(2)
else:
    for opt, val in opts:
        if opt in ('-h', '--Help'):
            Usage()
            sys.exit(2)
        elif opt in ('-w', '--Walltime'):
            Walltime = val
        elif opt in ('-i', '--Iterations'):
            Iterations = int(val)

# make a list of species
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Macaque', 'Marmoset',
           'Hedgehog', 'Shrew', 'Cat', 'Dog', 'Mouse', 'Cow', 'Horse', 'Sloth',
           'Armadillo', 'Opossum', 'Platypus']
# make a list of models to use
Models = ['randomization', 'extension']

# write shell script to launch the job
# loop over models
for model in Models:
    # lopp over species
    for species in Species:
        # build outputfile name
        outputfile = 'Perform' + model.title() + 'Simulations' + species + '.sh' 
        newfile = open(outputfile, 'w')
        newfile.write('cd /RQusagers/rjovelin/richard/Richard/Nested_Genes/')
        newfile.write('\n')
        newfile.write('python3 NeutralSimulations.py ' + model + ' ' + species + ' ' + str(Iterations))
        newfile.close()            
        MyCommand = 'qsub -l walltime={0}:00:00,nodes=1:ppn=12,mem=46g '.format(Walltime) + outputfile
        os.system(MyCommand)
        print(MyCommand)