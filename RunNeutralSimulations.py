# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:33:35 2018

@author: RJovelin
"""

import os

# use this script to generate shell scripts to run neutral simulations and launch jobs

# make a list of species
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Macaque', 'Marmoset',
           'Hedgehog', 'Shrew', 'Cat', 'Dog', 'Mouse', 'Cow', 'Horse', 'Sloth',
           'Armadillo', 'Opossum', 'Platypus']
# make a list of models to use
Models = ['randomization', 'extension']
# set number of simulations
Iterations = 1000

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
        MyCommand = 'qsub -l walltime=30:00:00,nodes=1:ppn=12,mem=46g ' + outputfile
        os.system(MyCommand)
        print(MyCommand)
        