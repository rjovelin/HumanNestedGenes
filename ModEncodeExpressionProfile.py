# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 17:13:09 2017

@author: RJovelin
"""

# use this script to write a file with median normalized FPKM in each tissue for each gene

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
import json
import random
import copy
import sys
import os
import math
import numpy as np
from scipy import stats
from HsaNestedGenes import *



# make a list of mouse tissues matching human tissues with expression
Tissues = ['Adipose_tissue', 'Subcutaneous_adipose_tissue	',
           'Adrenal_gland', 'Urinary_bladder', 'Brain', 'Cerebellum',
           'Cortical_plate', 'Frontal_cortex', 'Olfactory_bulb',
           	'Mammary_gland', 	'Colon', 'Sigmoid_colon', 'Large_intestine',
            'Heart', 'Kidney', 'Liver', 'Lung', 'Ovary', 'Pancreas',
            'Small_intestine', 'Spleen', 'Stomach', 'Testis']								
	
# create a dictionary with tissue as key and a list of dictionaries with gene expression as value
Expression = {}
 
# move to directory with sample files for each tissue
os.chdir('Mouse_Expression') 

# loop over tissues, compute expression for each sample, store the dictionary
for tissue in Tissues:
    L = ComputeTissueExpression(tissue)
    Expression[tissue] = L
    
# Match mouse tissues to human tissues
HumanMatches = {'Adipose Tissue': ['Adipose_tissue', 'Subcutaneous_adipose_tissue'],
                'Adrenal Gland': ['Adrenal_gland'],
                'Bladder': ['Urinary_bladder'], 
                'Brain': ['Brain', 'Cerebellum', 'Cortical_plate', 'Frontal_cortex', 'Olfactory_bulb'],
                'Breast': ['Mammary_gland'],							
                'Colon': ['Colon', 'Sigmoid_colon', 'Large_intestine'],					
                'Heart': ['Heart'],						
                'Kidney': ['Kidney'],								
                'Liver': ['Liver'],								
                'Lung': ['Lung'],								
                'Ovary': ['Ovary']	,							
                'Pancreas': ['Pancreas'],								
                'Small Intestine':	 ['Small_intestine'],								
                'Spleen': ['Spleen'],								
                'Stomach': ['Stomach'],							
                'Testis': ['Testis']}								
		
# create a dict with the tissue names matching human tissues and list with mouse gene expression
MouseExpression = {}
for tissue in HumanMatches:
    if len(HumanMatches[tissue]) == 1:
        # map list of mouse expression to human tissue name
        MouseExpression[tissue] = Expression[HumanMatches[tissue]]
    else:
        # merge lists of dictionaries in a single list
        MouseExpression[tissue] = []
        for subtissue in HumanMatches[tissue]:
            MouseExpression[tissue].extend(Expression[subtissue])

# Merge all samples for a given tissue
for tissue in MouseExpression:
    MouseExpression[tissue] = MergeSamples(MouseExpression[tissue])

# take the median expression for a given gene across the different samples
for tissue in MouseExpression:
    for gene in MouseExpression[tissue]:
        MouseExpression[tissue][gene] = np.median(MouseExpression[tissue][gene])

# do some qc
for tissue in MouseExpression:
    print(tissue, len(MouseExpression[tissue]))







#Human Tissue          Mouse Tissues
#Adipose Tissue		Adipose_tissue, Subcutaneous_adipose_tissue							
#Adrenal Gland		Adrenal_gland								
#Bladder			Urinary_bladder								
#Brain			Brain, Cerebellum, Cortical_plate, Frontal_cortex, Olfactory_bulb	
#Breast			Mammary_gland								
#Colon			Colon	Sigmoid_colon, Large_intestine					
#Heart			Heart								
#Kidney			Kidney								
#Liver			Liver								
#Lung				Lung								
#Ovary			Ovary								
#Pancreas			Pancreas								
#Small Intestine		Small_intestine								
#Spleen			Spleen								
#Stomach			Stomach								
#Testis			Testis								
										

##############################
