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




# get the UP-FPKM for each file separately
# remove the genes with 0 counts
# compute UP
# get the normalization



Adipose Tissue				Adipose_tissue	Subcutaneous_adipose_tissue							
Adrenal Gland				Adrenal_gland								
Bladder				Urinary_bladder								
Blood												
Blood Vessel												
Brain				Brain	Cerebellum		Cortical_plate		Frontal_cortex		Olfactory_bulb	
Breast				Mammary_gland								
Cervix Uteri												
Colon				Colon	Sigmoid_colon		Large_intestine					
Esophagus												
Fallopian Tube												
Heart				Heart								
Kidney				Kidney								
Liver				Liver								
Lung				Lung								
Muscle												
Nerve												
Ovary				Ovary								
Pancreas				Pancreas								
Pituitary												
Prostate												
Salivary Gland												
Skin												
Small Intestine				Small_intestine								
Spleen				Spleen								
Stomach				Stomach								
Testis				Testis								
Thyroid												
Uterus												
Vagina												


