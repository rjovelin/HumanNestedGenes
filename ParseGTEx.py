# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 18:37:02 2016

@author: RJovelin
"""


# use this script to parse the GTEX data


# match the sample ID in the metadata file with the sample ID in the normalized count file
# need to replace '.' with '-' in the sample ID name
# match tissue with sample iD
# merge tissue subtypes into types (eg, brain, heart, ...)
# for each gene get the sample size in each tissue
# for each gene get the median expression in each tissue
# this represent the exprefofile of the gene
# take the relative expression
# from the relative expression, can compute euclidian distance between genes
# remove 
