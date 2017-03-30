# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 23:42:40 2017

@author: Richard
"""



# use this script to download GFF files (for all species) and CDS files (for human, chimp and mouse) from ensembl
# Note that web links were correct as of March 2017



import os

# create a list of web links
WebLinks = []









# loop of urls and download the corresponding file
for link in WebLinks:
    os.system('wget ' + link)
    

