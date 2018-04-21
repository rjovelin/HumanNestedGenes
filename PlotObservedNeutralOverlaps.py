# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:13:23 2018

@author: rjovelin
"""

# use this script to plot numbers of observed and expected overlaps

# usage PlotObservedNeutralOverlaps.py NeutralSimsFile
# [NeutralSimsFile] file with results of neutral simulations

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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

#NeutralSimsFile = sys.argv[1]
NeutralSimsFile = 'NeutralSimulations.txt'

# read table with results of neutral simulations
Neutral = {}
infile = open(NeutralSimsFile)
header = infile.readline().rstrip().split('\t')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # get species name
        species = line[0]
        # get # observed overlapping genes
        observed = int(line[1])
        # get # expected genes and SEM for each model
        Shuffle, Extension = line[3].split(), line[6].split()
        MeanShuffle = float(Shuffle[0])
        SEMShuffle = float(Shuffle[1][Shuffle[1].index("(")+1:Shuffle[1].index(")")])
        MeanExtension = float(Extension[0])
        SEMExtension = float(Extension[1][Extension[1].index("(")+1:Extension[1].index(")")]) 
        Neutral[species] = [MeanShuffle, SEMShuffle, MeanExtension, SEMExtension, observed]
infile.close()

# create figure
figure = plt.figure(1, figsize = (3.5, 2))
# add a plot to figure (N row, N column, plot N)
ax = figure.add_subplot(1, 1, 1)

# set colors
colorscheme = ['#8856a7', '#f03b20', '#43a2ca', '#fee391', '#74c476']

# make a list of species               
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Macaque', 'Marmoset',
           'Mouse', 'Cat', 'Dog', 'Cow', 'Horse', 'Hedgehog', 'Shrew',
           'Sloth', 'Armadillo', 'Opossum', 'Platypus']               

# sort species according to observed values
SpObs =[]
for i in Neutral:
    SpObs.append([Neutral[i][4], i])
SpObs.sort()
print(SpObs)
Species = [i[1] for i in SpObs]               
      
# plot expected number of genes for shuffling model
ax.errorbar([i/10 for i in range(len(Species))], [Neutral[species][0] for species in Species],
             yerr = [[-Neutral[species][1] for species in Species], [Neutral[species][1] for species in Species]],
                     fmt = 'o', linewidth = 0.5, elinewidth = 0.5, color = '#43a2ca', ecolor = '#43a2ca', markersize = 5, marker = '^')
# plot expected number of genes for extension model
ax.errorbar([i/10 for i in range(len(Species))], [Neutral[species][2] for species in Species],
             yerr = [[-Neutral[species][3] for species in Species], [Neutral[species][3] for species in Species]],
                     fmt = 'o', linewidth = 0.5, elinewidth = 0.5, color = '#74c476', ecolor = '#74c476', markersize = 5, marker = 's')
# plot observed number of genes
ax.errorbar([i/10 for i in range(len(Species))], [Neutral[species][4] for species in Species],
             fmt = 'o', linewidth = 0.5, elinewidth = 0.5, color = '#f03b20', ecolor = '#f03b20', markersize = 5, marker = 'o')

# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Number of overlapping genes', color = 'black',  size = 7, ha = 'center', **FigFont)
plt.xticks([i/10 for i in range(17)], Species, rotation = 30, size = 7,
           color = 'black', ha = 'right', **FigFont)
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(True)  
# edit tick parameters    
plt.tick_params(axis='both', which='both', bottom='on', top='off',
                right = 'off', left = 'on', labelbottom='on',
                colors = 'black', labelsize = 7, direction = 'out')  
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')   

## add legend
O = mlines.Line2D([], [], color='#f03b20', marker='o', linestyle='None', markersize=4.5, label='observed')
S = mlines.Line2D([], [], color='#74c476', marker='s', linestyle='None', markersize=4.5, label='extension')
R = mlines.Line2D([], [], color='#43a2ca', marker='^', linestyle='None', markersize=4.5, label='shuffle')
ax.legend(handles=[O, S, R], loc = (0.02, 1.05), fontsize = 7, frameon = False, ncol = 3)

# save figure
for extension in ['.pdf', '.eps', '.png']:
    figure.savefig('NeutralSims' + extension, bbox_inches = 'tight')