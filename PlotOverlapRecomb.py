# -*- coding: utf-8 -*-
"""
Created on Wed May 31 16:51:25 2017

@author: RJovelin
"""

# use this script to test for enrichement of overlapping genes in regions of high recombination

# usage PlotOverlapRecomb.py

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


# load dictionaries of overlapping genes
JsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',
             'HumanPiggyBackGenes.json', 'HumanConvergentGenes.json',
             'HumanDivergentGenes.json']
# make a list of dictionaries
Overlap = []
# loop over files
for i in range(len(JsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(JsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    Overlap.append(overlapping)

# get GFF file
GFF = 'Homo_sapiens.GRCh38.88.gff3'
# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)

# generate gene sets
GeneSets = []
for i in range(len(Overlap)):
    GeneSets.append(MakeFullPartialOverlapGeneSet(Overlap[i]))
# make a set of non-overlapping genes
NonOverlappingGenes = MakeNonOverlappingGeneSet(Overlap[0], GeneCoord)


# get the coordinates of recombination hot and cold spots on each chromo
RecombSpots = GetColdHotRecomSport('CS_HRR_autosomes_final.bed')

# assign overlapping genes to hot and cold spots
OverlappingHotColdSpots = []
for i in range(len(GeneSets)):
    OverlappingHotColdSpots.append(AssignGenesToRecombSpots(GeneSets[i], RecombSpots, GeneCoord))
# assign non-overlapping genes to hot and cold spots
NonOverlappingHotColdSpots = AssignGenesToRecombSpots(NonOverlappingGenes, RecombSpots, GeneCoord)

# make lists of hot and cold spots [[overlap hot, non-overlap hot], [overlap cold, non-overlap cold]]
GetL = lambda x: len(x)
AllHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[0]), map(GetL, NonOverlappingHotColdSpots))) 
NestedHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[1]), map(GetL, NonOverlappingHotColdSpots)))
PbkHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[2]), map(GetL, NonOverlappingHotColdSpots)))
ConHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[3]), map(GetL, NonOverlappingHotColdSpots)))
DivHotColdSpots = list(zip(map(GetL, OverlappingHotColdSpots[4]), map(GetL, NonOverlappingHotColdSpots)))

#for L in [AllHotColdSpots, NestedHotColdSpots, PbkHotColdSpots, ConHotColdSpots, DivHotColdSpots]:
#    print(L[0][0]/sum(L[0]), L[1][0]/sum(L[1]), stats.fisher_exact(L)[1])

# compute odds ratios and standard error
OddsRatios, SEOdds = [], []
# Compute 
for L in [AllHotColdSpots, NestedHotColdSpots, PbkHotColdSpots, ConHotColdSpots, DivHotColdSpots]:
    OR = (L[0][0] * L[1][1]) / (L[0][1] * L[1][0])
    OddsRatios.append(OR)
    SE = math.sqrt((1 / L[0][0]) + (1/ L[0][1]) + (1/ L[1][0]) + (1/ L[1][1])) 
    SEOdds.append(SE)
# compute 95% confidence intervels
LowerConf, UpperConf = [], []
for i in range(len(OddsRatios)):
    # store distance between odds ratio and 95% CI
    LowerConf.append(OddsRatios[i] - math.exp(math.log(OddsRatios[i]) - (1.96 * SEOdds[i])))
    UpperConf.append(math.exp(math.log(OddsRatios[i]) + (1.96 * SEOdds[i])) - OddsRatios[i])

    
#names=data$method
#x = 1:13*2-1
#CI.up = as.numeric(data$mean)+as.numeric(data$ci)
#CI.dn = as.numeric(data$mean)-as.numeric(data$ci)
#plot(data$mean~x, cex=1.5,xaxt='n',ylim=c(0.3,0.40), xlab='',ylab='lalala!', main='blahblahblah',col='blue',pch=16)
#axis(1, at=x, labels=names)
#arrows(x,CI.dn,x,CI.up,code=3,length=0.2,angle=90,col='red')
#legend("bottomleft",paste(names,": S.E=",data$se),ncol=6,text.width=1)


#x = np.array([1, 2, 3, 4, 5])
#y = np.power(x, 2) # Effectively y = x**2
#e = np.array([1.5, 2.6, 3.7, 4.6, 5.5])
#plt.errorbar(x, y, e, linestyle='None', marker='^')


## example error bar values that vary with x-position
#error = 0.1 + 0.2 * x
## error bar values w/ different -/+ errors
#lower_error = 0.4 * error
#upper_error = error
#asymmetric_error = [lower_error, upper_error]


# create figure
figure = plt.figure(1, figsize = (1.5, 2))
# add a plot to figure (N row, N column, plot N)
ax = figure.add_subplot(1, 1, 1)

# set colors
colorscheme = ['#f03b20', '#43a2ca', '#fee391', '#74c476']

# plot odds ratio and confidence intervals
#ax.errorbar([0.1, 0.2, 0.3, 0.4, 0.5], OddsRatios, yerr = [LowerConf, UpperConf], fmt = 'o',
#            color = colorscheme, linewidth = 0.7, ecolor = colorscheme, elinewidth = 0.7, capsize = None,
#            capthick = None)
            
ax.errorbar([0.1, 0.2, 0.3, 0.4, 0.5], OddsRatios, yerr = [LowerConf, UpperConf], fmt = 'o',
            linewidth = 0.7, elinewidth = 0.7)
# set Y axis limits
ax.ylim([0, 1.2])


# set font for all text in figure
FigFont = {'fontname':'Arial'}   
# write y axis label
ax.set_ylabel('Odds ratios', color = 'black',  size = 7, ha = 'center', **FigFont)
# add ticks and lebels
plt.xticks([0.1, 0.2, 0.3, 0.4, 0.5], ['All', 'Nst', 'Pbk', 'Conv', 'Div'], rotation = 0, size = 7, color = 'black', ha = 'right', **FigFont)
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
   
   
  
#   
#       edgecolor = 'black', linewidth = 0.7,
#       error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth = 0.7))


# add margins
plt.margins(0.1)
    



# save figure
figure.savefig('truc.pdf', bbox_inches = 'tight')



#fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
#ax0.errorbar(x, y, yerr=error, fmt='-o')
#ax0.set_title('variable, symmetric error')
#
#ax1.errorbar(x, y, xerr=asymmetric_error, fmt='o')
#ax1.set_title('variable, asymmetric error')
#ax1.set_yscale('log')
#plt.show()
#    


#https://stackoverflow.com/questions/33282368/plotting-a-2d-heatmap-with-matplotlib
#http://seaborn.pydata.org/generated/seaborn.heatmap.html
#https://matplotlib.org/examples/pylab_examples/pcolor_small.html