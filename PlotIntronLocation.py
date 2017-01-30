# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 08:47:41 2017

@author: Richard
"""



# use this script to plot the distribution of intron location for external genes
# containing internal genes

# usage python3 PlotIntronLocation.py 

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
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

# load dictionaries of host and nested genes 
with open('HumanNestedGenes.json') as human_json_data:
    Nested = json.load(human_json_data)

# get GFF file
GFF = 'Homo_sapiens.GRCh38.86.gff3'

# get the coordinates of genes on each chromo
# {chromo: {gene:[chromosome, start, end, sense]}}
GeneChromoCoord = ChromoGenesCoord(GFF)
# map each gene to its mRNA transcripts
MapGeneTranscript = GeneToTranscripts(GFF)
# remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
# get the coordinates of each gene {gene:[chromosome, start, end, sense]}
GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
# Map Transcript names to gene names {transcript: gene}
MapTranscriptGene = TranscriptToGene(GFF)
# get the coordinates of all exons    
ExonCoord = GeneExonCoord(GFF)
ExonCoord = CleanGeneFeatureCoord(ExonCoord, MapTranscriptGene)
# get the intron coordinates of all transcripts {transcript: [[intron_start, intron_end]]}
IntronCoord = GeneIntronCoord(ExonCoord)
IntronCoord = CleanGeneFeatureCoord(IntronCoord, MapTranscriptGene)
# get the coordinates of each transcripts {transcript:[chromosome, start, end, sense]}
TranscriptCoordinates = TranscriptsCoord(GFF)
# map genes to their longest transcript {gene: longest_transcript}
GeneLongestTranscript = LongestTranscript(TranscriptCoordinates, MapGeneTranscript)
# match longest transcript of the nested genes to transcript of the host gene (longest transcript in priority)
Matches = MatchHostTranscriptWithNestedTranscript(Nested, MapGeneTranscript, GeneLongestTranscript, TranscriptCoordinates, IntronCoord)
# make a set of un-nested genes
UnNestedGenes = MakeNonOverlappingGeneSet(Nested, GeneCoord)
# list all host, nested transcript pairs [[host, nested]]
HostNestedPairs = GetHostNestedPairs(Matches)




# 1) plot the distribution of intron number per external gene 

# make a set of external transcripts
ExternalTS = list(Matches.keys())
# make a set of external transcripts with intronless internal genes
# make a set of external transcripts with intron-containing genes
ExternalIntronless, ExternalWithIntron = set(), set()
for i in range(len(HostNestedPairs)):
    # get the external and internal transcripts of the pair
    external, internal = HostNestedPairs[i][0], HostNestedPairs[i][1]
    # check if internal transcript has intron
    if internal in IntronCoord:
        ExternalWithIntron.add(external)
    else:
        ExternalIntronless.add(external)

# count the number of intron per external gene for all genes, external genes
# with intronless internal genes and external genes with intron-containing internal genes
TotalCount, IntronlessCount, WithIntronCount = [], [], []
for gene in ExternalTS:
    TotalCount.append(len(IntronCoord[gene]))
for gene in ExternalIntronless:
    IntronlessCount.append(len(IntronCoord[gene]))
for gene in ExternalWithIntron:
    WithIntronCount.append(len(IntronCoord[gene]))


# 2) plot the distribution of intron position of the external introns containing internal genes
# separately when internal genes are intronless or with introns

IntronlessPos, WithIntronPos = [], []
for i in range(len(HostNestedPairs)):
    # get the external and internal transcripts of the pair
    external, internal = HostNestedPairs[i][0], HostNestedPairs[i][1]
    assert TranscriptCoordinates[external][0] == TranscriptCoordinates[internal][0]
    # find the position of the internal gene in its host gene
    # get the start position of the internal gene
    InternalStart = TranscriptCoordinates[internal][1]
    # set boolean to check if intron is found
    FoundIntron = False    
    # check orientation of the external gene
    assert TranscriptCoordinates[external][-1] in ['-', '+']    
    if TranscriptCoordinates[external][-1] == '-':
        # 1st intron of the list is last intron of transcript
        # need to loop trhough introns in reversed order        
        for j in range(len(IntronCoord[external]) -1, -1, -1):
            # check if intron contains the internal gene
            if InternalStart in range(IntronCoord[external][j][0], IntronCoord[external][j][1]):
                # record position (1-based)
                # check if internal gene is ontronless or not
                if internal in IntronCoord:
                    WithIntronPos.append(len(IntronCoord[external]) -j)
                else:
                    IntronlessPos.append(len(IntronCoord[external]) -j)
            # update boolean
            FoundIntron = True
            # verify that internal gene is not in tron once intron has been found
            if FoundIntron == False:
                assert InternalStart not in range(IntronCoord[external][j][0], IntronCoord[external][j][1])
    if TranscriptCoordinates[external][-1] == '+':
        # loop over introns of the external gene    
        for j in range(len(IntronCoord[external])):
            # check if intron contains the internal gene
            if InternalStart in range(IntronCoord[external][j][0], IntronCoord[external][j][1]):
                # record position (1-based)
                # check if internal gene is ontronless or not
                if internal in IntronCoord:
                    WithIntronPos.append(j+1)
                else:
                    IntronlessPos.append(j+1)
                # update boolean
                FoundIntron = True
            # verify that internal gene is not in tron once intron has been found
            if FoundIntron == False:
                assert InternalStart not in range(IntronCoord[external][j][0], IntronCoord[external][j][1])
  

#avecintron = {}
#sansintron = {}
#for i in IntronlessPos:
#    if i in sansintron:
#        sansintron[i] += 1
#    else:
#        sansintron[i] = 1
#for i in WithIntronPos:
#    if i in avecintron:
#        avecintron[i] += 1
#    else:
#        avecintron[i] = 1
#a = [[key, val] for key, val in avecintron.items()]
#b = [[key, val] for key, val in sansintron.items()]        
#a.sort()
#b.sort()
#print(a)
#print(b)
#print('intronless', max(IntronlessPos))        
#print('with intron', max(WithIntronPos))        
        
    
# sort lists
WithIntronPos = np.sort(WithIntronPos)
IntronlessPos = np.sort(IntronlessPos)
# compute probabilities
PWithPos = np.array(range(len(WithIntronPos))) / len(WithIntronPos)
PNonePos = np.array(range(len(IntronlessPos))) / len(IntronlessPos) 


## create a list of lists with intron numbers/length for hosts, nested and un-nested genes
#IntronNumbers = [HostNum, NestedNum, OthersNum]
## create a list of lists with intron length for hosts, nested and un-nested genes
#IntronLength = [HostLength, NestedLength, OthersLength]
## create a list of lists with intron length of host genes with and without nested genes
#HostIntrons = [WithGene, WithoutGene]
#
#
#
## create a function to get the mean and SEM of items in a list
#def GetMeanSEM(L):
#    '''
#    (list) -> (list, list)
#    Take a list of inner lists of numbers and return a list with mean values
#    and a parallel list with SEM values for each item of the outter list
#    '''
#    # create lists of mean and SEM
#    MeanVal, SEMVal = [], []
#    # loop over the outter ist
#    for i in range(len(L)):
#        MeanVal.append(np.mean(L[i]))
#        SEMVal.append(np.std(L[i]) / math.sqrt(len(L[i])))
#    return MeanVal, SEMVal
#
#
#    
## create lists with means and with SEM
#NumMeans, NumSEM = GetMeanSEM(IntronNumbers)
#LengthMeans, LengthSEM = GetMeanSEM(IntronLength)
#HostIntronMeans, HostIntronSEM = GetMeanSEM(HostIntrons)
#
## perform statistical tests between gene categories
## create dict to store results
## {number or length: [P_host-nested, P_host-unnested, P_nested-unnested], host: [P_withgene_nogene]
#PValues = {}
## initialize dict with empty list
#PValues['number'] = []
## loop list of intron numbers
#for i in range(0, len(IntronNumbers) -1):
#    for j in range(i+1, len(IntronNumbers)):
#        P = stats.ranksums(IntronNumbers[i], IntronNumbers[j])[1]
#        PValues['number'].append(P)
#PValues['length'] = []
## loop list of intron length
#for i in range(0, len(IntronLength) -1):
#    for j in range(i+1, len(IntronLength)):
#        P = stats.ranksums(IntronLength[i], IntronLength[j])[1]
#        PValues['length'].append(P)
#PValues['host'] = []
## loop list of host intron length
#for i in range(0, len(HostIntrons) -1):
#    for j in range(i+1, len(HostIntrons)):
#        P = stats.ranksums(HostIntrons[i], HostIntrons[j])[1]
#        PValues['host'].append(P)
#
## print p values
#for i in PValues:
#    print(i, PValues[i])
#
#
## create a dict with significance level as stars
#Significance = {}
#for i in PValues:
#    # initialize dict with empty list
#    Significance[i] = [] 
#    # get the significance level
#    for pval in PValues[i]:
#        if pval >= 0.05:
#            Significance[i].append('')
#        elif pval < 0.05 and pval >= 0.01:
#            Significance[i].append('*')
#        elif pval < 0.01 and pval >= 0.001:
#            Significance[i].append('**')
#        elif pval < 0.001:
#            Significance[i].append('***')
#  
#


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, YLabel, Colors, YMax):
    '''
    (int, int, int, figure_object, list, list, list, list, list, str, str, int)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, 2 lists of data, a list of bar positions, the list of tick
    positions and their labels, a list of colors, a label for the Y axis,
    a maximum value for the Y axis and return an ax instance in the figure
    '''    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data    
    if type(Data[0]) == list:
        ax.hist(Data, bins = np.arange(0, max([max(Data[0]), max(Data[1])]) + 1, 1), linewidth = 0.7, histtype='bar', stacked=True)
    else:
        ax.hist(Data, bins = np.arange(0, max(Data) + 1, 1), facecolor= Colors, linewidth = 0.7)
    
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y
    ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    # add a range for the Y axis
    plt.ylim([0, YMax])
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(True)  
    # edit tick paramters
    plt.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
                    left = 'on', labelbottom='on', colors = 'black', labelsize = 8,
                    direction = 'out')  
    # add ticks on the x axis
    #plt.xticks(TickPos, Ticklabel)    
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    # create a margin around the x axis
    plt.margins(0.1)
    return ax      




## create figure
#fig = plt.figure(1, figsize = (4.5, 2.5))
#
## plot data for intron numner
#ax1 = CreateAx(3, 1, 1, fig, NumMeans, NumSEM, [0, 0.2, 0.4], [0.1, 0.3, 0.5], ['Hst', 'Nst', 'Un'], ['grey','black','white'], 'Number of introns per gene', 20)
#ax2 = CreateAx(3, 1, 2, fig, LengthMeans, LengthSEM, [0, 0.2, 0.4], [0.1, 0.3, 0.5], ['Hst', 'Nst', 'Un'], ['grey','black','white'], 'Intron length (Kbp)', 20)
#ax3 = CreateAx(3, 1, 3, fig, HostIntronMeans, HostIntronSEM, [0, 0.2], [0.1, 0.3], ['With', 'None'], ['grey','black'], 'Intron length (Kbp)', 70)
#
## make lists with bracket and star positions
#XPosNum = [[0.1, 0.28, 16, 0.2, 16.7], [0.1, 0.5, 17, 0.3, 18.2], [0.32, 0.5, 16, 0.4, 16.7]]
#XPosLength = [[0.1, 0.28, 16, 0.2, 16.7], [0.1, 0.5, 17, 0.3, 18.2], [0.32, 0.5, 16, 0.4, 16.7]]
#XPosHost = [[0.1, 0.28, 63, 0.2, 67]]
#    
## annotate figure to add significance
#for i in range(len(Significance['number'])):
#    if Significance['number'][i] != '':
#        ax1 = AddSignificance(ax1, Significance['number'][i], XPosNum[i][0], XPosNum[i][1], XPosNum[i][2], XPosNum[i][3], XPosNum[i][4])
#for i in range(len(Significance['length'])):
#    if Significance['length'][i]  != '':
#        ax2 = AddSignificance(ax2, Significance['length'][i], XPosLength[i][0], XPosLength[i][1], XPosLength[i][2], XPosLength[i][3], XPosLength[i][4])
#for i in range(len(Significance['host'])):
#    if Significance['host'][i] != '':
#        ax3 = AddSignificance(ax3, Significance['host'][i], XPosHost[i][0], XPosHost[i][1], XPosHost[i][2], XPosHost[i][3], XPosHost[i][4])
#
#
## add subplot labels
#ax1.text(-0.35, 21.5, 'A', horizontalalignment='center', verticalalignment='center',
#         color = 'black', fontname = 'Arial', size = 9)
#ax1.text(0.8, 21.5, 'B', horizontalalignment='center', verticalalignment='center',
#         color = 'black', fontname = 'Arial', size = 9)
#ax1.text(2.1, 21.5, 'C', horizontalalignment='center', verticalalignment='center',
#         color = 'black', fontname = 'Arial', size = 9)
#
## make sure subplots do not overlap
#plt.tight_layout()
#
## one can control padding between subplots with w_pad and h_pad 
##plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#
### padding between subplots can also be controlled with gridspec
##gs = gridspec.GridSpec(1, 3) # N rows and columns
##gs.update(wspace=0.3, hspace=0) # set the spacing between axes. 
#
## save figure
#fig.savefig('truc.pdf', bbox_inches = 'tight')






##################### cdf plots

# create figure
fig = plt.figure(1, figsize = (3.5, 4))


ax1 = CreateAx(1, 2, 1, fig, TotalCount, 'Number of introns per gene', 'lightgrey', max(TotalCount))
ax2 = CreateAx(1, 2, 2, fig, [WithIntronPos, IntronlessPos] , 'Position of gene-containing introns', ['black', 'lightgrey'], max([max(WithIntronPos), max(IntronlessPos)]))





def CombineHighValues(L, cutoff):
    '''
    (list, int) -> list
    Take a list of overlap length and return a modified list with values higher
    than cutoff equal to cutoff
    '''
    for i in range(len(L)):
        if L[i] >= cutoff:
            L[i] = cutoff
    return L



# make sure subplots do not overlap
plt.tight_layout()


    
fig.savefig('truc.pdf', bbox_inches = 'tight')

