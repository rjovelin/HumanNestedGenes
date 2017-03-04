# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:30:11 2017

@author: RJovelin
"""


# use this script to plot the distances between non-overlapping orthologs of human overlapping genes

# import modules
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
import matplotlib.gridspec as gridspec
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
jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json']

# make a list of dictionaries
AllOverlap = []
# loop over files
for i in range(len(jsonFiles)):
    # load dictionary of overlapping gene pairs
    json_data = open(jsonFiles[i])
    overlapping = json.load(json_data)
    json_data.close()
    AllOverlap.append(overlapping)

# get GFF file
GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3']

# make a list of gene coordinates       
AllCoordinates, AllOrdered = [], []
# loop over GFF files
for i in range(len(GFF)):
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    AllCoordinates.append(GeneCoord)
    AllOrdered.append(OrderedGenes)
HumanOrdered, ChimpOrdered = AllOrdered[0], AllOrdered[1]
HumanCoord, ChimpCoord = AllCoordinates[0], AllCoordinates[1]

# get 1:1 orthologs between human and chimp
Orthos = MatchOrthologPairs('HumanChimpOrthologs.txt')

# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# get the gene pairs
HumanPairs = AllPairs[:2]
ChimpPairs = AllPairs[2:]

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos or pair[1] not in Orthos:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
# remove order for chimp gene pairs
for i in range(len(ChimpPairs)):
    for j in range(len(ChimpPairs[i])):
        ChimpPairs[i][j] = set(ChimpPairs[i][j])



# 1) plot the proportions of orthologous gene pairs that are adjacent, separated
# and on different chromosomes and nested
# 2) plot a histogram of the intergenic distance between non-nested adjacent
# chimp orthologs of human nested genes

SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG, NstCons = 0, 0, 0, 0, 0, 0

# get the intergenic distance between non-nested adjacent orthologs of human nested genes
GeneDist = []  
 
# loop over human nested gene pairs
for i in range(len(HumanPairs[1])):
    # get the orthologs of the human genes
    ortho1, ortho2 = Orthos[HumanPairs[1][i][0]], Orthos[HumanPairs[1][i][1]]
    # check if human gene pairs is nested in chimp
    if set([ortho1, ortho2]) not in ChimpPairs[1]:
        # get chromosomes of chimp orthologs
        chromo1, chromo2 = ChimpCoord[ortho1][0], ChimpCoord[ortho2][0]        
        # check if chromosomes are different
        if chromo1 != chromo2:
            DiffLG += 1
        else:
            # get start positions of chimp orthologs
            S1, S2 = ChimpCoord[ortho1][1], ChimpCoord[ortho2][1]
            # get end positions of chimp orthologs
            E1, E2 = ChimpCoord[ortho1][2], ChimpCoord[ortho2][2]            
            # get indices of chimp orthologs in the ordred gene list
            P1, P2 = ChimpOrdered[chromo1].index(ortho1), ChimpOrdered[chromo2].index(ortho2)
            D = 'na'            
            # check if orthologs are overlapping in chimp
            if set([ortho1, ortho2]) in ChimpPairs[0]:
                # chimp genes are overlapping
                # check if they are adjcent
                if S1 < S2:
                    assert P1 < P2
                    if P2 != P1 + 1:
                        # overlapping non adjacent
                        SepOvlp += 1
                    elif P2 == P1 + 1:
                        # overlapping adjacent
                        AdjOvlp += 1
                        # get the distance between genes
                        D = S2 - E1
                elif S1 > S2:
                    assert P2 < P1
                    if P1 != P2 + 1:
                        # overlapping non adjacent
                        SepOvlp += 1
                    elif P1 == P2 + 1:
                        # overlapping adjacent
                        AdjOvlp += 1
                        # get the distance between genes              
                        D = S1 - E2
                elif S1 == S2:
                    print('merde')
            elif set([ortho1, ortho2]) not in ChimpPairs[0]:
                # chimp genes are not overlapping
                # check if they are adjacent
                if S1 < S2:
                    assert P1 < P2
                    if P2 != P1 + 1:
                        # non-overlapping non-adjacent
                        SepNonOvlp += 1
                    elif P2 == P1 + 1:
                        # non-overlapping adjacent
                        AdjNonOvlp += 1
                        # get the distance between genes
                        D = S2 - E1
                elif S1 > S2:
                    assert P2 < P1
                    if P1 != P2 + 1:
                        # non-overlapping non adjacent
                        SepNonOvlp += 1
                    elif P1 == P2 + 1:
                        # non-overlapping adjacent
                        AdjNonOvlp += 1
                        # get the distance between genes              
                        D = S1 - E2
                elif S1 == S2:
                    print('shit')
            if D != 'na':
                GeneDist.append(D)  
    elif set([ortho1, ortho2]) in ChimpPairs[1]:
        # orthologs of human nested genes are also nested
        NstCons += 1            

# make a list of counts
PairCounts = [SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG, NstCons]
# make a list of labels
Labels = ['Separated overlapping', 'Separated non-overlapping', 'Adjacent overlapping', 'Adjacent non-overlapping', 'Different chromosomes', 'Nested']
# remove counts and labels equal to 0
to_remove = []
for i in range(len(PairCounts)):
    if PairCounts[i] == 0:
        to_remove.append(PairCounts[i])
        to_remove.append(Labels[i])
for i in to_remove:
    if i in PairCounts:
        PairCounts.remove(i)
    elif i in Labels:
        Labels.remove(i)
# associate counts and corresponding labels
CountsLabels = list(zip(PairCounts, Labels))  
# sort count, labels according to counts    
CountsLabels.sort()
PairCounts = [i[0] for i in CountsLabels]
Labels = [i[1] for i in CountsLabels]
# get proportions of gene pairs
Proportions = [str(round((i / sum(PairCounts)*100),1)) for i in PairCounts]

print(PairCounts)
print(Labels)
print(Proportions)


            
# 1) plot the proportions of orthologous gene pairs that are adjacent, separated and on different chromosomes


# 2) plot a histogram of the intergenic distance between non-nested adjacent chimp orthologs of human nested genes

# collapse values for | distance | > 100 Kb
for i in range(len(GeneDist)):
    if GeneDist[i] < -100000:
        GeneDist[i] = -110000
    elif GeneDist[i] > 100000:
        GeneDist[i] = 110000


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, GraphType):
    '''
    return an ax object
    '''    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data    
    if GraphType == 'histo':
        # get current axis
        currentAxis = plt.gca()
        # add rectangle to current axis
        currentAxis.add_patch(mpl.patches.Rectangle((-110000, 0), 110000, 70, linewidth = 0.7, fill = True, edgecolor = '#f7fcb9', facecolor = '#f7fcb9', alpha = 1))      
        currentAxis.add_patch(mpl.patches.Rectangle((0, 0), 110000, 70, linewidth = 0.7, fill = True, edgecolor = '#ffeda0', facecolor = '#ffeda0', alpha = 1))      
        # plot data
        ax.hist(Data, bins = np.arange(min(Data), max(Data)+10000, 10000), linewidth = 0.7, histtype='bar', fill = True, facecolor = '#0c2c84', edgecolor = 'black', alpha = 1)    
        # set font for all text in figure
        FigFont = {'fontname':'Arial'}   
        # annotate graph with subplot labels
        ax.text(-90000, 55, 'Overlapping', color = 'black', size = 7, fontname = 'Arial')
        ax.text(20000, 55, 'Non-overlapping', color = 'black', size = 7, fontname = 'Arial')
        # write label axis
        ax.set_ylabel('Proportion of nested gene pairs in human', color = 'black',  size = 7, ha = 'center', **FigFont)
        ax.set_xlabel('Intergenic distance in chimp (Kb)', color = 'black',  size = 7, ha = 'center', **FigFont)
        # set a limit to y axis
        plt.ylim([0, 70])
        # add ticks on the x axis
        TickPos = ['', '< -100', '', '', '', '', '-50',
                   '', '', '', '', '0', '', '', '',
                   '', '50', '', '', '', '', '> 100', '']
        plt.xticks(np.arange(-110000, 120000, 10000), TickPos, size = 7, color = 'black', ha = 'center', **FigFont)
        # do not show lines around figure  
        ax.spines["top"].set_visible(False)    
        ax.spines["bottom"].set_visible(True)    
        ax.spines["right"].set_visible(False)    
        ax.spines["left"].set_visible(True)  
        # edit tick paramters
        plt.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
                    left = 'on', labelbottom='on', colors = 'black', labelsize = 7, direction = 'out')  
        # Set the tick labels font name
        for label in ax.get_yticklabels():
            label.set_fontname('Arial')
        plt.margins(0)
    elif GraphType == 'donut':
        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        sizes, labels = Data[0], Data[1]
        explode = [0] * len(sizes) # "explode" slices        
        colors = ['lightgreen', 'gold', 'lightskyblue', '#9e9ac8', 'lightcoral']
        '%1.1f%%'
        patches, texts, autotexts = ax.pie(sizes, explode=explode, labels=labels, colors = colors, autopct='%1.1f%%',
                                           shadow=False, startangle=90, pctdistance=0.6, labeldistance=1.1, 
                                           counterclock=True)
        for i in range(len(texts)):
            texts[i].set_fontsize(7)
        
            
               
               
        # draw a circle at the center of pie to make it look like a donut
        centre_circle = plt.Circle((0,0),0.65,color='black', fc='white',linewidth=1)
        # modify line parameters of pie chart
        mpl.rcParams['patch.linewidth'] = 1  
        mpl.rcParams['patch.edgecolor'] = 'white' 
        # Equal aspect ratio ensures that pie is drawn as a circle
        ax.axis('equal')  
        # add circle to pie chart
        fig.gca().add_artist(centre_circle)
        
        # annotate graph with legends
        # create patches 
        nst = mpatches.Patch(facecolor = 'lightcoral', edgecolor = 'black', linewidth = 0.5, label= 'Nested', alpha = 1)
        adjNon = mpatches.Patch(facecolor = '#9e9ac8', edgecolor = 'black', linewidth = 0.5, label= 'Adjacent non-overlapping', alpha = 1)
        sepNon = mpatches.Patch(facecolor = 'lightskyblue', edgecolor = 'black', linewidth = 0.5, label= 'Separated non-overlapping', alpha = 1)
        adjO = mpatches.Patch(facecolor = 'gold', edgecolor = 'black', linewidth = 0.5, label= 'Adjacent overlapping', alpha = 1)
        diffC = mpatches.Patch(facecolor = 'lightgreen', edgecolor = 'black', linewidth = 0.5, label= 'Different chromosomes', alpha = 1)
        ax.legend(handles = [nst, adjNon, sepNon, adjO, diffC], bbox_to_anchor=(0.8, 0.8), loc = 3, fontsize = 7, frameon = False, ncol = 1)
            
        
    return ax      



# create figure
fig = plt.figure(1, figsize = (5.5, 2.5))
ax1 = CreateAx(2, 1, 1, fig, [PairCounts, Labels], 'donut')
ax2 = CreateAx(2, 1, 2, fig, GeneDist, 'histo')

# annotate graph with subplot labels
ax2.text(-400000, 75, 'A', color = 'black', size = 7, fontname = 'Arial')
ax2.text(-120000, 75, 'B', color = 'black', size = 7, fontname = 'Arial')









# make sure subplots do not overlap
plt.tight_layout()

fig.savefig('truc.pdf', bbox_inches = 'tight')








#
### annotate graphs with legends
## create patches 
#Title = mpatches.Patch(facecolor = 'none', edgecolor = 'none', linewidth = 0, label= 'Internal genes:', alpha = 0)
#Yes = mpatches.Patch(facecolor = 'orange', edgecolor = 'none', linewidth = 0.5, label= 'with introns', alpha = 0.5)
#No = mpatches.Patch(facecolor = 'blue', edgecolor = 'black', linewidth = 0.5, label= 'intronless', alpha = 0.5)
#ax1.legend(handles = [Title, Yes, No], bbox_to_anchor=(0.2, 0.4), loc = 3, fontsize = 5, frameon = False, ncol = 1)
#ax2.legend(handles = [Title, Yes, No], bbox_to_anchor=(0.2, 0.4), loc = 3, fontsize = 5, frameon = False, ncol = 1)
## create lines
#orange_line = mlines.Line2D([], [], color='orange', marker='', markersize=15, label='with introns', alpha = 0.5)
#blue_line = mlines.Line2D([], [], color='blue', marker='', markersize=15, label='intronless', alpha = 0.5)
#ax3.legend(handles = [Title, orange_line, blue_line], bbox_to_anchor=(0.05, 0.5), loc = 3, fontsize = 5, frameon = False, ncol = 1)
#
#
