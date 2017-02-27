# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:46:17 2017

@author: RJovelin
"""

# use this script to plot the % of orthologous gene pairs with same topology between human and chimp

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
jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json', 'HumanPiggyBackGenes.json',
             'HumanConvergentGenes.json', 'HumanDivergentGenes.json', 
             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json', 'ChimpPiggyBackGenes.json',
             'ChimpConvergentGenes.json', 'ChimpDivergentGenes.json',
             'MouseOverlappingGenes.json', 'MouseNestedGenes.json', 'MousePiggyBackGenes.json',
             'MouseConvergentGenes.json', 'MouseDivergentGenes.json']

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
GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3', 'Mus_musculus.GRCm38.86.gff3']

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
HumanOrdered, ChimpOrdered, MouseOrdered = AllOrdered[0], AllOrdered[1], AllOrdered[2]
HumanCoord, ChimpCoord, MouseCoord = AllCoordinates[0], AllCoordinates[1], AllCoordinates[2]

# get 1:1 orthologs between human and chimp
HsaPtrOrthos = MatchOrthologPairs('HumanChimpOrthologs.txt')
HsaMmuOrthos = MatchOrthologPairs('HumanMouseOrthologs.txt')

# make a list of dictionaries with orthologs
Orthos = [HsaPtrOrthos, HsaMmuOrthos]

# make pairs of overlapping genes
AllPairs = []
for i in range(len(AllOverlap)):
    pairs = GetHostNestedPairs(AllOverlap[i])
    AllPairs.append(pairs)
# get the gene pairs
HumanPairs = AllPairs[:5]
# copy human gene pairs
HsaPairs = copy.deepcopy(HumanPairs)
ChimpPairs = AllPairs[5:10]
MousePairs = AllPairs[10:]

# remove human genes lacking orthologs
for i in range(len(HumanPairs)):
    to_remove = []
    for pair in HumanPairs[i]:
        if pair[0] not in Orthos[0] or pair[1] not in Orthos[0]:
            to_remove.append(pair)
    for pair in to_remove:
        HumanPairs[i].remove(pair)
for i in range(len(HsaPairs)):
    to_remove = []
    for pair in HsaPairs[i]:
        if pair[0] not in Orthos[1] or pair[1] not in Orthos[1]:
            to_remove.append(pair)
    for pair in to_remove:
        HsaPairs[i].remove(pair)        


# 1) plot proportions of gene pairs with varying distance conserved in chimp and mouse
# 2) plot proportions of overlapping gene pairs with conserved topology in chimp and mouse
# 3) plot proportions of overlapping gene pairs that are overlapping in chimp and mouse
# 4) plot differences between conservation of human overalapping in chimop and mouse
#    and overlapping genes in chimp and mouse conserved in human human

# create a dict with lists of gene pairs
AllGenePairs = {}
AllGenePairs['hsa_ptr'] = copy.deepcopy(HumanPairs)
AllGenePairs['hsa_mmu'] = copy.deepcopy(HsaPairs)
AllGenePairs['ptr'] = copy.deepcopy(ChimpPairs)
AllGenePairs['mmu'] = copy.deepcopy(MousePairs)

# replace human genes by their orthologs
species = ['hsa_ptr', 'hsa_mmu']
for i in range(len(species)):
    for j in range(len(AllGenePairs[species[i]])):
        for k in range(len(AllGenePairs[species[i]][j])):
            AllGenePairs[species[i]][j][k][0] = Orthos[i][AllGenePairs[species[i]][j][k][0]]
            AllGenePairs[species[i]][j][k][1] = Orthos[i][AllGenePairs[species[i]][j][k][1]]

# sort gene pairs
for sp in AllGenePairs:
    for i in range(len(AllGenePairs[sp])):
        for j in range(len(AllGenePairs[sp][i])):
            AllGenePairs[sp][i][j].sort()

# replace gene lists with strings
for sp in AllGenePairs:
    for i in range(len(AllGenePairs[sp])):
        for j in range(len(AllGenePairs[sp][i])):
            AllGenePairs[sp][i][j] = ':'.join(AllGenePairs[sp][i][j])

# do qc
for sp in AllGenePairs:
    for i in range(1, len(AllGenePairs[sp])):
        for pair in AllGenePairs[sp][i]:
            assert pair in AllGenePairs[sp][0]

# make lists of adjacent gene pairs with varying distance in human conserved in chimp and mouse
# [[gene1, gene2], ....[gene n, gene n+1]]
for i in range(len(species)):
    HsaPairsDist = []    
    for j in range(7):
        HsaPairsDist.append([])
    # loop over chromosomes
    for chromo in HumanOrdered:
        # loop over the list of ordered genes
        for k in range(len(HumanOrdered[chromo]) - 1):
            # get the end position of gene 1
            EndGene1 = HumanCoord[HumanOrdered[chromo][k]][2]                
            # get the start position of adjacent gene 2
            StartGene2 = HumanCoord[HumanOrdered[chromo][k+1]][1]
            # computance distance between genes
            D = StartGene2 - EndGene1
            # assign infinity value to k
            m = float('inf')
            if D < 0:
                m = 0
            elif D >= 0 and D < 1000:
                m = 1
            elif D >= 1000 and D < 10000:
                m = 2                
            elif D >= 10000 and D < 50000:
                m = 3
            elif D >= 50000 and D < 100000:
                m = 4
            elif D >= 100000 and D < 150000:
                m = 5
            elif D >= 150000:
                m = 6
            if HumanOrdered[chromo][k] in Orthos[i] and HumanOrdered[chromo][k+1] in Orthos[i] and m in range(7):
                # add the human orthologs
                pair = [Orthos[i][HumanOrdered[chromo][k]], Orthos[i][HumanOrdered[chromo][k+1]]]
                # sort gene pair
                pair.sort()
                # concert list to string
                pair = ':'.join(pair)
                HsaPairsDist[m].append(pair)
    # add the pairs of non-overlapping genes to the lists of gene pairs
    AllGenePairs[species[i]].extend(HsaPairsDist)


# make list of adjacent gene pairs regardless of distance in chimp and mouse [gene1:gene2, ....gene_n:gene_n+1]
sp2 = ['ptr', 'mmu']
Sp2Coord = [ChimpCoord, MouseCoord]
Sp2Ordered = [ChimpOrdered, MouseOrdered]
for i in range(len(sp2)):
    Sp2PairsDist = []
    # loop over chromosomes
    for chromo in Sp2Ordered[i]:
        # loop over the list of ordered genes
        for j in range(len(Sp2Ordered[i][chromo]) - 1):
            # get the end position of gene 1
            EndGene1 = Sp2Coord[i][Sp2Ordered[i][chromo][j]][2]                
            # get the start position of adjacent gene 2
            StartGene2 = Sp2Coord[i][Sp2Ordered[i][chromo][j+1]][1]
            # check that genes are not overlapping
            D = StartGene2 - EndGene1
            if D >= 0:
                # get gene pair
                pair = [Sp2Ordered[i][chromo][j], Sp2Ordered[i][chromo][j+1]]
                # sort pair
                pair.sort()
                # convert list to string
                pair = ':'.join(pair)
                Sp2PairsDist.append(pair)
    AllGenePairs[sp2[i]].append(Sp2PairsDist)


# convert lists to numpy arrays
for sp in AllGenePairs:
    for i in range(len(AllGenePairs[sp])):
        AllGenePairs[sp][i] = np.array(AllGenePairs[sp][i])


# 1) plot proportions of gene pairs with varying distance conserved in chimp and mouse
# get the proportions of non-overlapping adjacent gene pairs with conserved topology
NonOverlap = []
for i in range(5, len(AllGenePairs['hsa_ptr'])):
    # count the number of adjacent gene pairs gene pairs that are adjacent in chimp
    total = sum(np.in1d(AllGenePairs['hsa_ptr'][i], AllGenePairs['ptr'][5], invert = False))
    NonOverlap.append(total / len(AllGenePairs['hsa_ptr'][i]))
    # count the number of adjacent gene pairs gene pairs that are adjacent in mouse
    total = sum(np.in1d(AllGenePairs['hsa_mmu'][i], AllGenePairs['mmu'][5], invert = False))
    NonOverlap.append((total / len(AllGenePairs['hsa_mmu'][i])) * 100)    

# 2) plot proportions of overlapping gene pairs with conserved topology in chimp and mouse
# get the proportions of overlapping gene pairs with conserved topology in each overlapping category
OverlapCat = []
for i in range(5):
    # count the number of overlapping gene pairs that are also overlapping in the same category in chimp
    total = sum(np.in1d(AllGenePairs['hsa_ptr'][i], AllGenePairs['ptr'][i], invert = False))
    OverlapCat.append(total / len(AllGenePairs['hsa_ptr'][i]))
    # count the number of overlapping gene pairs that are also overlapping in the same category in mouse
    total = sum(np.in1d(AllGenePairs['hsa_mmu'][i], AllGenePairs['mmu'][i], invert = False))
    OverlapCat.append((total / len(AllGenePairs['hsa_mmu'][i])) * 100)

# 3) plot proportions of overlapping gene pairs that are overlapping in chimp and mouse
# get the proportions of overlapping gene pairs that are also overlapping (regadless of category)
OverlapAll = []
for i in range(5):
    # count the number of overlapping gene pairs that are also overlapping in chimp
    total = sum(np.in1d(AllGenePairs['hsa_ptr'][i], AllGenePairs['ptr'][0], invert = False))
    OverlapAll.append(total / len(AllGenePairs['hsa_ptr'][i]))
    # count the number of overlapping gene pairs that are also overlapping in mouse
    total = sum(np.in1d(AllGenePairs['hsa_mmu'][i], AllGenePairs['mmu'][0], invert = False))
    OverlapAll.append((total / len(AllGenePairs['hsa_mmu'][i])) * 100)
    

#CountPairs = {}
#for i in range(len(species)):
#    for j in range(len(AllGenePairs[species[i]])):
#        if j < 5:
#            total = sum(np.in1d(AllGenePairs[species[i]][j], AllGenePairs[sp2[i]][j], invert = False))
#        else:
#            total = sum(np.in1d(AllGenePairs[species[i]][j], AllGenePairs[sp2[i]][5], invert = False)) 
#        if species[i] in CountPairs:
#            CountPairs[species[i]].append([total, len(AllGenePairs[species[i]][j])])
#        elif species[i] not in CountPairs:
#            CountPairs[species[i]] = [[total, len(AllGenePairs[species[i]][j])]]
#            
#  
#
## create a list of overlapping gene categories parallel to the list of overlapping pairs
#GeneCats = ['overlapping', 'nested', 'piggyback', 'convergent', 'divergent',
#            'proximal', 'moderate', 'intermediate', 'distant']
#
#for sp in CountPairs:
#    for i in range(len(CountPairs[sp])):
#        print(sp + '\t' + str(i) + '\t' + str(CountPairs[sp][i][0] / CountPairs[sp][i][1]))



        

# 4) plot differences between conservation of human overalapping in chimop and mouse
#    and overlapping genes in chimp and mouse conserved in human human

# make a reverse dictionary of orthologs 
ChimpOrthos = {}
for gene in HsaPtrOrthos:
    assert HsaPtrOrthos[gene] not in ChimpOrthos
    ChimpOrthos[HsaPtrOrthos[gene]] = gene
MouseOrthos = {}
for gene in HsaMmuOrthos:
    assert HsaMmuOrthos[gene] not in MouseOrthos
    MouseOrthos[HsaMmuOrthos[gene]] = gene

# remove gene pairs lacking orthologs
for i in range(len(ChimpPairs)):
    to_remove = []
    for pair in ChimpPairs[i]:
        if pair[0] not in ChimpOrthos or pair[1] not in ChimpOrthos:
            to_remove.append(pair)
    for pair in to_remove:
        ChimpPairs[i].remove(pair)
for i in range(len(MousePairs)):
    to_remove = []
    for pair in MousePairs[i]:
        if pair[0] not in MouseOrthos or pair[1] not in MouseOrthos:
            to_remove.append(pair)
    for pair in to_remove:
        MousePairs[i].remove(pair)


# add dictionaries of orthologs to list
Orthos.append(ChimpOrthos)
Orthos.append(MouseOrthos)

# make a list of gene pair lists
GenePairs = [HumanPairs, HsaPairs, ChimpPairs, MousePairs]
# make a list of sets of gene pairs
GeneSets = []
for i in range(len(GenePairs)):
    spsets = []
    for j in range(len(GenePairs[i])):
        spsets.append([set(k) for k in GenePairs[i][j]])
    GeneSets.append(spsets)

# make a list of counts of conserved and non-conserved human overlapping genes in chimp and mouse
Conserved = []
for i in range(len(GenePairs)):
    if i < 2:
        j = i + 2
    else:
        j = i - 2
    conservation = []
    for k in range(len(GenePairs[i])):
        conserved, divergent = 0, 0
        for pair in GenePairs[i][k]:
            if set([Orthos[i][pair[0]], Orthos[i][pair[1]]]) in GeneSets[j][k]:
                conserved += 1
            else:
                divergent += 1
        conservation.append([conserved, divergent])
    Conserved.append(conservation) 


# create lists with proportions of overlapping genes with conserved topologies
# [[prop human overlapping genes conserved, prop chimp overlapping genes conserved in human, etc]]
Proportions = []
for i in range(2):
    prop = []
    for j in range(len(Conserved[i])):
        man = Conserved[i][j][0] / sum(Conserved[i][j])
        sp2 = Conserved[i+2][j][0] / sum(Conserved[i+2][j])
        diff = (sp2 - man) * 100
        prop.append(diff)
    Proportions.append(prop)

  
# create a single list with differences between human and chimp and between
# human and mouse for each overlapping gene category
Differences = []
for i in range(len(Proportions[0])):
    Differences.append(Proportions[0][i])
    Differences.append(Proportions[1][i])




# create a function to format the subplots
def CreateAx(Columns, Rows, Position, figure, Data, Colors, BarPos, YLabel,
             XTickpos, XTicklabels, YTicksRange, YMin, YMax):
    '''
    Returns a ax instance in figure
    '''    

    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data
    ax.bar(BarPos, Data, width = 0.2, color = Colors, edgecolor = 'black', linewidth = 0.7)
    # draw x axis line
    if Position == 4:
        ax.plot([0, 2.8], [0, 0], lw = 0.7, color = 'black')
    
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write y axis label
    ax.set_ylabel(YLabel, color = 'black',  size = 7, ha = 'center', **FigFont)
    # add ticks and lebels
    plt.xticks(XTickpos, XTicklabels, rotation = 0, size = 7, color = 'black', ha = 'center', **FigFont)
    # edit y axis ticks
    plt.yticks(YTicksRange)
    # add a range for the Y and X axes
    plt.ylim([YMin, YMax])    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(True)  
    
    # make sure the y axis crosses the x axis at 0
    ax.spines['left'].set_position('zero')

    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    right = 'off', left = 'on', labelbottom='on',
                    colors = 'black', labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    
    ## add margins
    #plt.margins(0.1)
    
    return ax




# 1) plot proportions of gene pairs with varying distance conserved in chimp and mouse
# 2) plot proportions of overlapping gene pairs with conserved topology in chimp and mouse
# 3) plot proportions of overlapping gene pairs that are overlapping in chimp and mouse
# 4) plot differences between conservation of human overalapping in chimop and mouse
#    and overlapping genes in chimp and mouse conserved in human human



# create figure
fig = plt.figure(1, figsize = (5, 4))


#ax1 = CreateAx
#ax2 = CreateAx
#ax3 = CreateAx





# 1) plot proportions of gene pairs with varying distance conserved in chimp and mouse
# get the proportions of non-overlapping adjacent gene pairs with conserved topology
print(NonOverlap)
print(OverlapCat)
print(OverlapAll)




[0.68611111111111112, 0.54855518711511131, 0.92485119047619047, 0.76479400749063675, 0.93409019236833801, 0.88601872780109781, 0.91476872505282936, 0.86489306207616068, 0.87984981226533165, 0.82042494859492798, 0.88163884673748105, 0.8197767145135566, 0.82381530984204132, 0.75690954773869346]
[0.24344262295081967, 0.40451745379876797, 0.39930555555555558, 0.46613545816733065, 0.078260869565217397, 0.37947494033412887, 0.30030959752321984, 0.437696335078534, 0.13007159904534607, 0.24320987654320989]
[0.24344262295081967, 0.40451745379876797, 0.44444444444444442, 0.60557768924302791, 0.1072463768115942, 0.3937947494033413, 0.30959752321981426, 0.46492146596858641, 0.15393794749403342, 0.27654320987654318]






BarPos = [[0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4, 2.7, 2.9, 3.2, 3.4],
          [0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4],
          [0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4],
          [0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4]]

Colors = [['black', 'lightgrey'] * 5, 
          ['black', 'lightgrey'] * 5,
          ['black', 'lightgrey'] * 5,
          ['black', 'lightgrey'] * 5]


YLabels = ['% with orthologous gene pairs',
           '',
           '',
           '% excess of orthologous gene pairs\nwith conserved topology']


XTickpos = [[0.4, 0.9, 1.4, 1.9, 2.4, 2.9, 3.4],
            [0.4, 0.9, 1.4, 1.9, 2.4],
            [0.4, 0.9, 1.4, 1.9, 2.4],
            [0.4, 0.9, 1.4, 1.9, 2.4]]

XTicklabels = [['< 0', '0-1000', '1000-10000', '10000-50000', '50000-100000', '100000-150000', '> 150000'],
               [],
               [],
               ['all', 'nst', 'pbk', 'conv', 'div']]

YTicksRange = [np.arange(0, 120, 20),
               np.arange(0, 120, 20),
               np.arange(0, 120, 20),
               np.arange(-20, 120, 20)]

YLimits = [[0, 100], [0, 100], [0, 100], [-20, 100]]



ax1 = CreateAx(2, 2, 1, fig, NonOverlap, Colors[0], BarPos[0], YLabels[0],
               XTickpos[0], XTicklabels[0], YTicksRange[0], YLimits[0][0], YLimits[0][1]) 







ax4 = CreateAx(2, 2, 4, fig, Differences, Colors[3], BarPos[3], YLabels[3],
               XTickpos[3], XTicklabels[3], YTicksRange[3], YLimits[3][0], YLimits[3][1])


# add legend
mouse = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'mouse')
chimp = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'chimp')
ax4.legend(handles = [chimp, mouse], loc = (0.2, 0.8), fontsize = 7, frameon = False, ncol = 2)

# save figure to file
fig.savefig('truc.pdf', bbox_inches = 'tight')




############################



## create figure
#fig = plt.figure(1, figsize = (3, 2))
#
## add a plot to figure (N row, N column, plot N)
#ax = fig.add_subplot(1, 1, 1)
## plot data
#Colors = ['black', 'lightgrey'] * 5    
#BarPos = [0.2, 0.4, 0.7, 0.9, 1.2, 1.4, 1.7, 1.9, 2.2, 2.4]
#ax.bar(BarPos, Differences, width = 0.2, color = Colors, edgecolor = 'black', linewidth = 0.7)
#ax.plot([0, 2.8], [0, 0], lw = 0.7, color = 'black')
#
#
## set font for all text in figure
#FigFont = {'fontname':'Arial'}   
## write y axis label
#ax.set_ylabel('% excess of orthologous gene pairs\nwith conserved topology', color = 'black',  size = 7, ha = 'center', **FigFont)
## add ticks and lebels
#XTickpos = [0.4, 0.9, 1.4, 1.9, 2.4]    
#XTicklabels = ['all', 'nst', 'pbk', 'conv', 'div']
#plt.xticks(XTickpos, XTicklabels, rotation = 0, size = 7, color = 'black', ha = 'center', **FigFont)
## edit y axis ticks
#plt.yticks(np.arange(-20, 120, 20))   
## add a range for the Y and X axes
#plt.ylim([-20, 100])    
## do not show lines around figure  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(False)    
#ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(True)  
#
## make sure the y axis crosses the x axis at 0
#ax.spines['left'].set_position('zero')
#
## edit tick parameters    
#plt.tick_params(axis='both', which='both', bottom='off', top='off',
#                right = 'off', left = 'on', labelbottom='on',
#                colors = 'black', labelsize = 7, direction = 'out')  
## Set the tick labels font name
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')   
#
## add legend
#mouse = mpatches.Patch(facecolor = 'lightgrey' , edgecolor = 'black', linewidth = 0.7, label= 'mouse')
#chimp = mpatches.Patch(facecolor = 'black' , edgecolor = 'black', linewidth = 0.7, label= 'chimp')
#ax.legend(handles = [chimp, mouse], loc = (0.2, 0.8), fontsize = 7, frameon = False, ncol = 2)
#
## save figure to file
#fig.savefig('truc.pdf', bbox_inches = 'tight')
