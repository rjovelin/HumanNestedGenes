# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:07:23 2017

@author: RJovelin
"""

#  use this script to plot conservation of overlapping pairs in mammals

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


# 1) plot the proportion of gene pairs in each overlapping gene category that are conserved in other species

# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Macaque', 'Orangutan', 'Marmoset',
           'Mouse', 'Cow', 'Dog', 'Sloth', 'Armadillo', 'Horse', 'Hedgehog',
           'Cat', 'Opossum', 'Platypus', 'Shrew']
# make a list of json files
JsonFiles = ['Overlapping', 'Nested', 'PiggyBack', 'Convergent', 'Divergent']

# make a list of lists of dictionaries for each type of overlapping gene in each species
AllOverlap = []
# loop over files
for i in range(len(Species)):
    ovlp = []
    for j in range(len(JsonFiles)):
        # load dictionary from json file
        json_data = open(Species[i] + JsonFiles[j] + 'Genes.json')
        ovlp.append(json.load(json_data))
        json_data.close()
    AllOverlap.append(ovlp)

    
# make a parallel list of ortholog files
OrthoFiles = ['Human' + i + 'Orthologs.txt' for i in Species[1:]]
## make a list of dicts with orthologs
#AllOrthologs = []
#for i in range(len(OrthoFiles)):
#    orthologs = MatchOrthologs(OrthoFiles[i])
#    AllOrthologs.append(orthologs)

# make a list of GFF files
GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
             'Gorilla_gorilla.gorGor3.1.88.gff3', 'Macaca_mulatta.Mmul_8.0.1.88.gff3',
             'Pongo_abelii.PPYG2.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
             'Mus_musculus.GRCm38.88.gff3', 'Bos_taurus.UMD3.1.88.gff3',
             'Canis_familiaris.CanFam3.1.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
             'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Equus_caballus.EquCab2.88.gff3',
             'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Felis_catus.Felis_catus_6.2.88.gff3',
             'Monodelphis_domestica.BROADO5.88.gff3', 'Ornithorhynchus_anatinus.OANA5.88.gff3',
             'Sorex_araneus.COMMON_SHREW1.88.gff3']

# make a list of gene coordinates       
AllCoordinates, AllOrdered = [], []
# loop over GFF files
for i in range(len(GFF_Files)):
    # get the coordinates of genes on each chromo
    # {chromo: {gene:[chromosome, start, end, sense]}}
    GeneChromoCoord = ChromoGenesCoord(GFF_Files[i])
    # map each gene to its mRNA transcripts
    MapGeneTranscript = GeneToTranscripts(GFF_Files[i])
    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
    AllCoordinates.append(GeneCoord)
    AllOrdered.append(OrderedGenes)

# make a list of lists of pairs of genes
AllGenePairs = []
for i in range(len(AllOverlap)):
    sppairs = []
    for j in range(len(AllOverlap[i])):
        sppairs.append(GetHostNestedPairs(AllOverlap[i][j]))
    AllGenePairs.append(sppairs)

# extract human coordinates, ordered genes, and overlapping gene pairs
HumanCoord, HumanOrdered, HumanGenePairs, HumanAllOverlap = AllCoordinates.pop(0), AllOrdered.pop(0), AllGenePairs.pop(0), AllOverlap.pop(0)
 
assert len(AllOverlap) == len(AllGenePairs) == len(OrthoFiles) == len(Species[1:])

# make a list with counts of conserved pairs for each class of overlapping genes in each species
ConservedAcrossSpecies = []
# loop over species
for i in range(len(AllGenePairs)):
    # make a list with coutns of conserved pairs for each type of overlapping genes
    ConservedPairs = [0] * len(AllGenePairs[i])
    # get the orthologs for that species
    Orthos = MatchOrthologs(OrthoFiles[i])
    # reverse dictionary
    SpeciesOrthos = {}
    for gene in Orthos:
        for ortho in Orthos[gene]:
            if ortho not in SpeciesOrthos:
                SpeciesOrthos[ortho] = [gene]
            else:
                SpeciesOrthos[ortho].append(gene)
    # loop over overlapping gene class
    for j in range(len(AllGenePairs[i])):
        # make pairs of human genes
        HumanPairs = GetHostNestedPairs(HsaAllOverlap[j])
        # remove pairs if any gene is lacking an ortholog
        to_remove = [L for L in HumanPairs if L[0] not in Orthos or L[1] not in Orthos]
        for L in to_remove:
            HumanPairs.remove(L)
        # make pairs of human orthologs for each species gene pairs
        PairsOrthos = []            
        for pair in AllGenePairs[i][j]:
            # check that both genes have corodinates
            assert pair[0] in AllCoordinates[i] and pair[1] in AllCoordinates
            # count only pairs in which both genes have orthos in human
            if pair[0] in SpeciesOrthos and pair[1] in SpeciesOrthos:
                for ortho1 in SpeciesOrthos[pair[0]]:
                    for ortho2 in SpeciesOrthos[pair[1]]:
                        # remove order
                        PairsOrthos.append(set([ortho1, ortho2]))
        # check if human pairs are conserved
        for pair in HumanPairs:
            if set(pair) in PairsOrthos:
                # update counter for the given gene class
                ConservedPairs[j] += 1
        # compute proportion
        ConservedPairs[j] = ConservedPairs[j] / len(humanPairs)
    # populate l;ist for the given species
    ConservedAcrossSpecies.append(ConservedPairs)
        
for i in range(len(ConservedAcrossSpecies)):
    print(i, Species[i+1], ConservedAcrossSpecies[i])
    
    
    
## convert list to numpy array
#Conserved = np.array(Conserved) 
    
## transpose array to get gene categories as columns and species as rows
#Conserved = np.transpose(Conserved)



    
## create a function to format the subplots
#def CreateAx(NColumns, NRows, Grid1, Grid2, RowPos,ColPos, figure, gs, Data, GraphType):
#    '''
#    Returns a ax instance in figure
#    '''    
#    # add a plot using gridspec
#    ax = plt.subplot2grid((Grid1,Grid2), (RowPos,ColPos), colspan=NColumns, rowspan=NRows)
#
#    if GraphType == 'heatmap':
#        # plot heatmap (use vmin and vmax to get the full range of values)
#        heatmap = ax.imshow(Data, interpolation = 'nearest', cmap = 'YlGn')
#        # add heatmap scale 
#        cbar = plt.colorbar(heatmap)
#        # edit tcik parameters of the heatmap scale
#        cbar.ax.tick_params(labelsize=7)
#        cbar.ax.tick_params(direction = 'out')
#        # edit xticks
#        plt.xticks([0,1,2,3,4,5,6], ['Not', 'Nst', 'Ext', 'Int', 'Pbk', 'Con', 'Div'])
#        plt.yticks([i for i in range(16)], ['Chimp', 'Gorilla', 'Orangutan', 'Macaque',
#                    'Marmoset', 'Hedgehog', 'Shrew', 'Cat', 'Dog', 'Mouse', 'Cow', 'Horse',
#                    'Sloth', 'Armadillo','Opossum', 'Platypus'])
#        # do not show lines around figure  
#        ax.spines["top"].set_visible(False)    
#        ax.spines["bottom"].set_visible(False)    
#        ax.spines["right"].set_visible(False)
#        ax.spines["left"].set_visible(False)  
#        # edit tick parameters    
#        plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                        right = 'off', left = 'on', labelbottom='on', labelleft = 'on',
#                        colors = 'black', labelsize = 7, direction = 'out')  
#        # Set the tick labels font name
#        for label in ax.get_yticklabels():
#            label.set_fontname('Arial')   
#    
#    elif GraphType == 'pairs':
#        a, b, c = [i[0] for i in Data], [i[1] for i in Data], [i[2] for i in Data]
#        # make a list for added values for a and b
#        d = [a[i] + b[i] for i in range(len(a))]
#        ## Create a bar plot for proportions of conserved gene pairs
#        ax.bar([0, 0.3, 0.6, 0.9], a, width = 0.2, label = '2 conserved', color= '#9e9ac8', linewidth = 0.7)
#        ax.bar([0, 0.3, 0.6, 0.9], b, width = 0.2, bottom = a, label = '1 conserved', color= '#fd8d3c', linewidth = 0.7)
#        ax.bar([0, 0.3, 0.6, 0.9], c, width = 0.2, bottom = d, label = '0 conserved', color= '#78c679', linewidth = 0.7)
#        LabelSize = 7
#        # set font for all text in figure
#        FigFont = {'fontname':'Arial'}   
#        # write label for y and x axis
#        ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = LabelSize, ha = 'center', **FigFont)
#        # write label for x axis
#        plt.xticks([0.1, 0.4, 0.7, 1], ['Nst', 'Pbk', 'Con', 'Div'], ha = 'center', fontsize = LabelSize, **FigFont)
#        # limit the y axis value range
#        plt.ylim([0, 1])   
#        # do not show lines around figure  
#        ax.spines["top"].set_visible(False)    
#        ax.spines["bottom"].set_visible(True)    
#        ax.spines["right"].set_visible(False)    
#        ax.spines["left"].set_visible(True)  
#        # do not show ticks
#        plt.tick_params(axis='both', which='both', bottom='on', top='off', right = 'off',
#                        left = 'on', labelbottom='on', colors = 'black', labelsize = LabelSize, direction = 'out')  
#        # Set the tick labels font name
#        for label in ax.get_yticklabels():
#            label.set_fontname('Arial')   
#        # create a margin around the x axis
#        plt.margins(0.1)
#        # add legend
#        Two = mpatches.Patch(facecolor = '#9e9ac8' , edgecolor = 'black', linewidth = 0.7, label= '2')
#        One = mpatches.Patch(facecolor = '#fd8d3c' , edgecolor = 'black', linewidth = 0.7, label= '1')
#        Zero = mpatches.Patch(facecolor = '#78c679' , edgecolor = 'black', linewidth = 0.7, label= '0')
#        ax.legend(handles = [Two, One, Zero], loc = (0.1, 1.05), fontsize = LabelSize, frameon = False, ncol = 3)
#    elif GraphType == 'distance':
#        # make a list of gene category names parallel to the list of gene pairs
#        GeneCats = ['Nst', 'Pbk', 'Con', 'Div', 'Prx', 'Mod', 'Int', 'Dst']
#        GeneTypes = ['Nested', 'PiggyBack', 'Convergent', 'Divergent', 'Proximal', 'Moderate', 'Intermediate', 'Distant']
#        # set colors
#        colorscheme = ['#f03b20', '#43a2ca', '#fee391', '#74c476', 'lightgrey', 'lightgrey', 'lightgrey', 'lightgrey']
#        # plot proportions of gene pairs
#        ax.bar([0.05, 0.35, 0.65, 0.95, 1.25, 1.55, 1.85, 2.15], [Data[i] for i in GeneTypes], 0.2, color = colorscheme,
#               edgecolor = 'black', linewidth = 0.7)
#        # set font for all text in figure
#        FigFont = {'fontname':'Arial'}   
#        # write y axis label
#        ax.set_ylabel('Proportion of gene pairs', color = 'black',  size = 7, ha = 'center', **FigFont)
#        # add ticks and lebels
#        plt.xticks([0.15, 0.45, 0.75, 1.05, 1.35, 1.65, 1.95, 2.25], GeneCats, size = 7, color = 'black', ha = 'center', **FigFont)
#
#        ## add a range for the Y and X axes
#        plt.ylim([0, 1])
#        # edit y axis ticks
#        plt.yticks(np.arange(0, 1.2, 0.2)) 
#        plt.xlim([0, 2.45])
#        # do not show lines around figure  
#        ax.spines["top"].set_visible(False)    
#        ax.spines["bottom"].set_visible(True)    
#        ax.spines["right"].set_visible(False)
#        ax.spines["left"].set_visible(True)  
#        # edit tick parameters    
#        plt.tick_params(axis='both', which='both', bottom='on', top='off',
#                        right = 'off', left = 'on', labelbottom='on',
#                        colors = 'black', labelsize = 7, direction = 'out')  
#        # Set the tick labels font name
#        for label in ax.get_yticklabels():
#            label.set_fontname('Arial')   
#      
#    return ax
#
## create figure
#figure = plt.figure(1, figsize = (6.5, 4))
## set up grid
#gs = gridspec.GridSpec(2, 3, width_ratios=[3, 1]) 
## plot data
#ax1 = CreateAx(2, 2, 2, 3, 0, 0, figure, gs, Conserved, 'heatmap')
#ax2 = CreateAx(1, 1, 2, 3, 0, 2, figure, gs, PairCounts, 'pairs')
#ax3 = CreateAx(1, 1, 2, 3, 1, 2, figure, gs, ConservedPairs, 'distance')
#
## add subplot labels
#ax2.text(-2.8, 1.1, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 8)
#ax2.text(-0.4, 1.1, 'B', horizontalalignment='center', verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#ax3.text(-0.4, 1.1, 'C', horizontalalignment='center', verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#
## make sure subplots do not overlap
#plt.tight_layout()
## save figure
#figure.savefig('truc.pdf', bbox_inches = 'tight')
