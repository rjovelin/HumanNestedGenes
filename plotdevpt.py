# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 18:58:58 2017

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


# make a parallel list of Species names
Species = ['Human', 'Chimp', 'Gorilla', 'Orangutan',
           'Macaque', 'Marmoset', 'Hedgehog', 'Shrew',
           'Cat', 'Dog', 'Mouse', 'Cow', 'Horse', 'Sloth',
           'Armadillo', 'Opossum', 'Platypus']

# make a list of json files
JsonFiles = ['Nested', 'PiggyBack', 'Convergent', 'Divergent']

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

# make a list of GFF files
GFF_Files = ['Homo_sapiens.GRCh38.88.gff3', 'Pan_troglodytes.CHIMP2.1.4.88.gff3',
             'Gorilla_gorilla.gorGor3.1.88.gff3', 'Pongo_abelii.PPYG2.88.gff3',
             'Macaca_mulatta.Mmul_8.0.1.88.gff3', 'Callithrix_jacchus.C_jacchus3.2.1.88.gff3',
             'Erinaceus_europaeus.HEDGEHOG.88.gff3', 'Sorex_araneus.COMMON_SHREW1.88.gff3',
             'Felis_catus.Felis_catus_6.2.88.gff3', 'Canis_familiaris.CanFam3.1.88.gff3',
             'Mus_musculus.GRCm38.88.gff3', 'Bos_taurus.UMD3.1.88.gff3',
             'Equus_caballus.EquCab2.88.gff3', 'Choloepus_hoffmanni.choHof1.88.gff3',
             'Dasypus_novemcinctus.Dasnov3.0.88.gff3', 'Monodelphis_domestica.BROADO5.88.gff3',
             'Ornithorhynchus_anatinus.OANA5.88.gff3']

# make a list of gene coordinates       
AllCoordinates = []
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
    AllCoordinates.append(GeneCoord)
    
# make a list of lists of pairs of genes
AllGenePairs = []
for i in range(len(AllOverlap)):
    sppairs = []
    for j in range(len(AllOverlap[i])):
        sppairs.append(GetHostNestedPairs(AllOverlap[i][j]))
    AllGenePairs.append(sppairs)

# extract human coordinates, ordered genes, and overlapping gene pairs
HumanCoord, HumanGenePairs, HumanAllOverlap = AllCoordinates.pop(0), AllGenePairs.pop(0), AllOverlap.pop(0)
 
assert len(AllOverlap) == len(AllGenePairs) == len(OrthoFiles) == len(Species[1:])


# 2) plot the proportion of nested gene pairs conserved in each species when human gene pairs have same or oppositte strand orientation

# make a list with counts of conserved pairs in each species for nested pairs separretely for same and possite orientation
ConservedNested = []
HumanNested = [copy.deepcopy(HumanAllOverlap[1]), copy.deepcopy(HumanAllOverlap[1])]
# make a list of nested pairs for each species
SpeciesNested = []
for i in range(len(AllGenePairs)):
    SpeciesNested.append([copy.deepcopy(AllGenePairs[i][1]), copy.deepcopy(AllGenePairs[i][1])])

# loop over species
for i in range(len(SpeciesNested)):
    # make a list with coutns of conserved pairs for each type of overlapping genes
    ConservedPairs = [0] * len(SpeciesNested[i])
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
    for j in range(len(SpeciesNested[i])):
        # make pairs of human genes
        HumanPairs = GetHostNestedPairs(HumanNested[j])
        # remove pairs if any gene is lacking an ortholog
        to_remove = [L for L in HumanPairs if L[0] not in Orthos or L[1] not in Orthos]
        for L in to_remove:
            HumanPairs.remove(L)
        if j == 0:
            # remove human nested pairs if genes have opposite orientation
            to_remove = [L for L in HumanPairs if len(set(GenePairOrientation(L, HumanCoord))) == 2]
            for L in to_remove:
                HumanPairs.remove(L)
        elif j == 1:
            # remove human nested pairs if genes have same orientation
            to_remove = [L for L in HumanPairs if len(set(GenePairOrientation(L, HumanCoord))) == 1]
            for L in to_remove:
                HumanPairs.remove(L)
        # make pairs of human orthologs for each species gene pairs
        PairsOrthos = []            
        for pair in SpeciesNested[i][j]:
            # check that both genes have corodinates
            assert pair[0] in AllCoordinates[i] and pair[1] in AllCoordinates[i]
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
        ConservedPairs[j] = ConservedPairs[j] / len(HumanPairs)
    # populate l;ist for the given species
    ConservedNested.append(ConservedPairs)
        
# convert list to numpy array
ConservedNested = np.array(ConservedNested) 



# 3) plot the proportion of gene pairs that are overlapping in each species 

# make a list with counts of conserved pairs in each species for overlapping genes defined by strand orientation
ConservedOrientation = []

# make a list of gene pairs by pooling same and opposite strand overlapping gene pairs [[same], [opposite]]
HumanOrientation = [[], []]
for i in range(1, len(HumanAllOverlap)):
    # make a list of human pairs
    pairs = GetHostNestedPairs(HumanAllOverlap[i])
    if i == 1:
        # check if nested pairs are same or opposite strand
        for L in pairs:
            if len(set(GenePairOrientation(L, HumanCoord))) == 2:
                HumanOrientation[1].append(L)
            elif len(set(GenePairOrientation(L, HumanCoord))) == 1:
                HumanOrientation[0].append(L)
    elif i == 2:
        # add all pbk pairs to same strand list
        for L in pairs:
            HumanOrientation[0].append(L)
    elif i == 3 or i == 4:
        # add convergent or divergent pairs to opposite strand list
        for L in pairs:
            HumanOrientation[1].append(L)

# make a list of all overlapping gene pairs for each species
SpeciesOverlap = []
for i in range(len(AllOverlap)):
    SpeciesOverlap.append(GetHostNestedPairs(AllOverlap[i][0]))
    
# loop over species
for i in range(len(SpeciesOverlap)):
    # make a list with coutns of conserved pairs for each orientation [same, opposite]
    ConservedPairs = [0, 0]
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
    # make pairs of human orthologs for each species gene pairs
    PairsOrthos = []            
    for pair in SpeciesOverlap[i]:
        # check that both genes have corodinates
        assert pair[0] in AllCoordinates[i] and pair[1] in AllCoordinates[i]
        # count only pairs in which both genes have orthos in human
        if pair[0] in SpeciesOrthos and pair[1] in SpeciesOrthos:
            for ortho1 in SpeciesOrthos[pair[0]]:
                for ortho2 in SpeciesOrthos[pair[1]]:
                    # remove order
                    PairsOrthos.append(set([ortho1, ortho2]))
    # loop over lists of human pairs
    for j in range(len(HumanOrientation)):
        # make a copy of the list iof human pairs
        HumanPairs = copy.deepcopy(HumanOrientation[j])
        # remove pairs if any gene is lacking an ortholog
        to_remove = [L for L in HumanPairs if L[0] not in Orthos or L[1] not in Orthos]
        for L in to_remove:
            HumanPairs.remove(L)
        # check if human pairs are conserved
        for pair in HumanPairs:
            if set(pair) in PairsOrthos:
                # update counter for the given gene class
                ConservedPairs[j] += 1
        # compute proportion
        ConservedPairs[j] = ConservedPairs[j] / len(HumanPairs)
    # populate l;ist for the given species
    ConservedOrientation.append(ConservedPairs)
        
# convert list to numpy array
ConservedOrientation = np.array(ConservedOrientation) 
    


# 4) plot the proportion of gene pairs overlapping in each otehr species

ConservedOverlap = []
# loop over species
for i in range(len(SpeciesOverlap)):
    # make a list with coutns of conserved pairs for each type of overlapping genes
    ConservedPairs = [0] * len(HumanGenePairs)
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
    # make pairs of human orthologs for each species gene pairs
    PairsOrthos = []            
    for pair in SpeciesOverlap[i]:
        # check that both genes have corodinates
        assert pair[0] in AllCoordinates[i] and pair[1] in AllCoordinates[i]
        # count only pairs in which both genes have orthos in human
        if pair[0] in SpeciesOrthos and pair[1] in SpeciesOrthos:
            for ortho1 in SpeciesOrthos[pair[0]]:
                for ortho2 in SpeciesOrthos[pair[1]]:
                    # remove order
                    PairsOrthos.append(set([ortho1, ortho2]))
    # loop over lists of human pairs
    for j in range(len(HumanGenePairs)):
        # make a copy of the list iof human pairs
        HumanPairs = copy.deepcopy(HumanGenePairs[j])
        # remove pairs if any gene is lacking an ortholog
        to_remove = [L for L in HumanPairs if L[0] not in Orthos or L[1] not in Orthos]
        for L in to_remove:
            HumanPairs.remove(L)
        # check if human pairs are conserved
        for pair in HumanPairs:
            if set(pair) in PairsOrthos:
                # update counter for the given gene class
                ConservedPairs[j] += 1
        # compute proportion
        ConservedPairs[j] = ConservedPairs[j] / len(HumanPairs)
    # populate l;ist for the given species
    ConservedOverlap.append(ConservedPairs)
        
# convert list to numpy array
ConservedOverlap = np.array(ConservedOverlap) 


























## make a list with counts of conserved pairs in each species for overlapping genes defined by strand orientation
#ConservedOrientation = []
#
## make a list of dicts in human with [nested, nested, pbk, con, div] 
## remove dict of non-overlapping gene
#HumanAllOverlap.remove(HumanAllOverlap[0])
## insert a copy of dict with nested genes
#HumanAllOverlap.insert(0, copy.deepcopy(HumanAllOverlap[0]))
#
## for each species, remove the pairs of non-overlapping genes and make a copy of the nested pairs
#for i in range(len(AllGenePairs)):
#    # remove the pairs of non-overlapping genes
#    AllGenePairs[i].remove(AllGenePairs[i][0])
#    # make a copy of the nested gene pairs
#    AllGenePairs[i].insert(0, copy.deepcopy(AllGenePairs[i][0]))
#
## loop over species
#for i in range(len(AllGenePairs)):
#    # make a list with coutns of conserved pairs for each type of overlapping genes
#    ConservedPairs = [0] * len(AllGenePairs[i])
#    # get the orthologs for that species
#    Orthos = MatchOrthologs(OrthoFiles[i])
#    # reverse dictionary
#    SpeciesOrthos = {}
#    for gene in Orthos:
#        for ortho in Orthos[gene]:
#            if ortho not in SpeciesOrthos:
#                SpeciesOrthos[ortho] = [gene]
#            else:
#                SpeciesOrthos[ortho].append(gene)
#    # loop over overlapping gene class
#    for j in range(len(AllGenePairs[i])):
#        # make pairs of human genes
#        HumanPairs = GetHostNestedPairs(HumanAllOverlap[j])
#        # remove pairs if any gene is lacking an ortholog
#        to_remove = [L for L in HumanPairs if L[0] not in Orthos or L[1] not in Orthos]
#        for L in to_remove:
#            HumanPairs.remove(L)
#        if j == 0:
#            # remove human nested pairs if genes have opposite orientation
#            to_remove = [L for L in HumanPairs if len(set(GenePairOrientation(L, HumanCoord))) == 2]
#            for L in to_remove:
#                HumanPairs.remove(L)
#        elif j == 1:
#            # remove human nested pairs if genes have same orientation
#            to_remove = [L for L in HumanPairs if len(set(GenePairOrientation(L, HumanCoord))) == 1]
#            for L in to_remove:
#                HumanPairs.remove(L)
#        # make pairs of human orthologs for each species gene pairs
#        PairsOrthos = []            
#        for pair in AllGenePairs[i][j]:
#            # check that both genes have corodinates
#            assert pair[0] in AllCoordinates[i] and pair[1] in AllCoordinates[i]
#            # count only pairs in which both genes have orthos in human
#            if pair[0] in SpeciesOrthos and pair[1] in SpeciesOrthos:
#                for ortho1 in SpeciesOrthos[pair[0]]:
#                    for ortho2 in SpeciesOrthos[pair[1]]:
#                        # remove order
#                        PairsOrthos.append(set([ortho1, ortho2]))
#        # check if human pairs are conserved
#        for pair in HumanPairs:
#            if set(pair) in PairsOrthos:
#                # update counter for the given gene class
#                ConservedPairs[j] += 1
#        # compute proportion
#        ConservedPairs[j] = ConservedPairs[j] / len(HumanPairs)
#    # populate l;ist for the given species
#    ConservedOrientation.append(ConservedPairs)
#        
## convert list to numpy array
#ConservedOrientation = np.array(ConservedOrientation) 
    
## transpose array to get gene categories as columns and species as rows
#Conserved = np.transpose(Conserved)

   
# create a function to format the subplots
def CreateAx(NColumns, NRows, Grid1, Grid2, RowPos,ColPos, figure, gs, Data, GraphType):
    '''
    Returns a ax instance in figure
    '''    
    # add a plot using gridspec
    ax = plt.subplot2grid((Grid1,Grid2), (RowPos,ColPos), colspan=NColumns, rowspan=NRows)

    if GraphType == 'pairs':
        # plot heatmap (use vmin and vmax to get the full range of values)
        heatmap = ax.imshow(Data, interpolation = 'nearest', cmap = 'Purples')
        # add heatmap scale 
        cbar = plt.colorbar(heatmap)
        # edit tcik parameters of the heatmap scale
        cbar.ax.tick_params(labelsize=7)
        cbar.ax.tick_params(direction = 'out')
        # edit xticks
        plt.xticks([0,1,2,3,4], ['Not', 'Nst', 'Pbk', 'Con', 'Div'])
    elif GraphType == 'orientation':
        # plot heatmap (use vmin and vmax to get the full range of values)
        heatmap = ax.imshow(Data, interpolation = 'nearest', cmap = 'Oranges')
        # add heatmap scale 
        cbar = plt.colorbar(heatmap)
        # edit tcik parameters of the heatmap scale
        cbar.ax.tick_params(labelsize=7)
        cbar.ax.tick_params(direction = 'out')
        # edit xticks
        plt.xticks([0,1,2,3,4], ['NstS', 'NstO', 'Pbk', 'Con', 'Div'])
   
    elif GraphType == 'nested':
        # plot heatmap (use vmin and vmax to get the full range of values)
        heatmap = ax.imshow(Data, interpolation = 'nearest', cmap = 'Oranges')
        # add heatmap scale 
        cbar = plt.colorbar(heatmap)
        # edit tcik parameters of the heatmap scale
        cbar.ax.tick_params(labelsize=7)
        cbar.ax.tick_params(direction = 'out')
        # edit xticks
        plt.xticks([0,1], ['NstS', 'NstO'])
   
    plt.yticks([i for i in range(16)], ['Chimp', 'Gorilla', 'Orangutan', 'Macaque',
               'Marmoset', 'Hedgehog', 'Shrew', 'Cat', 'Dog', 'Mouse', 'Cow', 'Horse',
               'Sloth', 'Armadillo','Opossum', 'Platypus'])
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)  
    # edit tick parameters    
    plt.tick_params(axis='both', which='both', bottom='on', top='off',
                    right = 'off', left = 'on', labelbottom='on', labelleft = 'on',
                    colors = 'black', labelsize = 7, direction = 'out')  
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')   
    return ax

# create figure
figure = plt.figure(1, figsize = (6.5, 4))
# set up grid
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
# plot data
ax1 = CreateAx(1, 1, 1, 2, 0, 0, figure, gs, ConservedAcrossSpecies, 'pairs')
#ax2 = CreateAx(1, 1, 1, 2, 0, 1, figure, gs, ConservedOrientation, 'orientation')


#ax2 = CreateAx(1, 1, 1, 2, 0, 1, figure, gs, ConservedNested, 'nested')




#ax2 = CreateAx(1, 1, 1, 2, 0, 1, figure, gs, ConservedOrientation, 'nested')

ax2 = CreateAx(1, 1, 1, 2, 0, 1, figure, gs, ConservedOverlap, 'pairs')

## add subplot labels
#ax2.text(-2.8, 1.1, 'A', ha='center', va='center', color = 'black', fontname = 'Arial', size = 8)
#ax2.text(-0.4, 1.1, 'B', horizontalalignment='center', verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)
#ax3.text(-0.4, 1.1, 'C', horizontalalignment='center', verticalalignment='center', color = 'black', fontname = 'Arial', size = 8)

# make sure subplots do not overlap
plt.tight_layout()
# save figure
figure.savefig('truc.pdf', bbox_inches = 'tight')
