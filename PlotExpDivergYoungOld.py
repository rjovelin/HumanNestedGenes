# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 11:38:39 2016

@author: RJovelin
"""


# use this script to plot expression divergence between external their
# un-nested orthologs and between internal and their un-nested orthologs

# usage python3 PlotExpDivergYoungOld.py 



########################################################################


## import modules
## use Agg backend on server without X server
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#from matplotlib import rc
#rc('mathtext', default='regular')
#import json
#import random
#import copy
#import sys
#import os
#import math
#import numpy as np
#from scipy import stats
#from HsaNestedGenes import *
#
#
## load dictionaries of overlapping genes
#jsonFiles = ['HumanOverlappingGenes.json', 'HumanNestedGenes.json',  
#             'ChimpOverlappingGenes.json', 'ChimpNestedGenes.json',
#             'GorillaOverlappingGenes.json', 'GorillaNestedGenes.json']
## make a list of dictionaries
#AllOverlap = []
## loop over files
#for i in range(len(jsonFiles)):
#    # load dictionary of overlapping gene pairs
#    json_data = open(jsonFiles[i])
#    overlapping = json.load(json_data)
#    json_data.close()
#    AllOverlap.append(overlapping)
#
#
#
#
#
#
## get GFF file
#GFF = ['Homo_sapiens.GRCh38.86.gff3', 'Pan_troglodytes.CHIMP2.1.4.86.gff3', 'Gorilla_gorilla.gorGor3.1.86.gff3']
#
## make a list of gene coordinates       
#AllCoordinates, AllOrdered = [], []
## loop over GFF files
#for i in range(len(GFF)):
#    # get the coordinates of genes on each chromo
#    # {chromo: {gene:[chromosome, start, end, sense]}}
#    GeneChromoCoord = ChromoGenesCoord(GFF[i])
#    # map each gene to its mRNA transcripts
#    MapGeneTranscript = GeneToTranscripts(GFF[i])
#    # remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#    GeneChromoCoord = FilterOutGenesWithoutValidTranscript(GeneChromoCoord, MapGeneTranscript)
#    # get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#    GeneCoord = FromChromoCoordToGeneCoord(GeneChromoCoord)
#    # Order genes along chromo {chromo: [gene1, gene2, gene3...]} 
#    OrderedGenes = OrderGenesAlongChromo(GeneChromoCoord)
#    AllCoordinates.append(GeneCoord)
#    AllOrdered.append(OrderedGenes)
#
#
## make pairs of overlapping genes
#AllPairs = []
#for i in range(len(AllOverlap)):
#    pairs = GetHostNestedPairs(AllOverlap[i])
#    print(len(pairs))
#    AllPairs.append(pairs)
## make pairs of overlapping genes
#HumanPairs = AllPairs[:2]
#ChimpPairs = AllPairs[2:4]
#GorillaPairs = AllPairs[4:]
#
#
#
#
#
#inthuman, exthuman, intchimp, extchimp, intgorilla, extgorilla = [], [] , [] ,[], [], []
#for pair in HumanPairs[1]:
#    exthuman.append(pair[0])
#    inthuman.append(pair[1])
#for pair in ChimpPairs[1]:
#    extchimp.append(pair[0])
#    intchimp.append(pair[1])
#for pair in GorillaPairs[1]:
#    extgorilla.append(pair[0])
#    intgorilla.append(pair[1])
#
#
#
#
#
## make list with sets of non-overlapping genes
#NonOverlappingSets = []
#for i in range(3):
#    j = i * 2
#    print(i, j, jsonFiles[j])
#    # make a set of non-overlapping gene
#    nonoverlap = MakeNonOverlappingGeneSet(AllOverlap[j], AllCoordinates[i])
#    NonOverlappingSets.append(nonoverlap)    
#
## make sets of host and nested nested genes
#NestedSets = []
#for i in range(1, len(AllOverlap), 2):
#    nestedset = MakeFullPartialOverlapGeneSet(AllOverlap[i])
#    NestedSets.append(nestedset)
#
## get 1:1 orthologs between human, chimp and gorilla {human:[chimp,gorilla]}
#Orthos = MatchOrthologTrios('HumanChimpGorillaOrthologs.txt')
#
## make pairs of overlapping genes
#machin = []
#for i in range(len(AllOverlap)):
#    pairs = GetHostNestedPairs(AllOverlap[i])
#    machin.append(pairs)
## make pairs of overlapping genes
#Humanmachin = machin[1]
#print(len(Humanmachin))
#to_remove = []
#for pair in Humanmachin:
#    if pair[0] not in Orthos and pair[1] not in Orthos:
#        to_remove.append(pair)
#for pair in to_remove:
#    Humanmachin.remove(pair)
#print(len(Humanmachin))
#
#assert 0 > 1
#
#
#
#
#
#
#
#print(len(HumanPairs[0]), len(HumanPairs[1]), len(ChimpPairs[0]), len(ChimpPairs[1]), len(GorillaPairs[0]), len(GorillaPairs[1]))
#
## make sets of chimp and gorilla genes  with orthologs
#chimporthos, gorillaorthos = set(), set()
#for gene in Orthos:
#    chimporthos.add(Orthos[gene][0])
#    gorillaorthos.add(Orthos[gene][1])
#
## remove gene pairs lacking orthos
#for i in range(len(HumanPairs)):
#    to_remove = []
#    for pair in HumanPairs[i]:
#        if pair[0] not in Orthos or pair[1] not in Orthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        HumanPairs[i].remove(pair)
## remove chimp pairs lacking orthologs
#for i in range(len(ChimpPairs)):
#    to_remove = []
#    for pair in ChimpPairs[i]:
#        if pair[0] not in chimporthos or pair[1] not in chimporthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        ChimpPairs[i].remove(pair)
## remove gorilla pairs lacking orthologs
#for i in range(len(GorillaPairs)):
#    to_remove = []
#    for pair in GorillaPairs[i]:
#        if pair[0] not in gorillaorthos or pair[1] not in gorillaorthos:
#            to_remove.append(pair)
#    for pair in to_remove:
#        GorillaPairs[i].remove(pair)
#
#print(len(HumanPairs[0]), len(HumanPairs[1]), len(ChimpPairs[0]), len(ChimpPairs[1]), len(GorillaPairs[0]), len(GorillaPairs[1]))
#
#
## create a dict with chimp genes as key {chimp_gene: [human_orto, gorilla_ortho]}
#ChimpOrthologs = {}
#for gene in Orthos:
#    ChimpOrthologs[Orthos[gene][0]] = [gene, Orthos[gene][1]]
#
#
## infer young and old nesting events in human and chimp
#HumanOld, HumanYoung = InferYoungOldNestingEvents(Orthos, ChimpPairs[1], GorillaPairs[1], HumanPairs[1])
#ChimpOld, ChimpYoung = InferYoungOldNestingEvents(ChimpOrthologs, HumanPairs[1], GorillaPairs[1], ChimpPairs[1])
#
#print(len(HumanOld), len(HumanYoung))
#print(len(ChimpOld), len(ChimpYoung))
#
#
## do QC
#for pair in HumanOld:
#    ortho1, ortho2 = Orthos[pair[0]][0], Orthos[pair[1]][0]
#    assert [ortho1, ortho2] in ChimpOld
#for pair in ChimpOld:
#    ortho1, ortho2 = ChimpOrthologs[pair[0]][0], ChimpOrthologs[pair[1]][0]
#    assert [ortho1, ortho2] in HumanOld
#
#    
## get expression profile in human and chimp
#HumanExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'human')
#ChimpExpression = ParsePrimateExpressionData('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt', 'chimp')
## remove genes without expression
#HumanExpression = RemoveGenesLackingExpression(HumanExpression)
#ChimpExpression = RemoveGenesLackingExpression(ChimpExpression)
## get relative expression
#HumanExpression = TransformRelativeExpression(HumanExpression)
#ChimpExpression = TransformRelativeExpression(ChimpExpression)
#
#
## remove gene pairs with genes lacking expression
#HumanOld = FilterGenePairsWithoutExpression(HumanOld, HumanExpression)
#HumanYoung = FilterGenePairsWithoutExpression(HumanYoung, HumanExpression)
#ChimpOld = FilterGenePairsWithoutExpression(ChimpOld, ChimpExpression)
#ChimpYoung = FilterGenePairsWithoutExpression(ChimpYoung, ChimpExpression)
#
#InferredPairs = [HumanOld, HumanYoung, ChimpOld, ChimpYoung]
#
#
#
#
#for i in InferredPairs:
#    print(len(i))
#
## compare expression divergence between young nested genes and their un-nested orthologs
## compare expression divergence between young host genes and their un-nested orthologs
#
#
## make sets of external and internal genes
#ExtIntGenes = []
#for i in range(len(InferredPairs)):
#    external, internal = set(), set()
#    for pair in InferredPairs[i]:
#        external.add(pair[0])
#        internal.add(pair[1])
#    ExtIntGenes.append(external)
#    ExtIntGenes.append(internal)
#
#a = []
#for i in ExtIntGenes:
#    a.append(len(i))
#
##[hsaextold, hsaintold, hsaextyoung, hsaintyoung, ptrextold, ptrintold, ptrextyoung, ptrextyoung]
#
#c = []
#
#
## for young external and internal, remove genes if ortholog is nested
#for i in range(len(ExtIntGenes)):
#    if i in [2, 3, 6, 7]:
#        to_remove = []
#        for gene in ExtIntGenes[i]:
#            if i in [2, 3]:
#                if Orthos[gene][0] in NestedSets[1] or Orthos[gene][1] in NestedSets[2]:
#                    # gene is nested in other species, cannot be a young nested gene in human
#                    to_remove.append(gene)
#                    c.append(gene)
#                    
#                    
#                    
#            elif i in [6, 7]:
#                if ChimpOrthologs[gene][0] in NestedSets[0] or ChimpOrthologs[gene][1] in NestedSets[2]:
#                    # gene is nested in other species, cannot be a young nested gene in chimp
#                    to_remove.append(gene)
#        for gene in to_remove:
#            ExtIntGenes[i].remove(gene)
#    else:
#        for gene in ExtIntGenes[i]:
#            if i in [0, 1]:
#                assert Orthos[gene][0] in NestedSets[1] or Orthos[gene][1] in NestedSets[2]
#            elif i in [4, 5]:
#                assert ChimpOrthologs[gene][0] in NestedSets[0] or ChimpOrthologs[gene][1] in NestedSets[2]
#
#
#
#
#
#
#
#
#
#b = []
#for i in ExtIntGenes:
#    b.append(len(i))
#
#
#
#
#for i in range(len(a)):
#    print(i, a[i], b[i], a[i] == b[i])
#
#hsaextptr, hsaintptr, hsaextgo, hsaintgo = [], [] , [] , []
#for gene in c:
#    if Orthos[gene][0] in intchimp:
#        hsaintptr.append(gene)
#    if Orthos[gene][0] in extchimp:
#        hsaextptr.append(gene)
#    if Orthos[gene][1] in intgorilla:
#        hsaintgo.append(gene)
#    if Orthos[gene][1] in extgorilla:
#        hsaextgo.append(gene)
#
#print(set(inthuman).intersection(set(hsaintptr)))
#print(set(inthuman).intersection(set(hsaextptr)))
#
#print(set(exthuman).intersection(set(hsaintptr)))
#print(set(exthuman).intersection(set(hsaextptr)))
#
#
#
## make pairs of overlapping genes
#truc = []
#for i in range(len(AllOverlap)):
#    pairs = GetHostNestedPairs(AllOverlap[i])
#    truc.append(pairs)
## make pairs of overlapping genes
#Humantruc = truc[1]
#Chimptruc = truc[3]
#Gorillatruc = truc[5]
#
#
#
#
#for gene in {'ENSG00000197757', 'ENSG00000145936', 'ENSG00000152254', 'ENSG00000187918', 'ENSG00000204010'}:
#    for pair in Humantruc:
#        if gene == pair[1]:
#            for doublet in Chimptruc:
#                if Orthos[gene][0] == doublet[1] and pair[0] in Orthos and pair[1] in Orthos:
#                    print(pair, [Orthos[pair[0]][0], Orthos[pair[1]][0]], doublet)
#
#
#
#
#{'ENSG00000119917'}
#{'ENSG00000197757', 'ENSG00000152253', 'ENSG00000159917'}
#
#
#
##[hsaextold, hsaintold, hsaextyoung, hsaintyoung, ptrextold, ptrintold, ptrextyoung, ptrextyoung]
#for gene in ExtIntGenes[2]:
#    assert Orthos[gene][0] not in ExtIntGenes[6]
#    assert Orthos[gene][0] not in ExtIntGenes[7]
#for gene in ExtIntGenes[3]:
#    assert Orthos[gene][0] not in ExtIntGenes[6]
#    assert Orthos[gene][0] not in ExtIntGenes[7]
#for gene in ExtIntGenes[0]:
#    assert Orthos[gene][0] in ExtIntGenes[4]
#    assert Orthos[gene][0] not in ExtIntGenes[5]
#for gene in ExtIntGenes[1]:
#    assert Orthos[gene][0] in ExtIntGenes[5]
#    assert Orthos[gene][0] not in ExtIntGenes[4] 
#
#
#
## make pairs of overlapping genes
#machin = []
#for i in range(len(AllOverlap)):
#    pairs = GetHostNestedPairs(AllOverlap[i])
#    machin.append(pairs)
## make pairs of overlapping genes
#Humanmachin = machin[1]
#print(len(Humanmachin))
#to_remove = []
#for pair in Humanmachin:
#    if pair[0] not in Orthos and pair[1] not in Orthos:
#        to_remove.append(pair)
#for pair in to_remove:
#    Humanmachin.remove(pair)
#print(len(Humanmachin))


##############################################################





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



extmissing, intmissing = set(), set()
for pair in HumanPairs[1]:
    if pair[0] not in Orthos:
        extmissing.add(pair[0])
    if pair[1] not in Orthos:
        intmissing.add(pair[1])
print(len(extmissing), len(intmissing))


print(len(HumanPairs[0]), len(HumanPairs[1]))
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
print(len(HumanPairs[0]), len(HumanPairs[1]))


SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG, NstCons = [], [], [], [], [], []
Conserved, NotConserved = 0, 0


# loop over human nested gene pairs
for pair in HumanPairs[1]:
    # get the orthologs of the human genes
    ortho1, ortho2 = Orthos[pair[0]], Orthos[pair[1]]
    # check if human gene pairs is nested in chimp
    if set([ortho1, ortho2]) not in ChimpPairs[1]:
        NotConserved += 1
        # get chromosomes of chimp orthologs
        chromo1, chromo2 = ChimpCoord[ortho1][0], ChimpCoord[ortho2][0]        
        if chromo1 == chromo2:
            # get start positions of chimp orthologs
            S1, S2 = ChimpCoord[ortho1][1], ChimpCoord[ortho2][1]
            # get end positions of chimp orthologs
            E1, E2 = ChimpCoord[ortho1][2], ChimpCoord[ortho2][2]            
            # get indices of chimp orthologs in the ordred gene list
            P1, P2 = ChimpOrdered[chromo1].index(ortho1), ChimpOrdered[chromo2].index(ortho2)
            # check if orthologs are overlapping in chimp
            if set([ortho1, ortho2]) in ChimpPairs[0]:
                # chimp genes are overlapping
                # check if they are adjcent
                if S1 < S2:
                    assert P1 < P2
                    if P2 != P1 + 1:
                        # overlapping non adjacent
                        SepOvlp.append(pair)
                    elif P2 == P1 + 1:
                        # overlapping adjacent
                        AdjOvlp.append(pair)
                elif S1 > S2:
                    assert P2 < P1
                    if P1 != P2 + 1:
                        # overlapping non adjacent
                        SepOvlp.append(pair)
                    elif P1 == P2 + 1:
                        # overlapping adjacent
                        AdjOvlp.append(pair)
                elif S1 == S2:
                    print('merde')
            elif set([ortho1, ortho2]) not in ChimpPairs[0]:
                # chimp genes are not overlapping
                # check if they are adjacent
                if S1 < S2:
                    assert P1 < P2
                    if P2 != P1 + 1:
                        # non-overlapping non-adjacent
                        SepNonOvlp.append(pair)
                    elif P2 == P1 + 1:
                        # non-overlapping adjacent
                        AdjNonOvlp.append(pair)
                elif S1 > S2:
                    assert P2 < P1
                    if P1 != P2 + 1:
                        # non-overlapping non adjacent
                        SepNonOvlp.append(pair)
                    elif P1 == P2 + 1:
                        # non-overlapping adjacent
                        AdjNonOvlp.append(pair)
                elif S1 == S2:
                    print('shit')
              
    elif set([ortho1, ortho2]) in ChimpPairs[1]:
        # orthologs of human nested genes are also nested
        Conserved += 1            

# make a list of counts
Pairs = [SepOvlp, SepNonOvlp, AdjOvlp, AdjNonOvlp, DiffLG, NstCons]

# make a list of labels
Labels = ['Separated overlapping', 'Separated non-overlapping', 'Adjacent overlapping', 'Adjacent non-overlapping', 'Different chromosomes', 'Nested']
# remove counts and labels equal to 0
to_remove = []
for i in range(len(Pairs)):
    if len(Pairs[i]) != 0:
        print(Labels[i])
        print(Pairs[i][0][0], Pairs[i][0][1], Orthos[Pairs[i][0][0]], Orthos[Pairs[i][0][1]])
        print(HumanCoord[Pairs[i][0][0]], HumanCoord[Pairs[i][0][1]], ChimpCoord[Orthos[Pairs[i][0][0]]], ChimpCoord[Orthos[Pairs[i][0][1]]])
        
