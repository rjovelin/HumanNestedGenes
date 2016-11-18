# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 17:19:12 2016

@author: RJovelin
"""




import os
import sys
import numpy as np
import random
import json
import math
from HsaNestedGenes import *




# load dictionaries of host and nested genes 
# gene names have wormbase ID for cel and cbr, but transcript names for cr
with open('HumanHostNestedGenes.json') as human_json_data:
    HumanHostGenes = json.load(human_json_data)
with open('MouseHostNestedGenes.json') as mouse_json_data:
    MouseHostGenes = json.load(mouse_json_data)
with open('DogHostNestedGenes.json') as dog_json_data:
    DogHostGenes = json.load(dog_json_data)
with open('ChimpHostNestedGenes.json') as chimp_json_data:
    ChimpHostGenes = json.load(chimp_json_data)
with open('MacaqueHostNestedGenes.json') as macaque_json_data:
    MacaqueHostGenes = json.load(macaque_json_data)



with open('HumanContainedGenes.json') as human_json_data:
    HumanContained = json.load(human_json_data)
with open('ChimpContainedGenes.json') as chimp_json_data:
    ChimpContained = json.load(chimp_json_data)
with open('MouseContainedGenes.json') as mouse_json_data:
    MouseContained = json.load(mouse_json_data)
with open('MacaqueContainedGenes.json') as macaque_json_data:
    MacaqueContained = json.load(macaque_json_data)
with open('DogContainedGenes.json') as dog_json_data:
    DogContained = json.load(dog_json_data)


# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
MmuGFF = 'Mus_musculus.GRCm38.86.gff3'    
CfaGFF = 'Canis_familiaris.CanFam3.1.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
McaGFF = 'Macaca_mulatta.Mmul_8.0.1.86.gff3' 


# make lists of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
print('host-gene pairs in human', len(HumanHostNestedPairs))
MouseHostNestedPairs = GetHostNestedPairs(MouseHostGenes)
print('host-gene pairs in mouse', len(MouseHostNestedPairs))
DogHostNestedPairs = GetHostNestedPairs(DogHostGenes)
print('host-gene pairs in dog', len(DogHostNestedPairs))
ChimpHostNestedPairs = GetHostNestedPairs(ChimpHostGenes)
print('host-gene pairs in chimp', len(ChimpHostNestedPairs))
MacaqueHostNestedPairs = GetHostNestedPairs(MacaqueHostGenes)
print('host-gene pairs in macaque', len(MacaqueHostNestedPairs))




## make lists of host-nested gene pairs
#HumanHostNestedPairs = GetHostNestedPairs(HumanContained)
#print('host-gene pairs in human', len(HumanHostNestedPairs))
#MouseHostNestedPairs = GetHostNestedPairs(MouseContained)
#print('host-gene pairs in mouse', len(MouseHostNestedPairs))
#DogHostNestedPairs = GetHostNestedPairs(DogContained)
#print('host-gene pairs in dog', len(DogHostNestedPairs))
#ChimpHostNestedPairs = GetHostNestedPairs(ChimpContained)
#print('host-gene pairs in chimp', len(ChimpHostNestedPairs))
#MacaqueHostNestedPairs = GetHostNestedPairs(MacaqueContained)
#print('host-gene pairs in macaque', len(MacaqueHostNestedPairs))


infile = open('NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt')
infile.readline()

orthos = {}
for line in infile:
    if 'ENS' in line:
        line = line.rstrip().split('\t')
        gene1, gene2, gene3 = line[0], line[1], line[4]
        assert gene1 not in orthos
        orthos[gene1] = [gene2, gene3]
infile.close()

chimporthos = set()
macaqueorthos = set()
for i in orthos:
    chimporthos.add(orthos[i][0])
    macaqueorthos.add(orthos[i][1])




hsaorthos, ptrorthos, mmlorthos = 0, 0, 0
for pair in HumanHostNestedPairs:
    if pair[0] in orthos and pair[1] in orthos:
        hsaorthos += 1
for pair in ChimpHostNestedPairs:
    if pair[0] in chimporthos and pair[1] in chimporthos:
        ptrorthos += 1
for pair in MacaqueHostNestedPairs:
    if pair[0] in macaqueorthos and pair[1] in macaqueorthos:
        mmlorthos += 1

print('orthos pair', hsaorthos, ptrorthos, mmlorthos)




orthologs = ParseOrthologFile('HumanChimpMacaqueOrthologs.txt')
print(len(orthos), len(orthologs))

missing = []
different = 0
for gene in orthos:
    if gene in orthologs:
        if orthologs[gene] != orthos[gene]:
            print(gene, orthologs[gene][0], orthos[gene][0])
            print(gene, orthologs[gene][1], orthos[gene][1])            
            break

#
#
#
## get the coordinates of human genes on each chromo
## {chromo: {gene:[chromosome, start, end, sense]}}
#HumanGeneChromoCoord = ChromoGenesCoord(HsaGFF)
## map each gene to its mRNA transcripts
#HumanMapGeneTranscript = GeneToTranscripts(HsaGFF)
## remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#HumanGeneChromoCoord = FilterOutGenesWithoutValidTranscript(HumanGeneChromoCoord, HumanMapGeneTranscript)
## get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#HumanGeneCoord = FromChromoCoordToGeneCoord(HumanGeneChromoCoord)
#
## get the coordinates of human genes on each chromo
## {chromo: {gene:[chromosome, start, end, sense]}}
#ChimpGeneChromoCoord = ChromoGenesCoord(PtrGFF)
## map each gene to its mRNA transcripts
#ChimpMapGeneTranscript = GeneToTranscripts(PtrGFF)
## remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#ChimpGeneChromoCoord = FilterOutGenesWithoutValidTranscript(ChimpGeneChromoCoord, ChimpMapGeneTranscript)
## get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#ChimpGeneCoord = FromChromoCoordToGeneCoord(ChimpGeneChromoCoord)
#
## get the coordinates of human genes on each chromo
## {chromo: {gene:[chromosome, start, end, sense]}}
#MacaqueGeneChromoCoord = ChromoGenesCoord(McaGFF)
## map each gene to its mRNA transcripts
#MacaqueMapGeneTranscript = GeneToTranscripts(McaGFF)
## remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#MacaqueGeneChromoCoord = FilterOutGenesWithoutValidTranscript(MacaqueGeneChromoCoord, MacaqueMapGeneTranscript)
## get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#MacaqueGeneCoord = FromChromoCoordToGeneCoord(MacaqueGeneChromoCoord)
#
## get the coordinates of human genes on each chromo
## {chromo: {gene:[chromosome, start, end, sense]}}
#DogGeneChromoCoord = ChromoGenesCoord(CfaGFF)
## map each gene to its mRNA transcripts
#DogMapGeneTranscript = GeneToTranscripts(CfaGFF)
## remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#DogGeneChromoCoord = FilterOutGenesWithoutValidTranscript(DogGeneChromoCoord, DogMapGeneTranscript)
## get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#DogGeneCoord = FromChromoCoordToGeneCoord(DogGeneChromoCoord)
#
## get the coordinates of human genes on each chromo
## {chromo: {gene:[chromosome, start, end, sense]}}
#MouseGeneChromoCoord = ChromoGenesCoord(MmuGFF)
## map each gene to its mRNA transcripts
#MouseMapGeneTranscript = GeneToTranscripts(MmuGFF)
## remove genes that do not have a mRNA transcripts (may have abberant transcripts, NMD processed transcripts, etc)
#MouseGeneChromoCoord = FilterOutGenesWithoutValidTranscript(MouseGeneChromoCoord, MouseMapGeneTranscript)
## get the coordinates of each gene {gene:[chromosome, start, end, sense]}
#MouseGeneCoord = FromChromoCoordToGeneCoord(MouseGeneChromoCoord)
#
#hsasame, hsaopposite = 0, 0
#for pair in HumanHostNestedPairs:
#    orientation = set(GenePairOrientation(pair, HumanGeneCoord))
#    if len(orientation) == 1:
#        hsasame += 1
#    elif len(orientation) == 2:
#        hsaopposite += 1
#ptrsame, ptropposite = 0, 0
#for pair in ChimpHostNestedPairs:
#    oientation = set(GenePairOrientation(pair, ChimpGeneCoord))
#    if len(orientation) == 1:
#        ptrsame += 1
#    elif len(orientation) == 2:
#        ptropposite += 1
#mmlsame, mmlopposite = 0, 0
#for pair in MacaqueHostNestedPairs:
#    orientation = set(GenePairOrientation(pair, MacaqueGeneCoord))
#    if len(orientation) == 1:
#        mmlsame += 1
#    elif len(orientation) == 2:
#        mmlopposite += 1    
#mmusame, mmuopposite = 0, 0
#for pair in MouseHostNestedPairs:
#    orientation = set(GenePairOrientation(pair, MouseGeneCoord))
#    if len(orientation) == 1:
#        mmusame += 1
#    elif len(orientation) == 2:
#        mmuopposite += 1 
#cfasame, cfaopposite = 0,0
#for pair in DogHostNestedPairs:
#    orientation = set(GenePairOrientation(pair, DogGeneCoord))
#    if len(orientation) == 1:
#        cfasame += 1
#    elif len(orientation) == 2:
#        cfaopposite += 1
#
#
#print('human: same {0}, opposite {1}'.format(hsasame, hsaopposite))
#
#print('chimp: same {0}, opposite {1}'.format(ptrsame, ptropposite))
#
#print('macaque: same {0}, opposite {1}'.format(mmlsame, mmlopposite))
#
#print('mouse: same {0}, opposite {1}'.format(mmusame, mmuopposite))
#
#print('dog: same {0}, opposite {1}'.format(cfasame, cfaopposite))
#
##
### get the 1:1 orthologs between human and mouse
##HumanMouseOrthos = ParseOrthologFile('HumanMouseOrthologs.txt')
##print('mapped human genes to their orthologs in mouse', len(HumanMouseOrthos))
##HumanChimpOrthos = ParseOrthologFile('HumanChimpOrthologs.txt')
##ChimpMouseOrthos = ParseOrthologFile('ChimpMouseOrthologs.txt')
##HumanDogOrthos = ParseOrthologFile('HumanDogOrthologs.txt')
##MouseDogOrthos = ParseOrthologFile('MouseDogOrthologs.txt')
##ChimpMacaqueOrthos = ParseOrthologFile('ChimpMacaqueOrthologs.txt')
##HumanMacaqueOrthos = ParseOrthologFile('HumanMacaqueOrthologs.txt')
##
##
##
##
##
##
##
##
##
##
##HumanChimpMacaqueOrthos = ParseOrthologFile('HumanChimpMacaqueOrthologs.txt')
##print('human-chimp-macaque', len(HumanChimpMacaqueOrthos))
##
##HsaOrthoPtr, HsaOrthoMmu, HsaOrthoCfa, PtrOrthoMmu, HsaOrthoMca, PtrOrthoMca = 0, 0, 0, 0, 0, 0
##HsaMmuPtrOrtho, HsaMmuCfaOrtho, HsaMcaPtrOrtho = 0, 0, 0
##
##HsaPtrMcaOrthos = 0
##
##NestedHomolgy = []
##
##for pair in HumanHostNestedPairs:
##    if pair[0] in HumanChimpMacaqueOrthos and pair[1] in HumanChimpMacaqueOrthos:
##        HsaPtrMcaOrthos += 1
##        NestedHomolgy.append(pair)
##        
##    
##    
###    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos:
###        HsaOrthoMmu += 1
###    if pair[0] in HumanDogOrthos and pair[1] in HumanDogOrthos:
###        HsaOrthoCfa += 1
###    if pair[0] in HumanChimpOrthos and pair[1] in HumanChimpOrthos:
###        HsaOrthoPtr += 1
###    if pair[0] in HumanMacaqueOrthos and pair[1] in HumanMacaqueOrthos:
###        HsaOrthoMca += 1
###    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos and pair[0] in HumanChimpOrthos and pair[1] in HumanChimpOrthos:
###        HsaMmuPtrOrtho += 1
###    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos and pair[0] in HumanDogOrthos and pair[1] in HumanDogOrthos:
###        HsaMmuCfaOrtho += 1
###    if pair[0] in HumanMacaqueOrthos and pair[1] in HumanMacaqueOrthos and pair[0] in HumanChimpOrthos and pair[1] in HumanChimpOrthos:
###        HsaMcaPtrOrtho += 1
##    
##
###print('ortho human-mouse', HsaOrthoMmu)
###print('ortho human-dog', HsaOrthoCfa)
###print('ortho human-chimp', HsaOrthoPtr)
###print('ortho human-mouse-dog', HsaMmuCfaOrtho)
###print('ortho human-mouse-chimp', HsaMmuPtrOrtho)
###print('ortho human-chimp-macaque', HsaMcaPtrOrtho)
##
##print('ortho human-chimp-macaque', HsaPtrMcaOrthos)
##
##
##
##
###PtrOrthoMmu = 0
###for pair in ChimpHostNestedPairs:
###    if pair[0] in ChimpMouseOrthos and pair[1] in ChimpMouseOrthos:
###        PtrOrthoMmu += 1
###    if pair[0] in ChimpMacaqueOrthos and pair[1] in ChimpMacaqueOrthos:
###        PtrOrthoMca += 1
###print('orthos chimp-mouse', PtrOrthoMmu)
###print('orthos chimp-macaque', PtrOrthoMca)
##
### make set of host and nested genes
##chimphostcontained = MakeHostNestedGeneSet(ChimpContained)
##mousehostcontained = MakeHostNestedGeneSet(MouseContained)
##chimphostnested = MakeHostNestedGeneSet(ChimpHostGenes)
##mousehostnested = MakeHostNestedGeneSet(MouseHostGenes)
##macaquehostnested = MakeHostNestedGeneSet(MacaqueHostGenes)
##macaquehostcontained = MakeHostNestedGeneSet(MacaqueContained)
##
##hsaspecific = 0
##hsaspnotcontained = 0
##hsaspecontained = 0
##
### get the human specific nested pairs
##for pair in NestedHomolgy:
##    if HumanChimpMacaqueOrthos[pair[0]][0][0] not in chimphostcontained and HumanChimpMacaqueOrthos[pair[1]][0][0] not in chimphostcontained and HumanChimpMacaqueOrthos[pair[0]][1][0] not in macaquehostcontained and HumanChimpMacaqueOrthos[pair[1]][1][0] not in macaquehostcontained:
##        hsaspnotcontained += 1
##    
##print('human specific not contained', hsaspnotcontained)
##            
##                   
##        
