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

with open('ChimpContainedGenes.json') as chimp_json_data:
    ChimpContained = json.load(chimp_json_data)
with open('MouseContainedGenes.json') as mouse_json_data:
    MouseContained = json.load(mouse_json_data)
with open('MacaqueContainedGenes.json') as macaque_json_data:
    MacaqueContained = json.load(macaque_json_data)

# get GFF file
HsaGFF = 'Homo_sapiens.GRCh38.86.gff3'
MmuGFF = 'Mus_musculus.GRCm38.86.gff3'    
CfaGFF = 'Canis_familiaris.CanFam3.1.86.gff3'
PtrGFF = 'Pan_troglodytes.CHIMP2.1.4.86.gff3'
Mca = 'Macaca_mulatta.Mmul_8.0.1.86.gff3' 

# make lists of host-nested gene pairs
HumanHostNestedPairs = GetHostNestedPairs(HumanHostGenes)
print('host-gene pairs in human', len(HumanHostNestedPairs))
MouseHostNestedPairs = GetHostNestedPairs(MouseHostGenes)
print('host-gene pairs in mouse', len(MouseHostNestedPairs))
DogHostNestedPairs = GetHostNestedPairs(DogHostGenes)
print('host-gene pairs in dog', len(DogHostNestedPairs))
ChimpHostNestedPairs = GetHostNestedPairs(ChimpHostGenes)
print('host-gene pairs in chimp', len(ChimpHostNestedPairs))
MacaqueHostnestedPairs = GetHostNestedPairs(MacaqueHostGenes)
print('host-gene pairs in macaque', len(MacaqueHostNestedPairs))

# get the 1:1 orthologs between human and mouse
HumanMouseOrthos = ParseOrthologFile('HumanMouseOrthologs.txt')
print('mapped human genes to their orthologs in mouse', len(HumanMouseOrthos))
HumanChimpOrthos = ParseOrthologFile('HumanChimpOrthologs.txt')
ChimpMouseOrthos = ParseOrthologFile('ChimpMouseOrthologs.txt')
HumanDogOrthos = ParseOrthologFile('HumanDogOrthologs.txt')
MouseDogOrthos = ParseOrthologFile('MouseDogOrthologs.txt')
ChimpMacaqueOrthos = ParseOrthologFile('ChimpMacaqueOrthologs.txt')
HumanMacaqueOrthos = ParseOrthologFile('HumanMacaqueOrthologs.txt')



HsaOrthoPtr, HsaOrthoMmu, HsaOrthoCfa, PtrOrthoMmu, HsaOrthoMca, PtrOrthoMca = 0, 0, 0, 0, 0, 0
HsaMmuPtrOrtho, HsaMmuCfaOrtho, HsaMcaPtrOrtho = 0, 0, 0

for pair in HumanHostNestedPairs:
    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos:
        HsaOrthoMmu += 1
    if pair[0] in HumanDogOrthos and pair[1] in HumanDogOrthos:
        HsaOrthoCfa += 1
    if pair[0] in HumanChimpOrthos and pair[1] in HumanChimpOrthos:
        HsaOrthoPtr += 1
    if pair[0] in HumanMacaqueOrthos and pair[1] in HumanMacaqueOrthos:
        HsaOrthoMca += 1
    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos and pair[0] in HumanChimpOrthos and pair[1] in HumanChimpOrthos:
        HsaMmuPtrOrtho += 1
    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos and pair[0] in HumanDogOrthos and pair[1] in HumanDogOrthos:
        HsaMmuCfaOrtho += 1
    if pair[0] in HumanMacaqueOrthos and pair[1] in HumanMacaqueOrthos and pair[0] in HumanChimpOrthos and pair[1] in HumanChimpOrthos:
        HsaMcaPtrOrtho += 1
    

print('ortho human-mouse', HsaOrthoMmu)
print('ortho human-dog', HsaOrthoCfa)
print('ortho human-chimp', HsaOrthoPtr)
print('ortho human-mouse-dog', HsaMmuCfaOrtho)
print('ortho human-mouse-chimp', HsaMmuPtrOrtho)
print('ortho human-chimp-macaque', HsaMcaPtrOrtho)

PtrOrthoMmu = 0
for pair in ChimpHostNestedPairs:
    if pair[0] in ChimpMouseOrthos and pair[1] in ChimpMouseOrthos:
        PtrOrthoMmu += 1
    if pair[0] in ChimpMacaqueOrthos and pair[1] in ChimpMacaqueOrthos:
        PtrOrthoMca += 1
print('orthos chimp-mouse', PtrOrthoMmu)
print('orthos chimp-macaque', PtrOrthoMca)

# make set of host and nested genes
chimphostcontained = MakeHostNestedGeneSet(ChimpContained)
mousehostcontained = MakeHostNestedGeneSet(MouseContained)
chimphostnested = MakeHostNestedGeneSet(ChimpHostGenes)
mousehostnested = MakeHostNestedGeneSet(MouseHostGenes)


hsaspecific = 0
hsaspnotcontained = 0
hsaspecontained = 0

# get the human specific nested pairs
for pair in HumanHostNestedPairs:
    if pair[0] in HumanMouseOrthos and pair[1] in HumanMouseOrthos and pair[0] in HumanChimpOrthos and pair[1] in HumanChimpOrthos:
        if HumanMouseOrthos[pair[0]][0] not in mousehostnested and HumanMouseOrthos[pair[1]][0] not in mousehostnested:
            if HumanChimpOrthos[pair[0]][0] not in chimphostnested and HumanChimpOrthos[pair[1]][0] not in chimphostnested:
                hsaspecific += 1
                if (HumanMouseOrthos[pair[0]][0] in mousehostcontained and HumanMouseOrthos[pair[1]][0] in mousehostcontained) or (HumanChimpOrthos[pair[0]][0] in chimphostcontained and HumanChimpOrthos[pair[1]][0] in chimphostcontained):
                    hsaspecontained += 1
                elif HumanMouseOrthos[pair[0]][0] not in mousehostcontained and HumanMouseOrthos[pair[1]][0] not in mousehostcontained and HumanChimpOrthos[pair[0]][0] not in chimphostcontained and HumanChimpOrthos[pair[1]][0] not in chimphostcontained:
                    hsaspnotcontained += 1
print('human specific', hsaspecific)
print('human specific not contained', hsaspnotcontained)
print('human specific contained', hsaspecontained)            
            
                   
        
