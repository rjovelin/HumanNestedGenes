# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 23:42:40 2017

@author: Richard
"""



# use this script to download files from Ensembl

# Note that web links were correct as of March 2017

# import modules
import os

# download GFF files for all species of interest 
# create a list of web links to GFF files
GFFLinks = ['ftp://ftp.ensembl.org/pub/release-88/gff3/homo_sapiens/Homo_sapiens.GRCh38.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/mus_musculus/Mus_musculus.GRCm38.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/felis_catus/Felis_catus.Felis_catus_6.2.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/bos_taurus/Bos_taurus.UMD3.1.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/canis_familiaris/Canis_familiaris.CanFam3.1.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/gorilla_gorilla/Gorilla_gorilla.gorGor3.1.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/erinaceus_europaeus/Erinaceus_europaeus.HEDGEHOG.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/equus_caballus/Equus_caballus.EquCab2.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/macaca_mulatta/Macaca_mulatta.Mmul_8.0.1.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/callithrix_jacchus/Callithrix_jacchus.C_jacchus3.2.1.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/monodelphis_domestica/Monodelphis_domestica.BROADO5.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/pongo_abelii/Pongo_abelii.PPYG2.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/ornithorhynchus_anatinus/Ornithorhynchus_anatinus.OANA5.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/choloepus_hoffmanni/Choloepus_hoffmanni.choHof1.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/sorex_araneus/Sorex_araneus.COMMON_SHREW1.88.gff3.gz',
'ftp://ftp.ensembl.org/pub/release-88/gff3/dasypus_novemcinctus/Dasypus_novemcinctus.Dasnov3.0.88.gff3.gz']

# loop of urls and download the corresponding file
for link in GFFLinks:
    os.system('wget ' + link)
    
# download CDS files for human, chimp and mouse
# create a link of web links to CDS files
CDSLinks = ['ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz',
'ftp://ftp.ensembl.org/pub/release-88/fasta/pan_troglodytes/cds/Pan_troglodytes.CHIMP2.1.4.cds.all.fa.gz',
'ftp://ftp.ensembl.org/pub/release-88/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz']

for link in CDSLinks:
    os.system('wget ' + link)


