#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:16:43 2020

@author: baptiste
"""
#Ce script doit être placé dans le dossier de sortie des résultats de ppanggolin. Les fichiers txt contenant les noms des séquences cloud / persistent / shell doivent se trouver dans un dossier "partitions", et le fasta "all_genes.fna" obtenu par la commande write avec l'option --all_genes doit se trouver dans un autre dossier intitulé "all_genes"

import os

from Bio import SeqIO
from Bio import SeqFeature as sf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna




#On classe les séquences listées dans le fichier fasta en 3 strings, selon qu'on retrouve leur identifiant dans le fichier cloud, persistent ou shell

cloudstring = ''
persistentstring = ''
shellstring = ''

#On crée une liste de tuples à partir des résultats d'eggnog contenant les identifiants des gènes et les EC qui leurs sont associés

with open('fasta_gene_annotation/persistent_seq_methanosarcina_complete/persistent_seq_methanosarcina_complete.emapper.annotations') as f:
    persistent_ecidlist_avec_headers = [tuple(map(str, i.split('\t'))) for i in f]
    persistent_ecidlist = []
    for k in persistent_ecidlist_avec_headers[4:-3]:
        persistent_ecidlist.append([k[0],k[7]])

with open('fasta_gene_annotation/cloud_seq_methanosarcina_complete/cloud_seq_methanosarcina_complete.emapper.annotations') as f:
    cloud_ecidlist_avec_headers = [tuple(map(str, i.split('\t'))) for i in f]
    cloud_ecidlist = []
    for k in cloud_ecidlist_avec_headers[4:-3]:
        cloud_ecidlist.append([k[0],k[7]])

with open('fasta_gene_annotation/shell_seq_methanosarcina_complete/shell_seq_methanosarcina_complete.emapper.annotations') as f:
    shell_ecidlist_avec_headers = [tuple(map(str, i.split('\t'))) for i in f]
    shell_ecidlist = []
    for k in shell_ecidlist_avec_headers[4:-3]:
        shell_ecidlist.append([k[0],k[7]])
        

#On liste les EC et les séquences, ainsi que leur position dans le pseudo-contig, en fonction de leur appartenance au cloud, persistent ou shell. On crée ainsi un pseudo-contig  pour chaque partition, et on répertorie les EC et les positions de chaque séquence dans des listes
        
cloud_eclist=[]
shell_eclist=[]
persistent_eclist=[]        

cloud_len = 0
persistent_len = 0
shell_len = 0

for seq_record in SeqIO.parse("all_genes/all_genes.fna", "fasta"):
    trouve = False
    for idec in cloud_ecidlist:
        if (seq_record.id in idec):
            cloudstring+=str(seq_record.seq)
            cloud_eclist.append([cloud_len, cloud_len + len(str(seq_record.seq)), idec[1]])
            cloud_len += len(str(seq_record.seq))
            trouve = True
            
    if not trouve:
        for idec in shell_ecidlist:
            if (seq_record.id in idec):
                shellstring+=str(seq_record.seq)
                shell_eclist.append([shell_len, shell_len + len(str(seq_record.seq)), idec[1]])
                shell_len += len(str(seq_record.seq))
                trouve = True
                
        if not trouve:
            for idec in persistent_ecidlist:
                if (seq_record.id in idec):
                    persistentstring+=str(seq_record.seq)
                    persistent_eclist.append([persistent_len, persistent_len + len(str(seq_record.seq)), idec[1]])
                    persistent_len += len(str(seq_record.seq))
                    trouve = True

 
cloudlist = Seq(cloudstring, generic_dna)
persistentlist= Seq(persistentstring, generic_dna)
shelllist = Seq(shellstring, generic_dna)


#On crée les objets SeqRecord correspondants

cloud_record = SeqRecord(cloudlist, id="Methanosarcina_Cloud_Pseudocontig", name="cloud", description="Pseudo_contig for the methanosarcina cloud genome")
persistent_record = SeqRecord(persistentlist, id="Methanosarcina_Persistent_Pseudocontig", name="persistent", description="Pseudo_contig for the methanosarcina persistent genome")
shell_record = SeqRecord(shelllist, id="Methanosarcina_Shell_Pseudocontig", name="shell", description="Pseudo_contig for the methanosarcina shell genome")

#On crée les Features contenant les EC correspondant à chaque partition, puis on les ajoute aux SeqRecords


            
for i in cloud_eclist:
    cloud_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1), type = 'CDS')
    if i[2] != '':
        eclist = i[2].split(',')
        cloud_feature_cds.qualifiers['EC_number'] = eclist
    cloud_record.features.append(cloud_feature_cds)

for i in persistent_eclist:
    persistent_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1) , type = 'CDS')
    if i[2] != '':
        eclist = i[2].split(',')
        persistent_feature_cds.qualifiers['EC_number'] = eclist
    persistent_record.features.append(persistent_feature_cds)
        
for i in shell_eclist:
    shell_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1), type = 'CDS')
    if i[2] != '':
        eclist = i[2].split(',')
        shell_feature_cds.qualifiers['EC_number'] = eclist
    shell_record.features.append(shell_feature_cds)
        


#On crée les fichiers GenBank dans un dossier GenBank_data

path = os.path.join(os.getcwd(), "GenBank_data/eggnog")
os.mkdir(path)
os.chdir(path)
        
SeqIO.write(cloud_record, "cloud_genbank.gbk", "genbank")
SeqIO.write(persistent_record, "persistent_genbank.gbk", "genbank")
SeqIO.write(shell_record, "shell_genbank.gbk", "genbank")
