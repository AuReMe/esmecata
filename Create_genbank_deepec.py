#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:47:33 2020

@author: baptiste
"""

#Ce script doit être placé dans le dossier de sortie des résultats de ppanggolin. Les fichiers txt contenant les noms des séquences cloud / persistent / shell doivent se trouver dans un dossier "partitions", et le fasta "all_genes.fna" obtenu par la commande write avec l'option --all_genes doit se trouver dans un autre dossier intitulé "all_genes"

import os

from Bio import SeqIO
from Bio import SeqFeature as sf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


#Ouvrir chaque fichier txt contenant les listes de séquences d'intérêt

cloud = open("partitions/cloud.txt", "r")
persistent = open("partitions/persistent.txt", "r")
shell = open("partitions/shell.txt", "r")

#Importer la liste de numeros d'EC sous forme d'une liste de tuples

with open('deepec_output/DeepEC_Result.txt') as f:
    ec_list = [tuple(map(str, i.split('\t'))) for i in f]

#Exporter le contenu des txt pour référence

cloudlines = cloud.readlines()
persistentlines = persistent.readlines()
shelllines = shell.readlines()

#On classe les séquences listées dans le fichier fasta en 3 strings, selon qu'on retrouve leur identifiant dans le fichier cloud, persistent ou shell

cloudstring = ''
persistentstring = ''
shellstring = ''

#On initialise aussi les listes des EC associés:

cloud_eclist = []
persistent_eclist = []
shell_eclist = []

#On répertorie les longueurs des séquences

cloud_len = 0
persistent_len = 0
shell_len = 0

for seq_record in SeqIO.parse("all_genes/all_genes.fna", "fasta"):
    if (seq_record.id + "\n" in cloudlines):
        cloudstring+=str(seq_record.seq)
        found_ec = [ec for (identifiant, ec) in ec_list if identifiant==seq_record.id]
        cloud_eclist.append([cloud_len + 1, cloud_len + len(str(seq_record.seq)), found_ec])
        cloud_len += len(str(seq_record.seq))
           
    elif (seq_record.id + "\n" in persistentlines):
        persistentstring+=str(seq_record.seq)
        found_ec = [ec for (identifiant, ec) in ec_list if identifiant==seq_record.id]
        persistent_eclist.append([persistent_len + 1, persistent_len + len(str(seq_record.seq)), found_ec ])
        persistent_len += len(str(seq_record.seq))
    
    elif (seq_record.id + "\n" in shelllines):
        shellstring+=str(seq_record.seq)
        found_ec = [ec for (identifiant, ec) in ec_list if identifiant==seq_record.id]
        shell_eclist.append([shell_len + 1, shell_len + len(str(seq_record.seq)), found_ec])
        shell_len += len(str(seq_record.seq))

       
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
    if i[2] != []:
        cloud_feature_cds.qualifiers['EC_number'] = [(ec.replace('ec:', '')).replace('\n','') for ec in i[2]]
    cloud_record.features.append(cloud_feature_cds)

for i in persistent_eclist:
    persistent_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1) , type = 'CDS')
    if i[2] != []:
        persistent_feature_cds.qualifiers['EC_number'] = [(ec.replace('ec:', '')).replace('\n','') for ec in i[2]]
    persistent_record.features.append(persistent_feature_cds)
        
for i in shell_eclist:
    shell_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1), type = 'CDS')
    if i[2] != []:
        shell_feature_cds.qualifiers['EC_number'] = [(ec.replace('ec:', '')).replace('\n','') for ec in i[2]]
    shell_record.features.append(shell_feature_cds)
        


#On crée les fichiers GenBank dans un dossier GenBank_data

path = os.path.join(os.getcwd(), "GenBank_data")
os.mkdir(path)
os.chdir(path)
        
SeqIO.write(cloud_record, "cloud_genbank.gbk", "genbank")
SeqIO.write(persistent_record, "persistent_genbank.gbk", "genbank")
SeqIO.write(shell_record, "shell_genbank.gbk", "genbank")

