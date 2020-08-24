#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 11:47:33 2020

@author: baptiste
"""

#Ce script doit être placé dans le dossier de sortie des résultats de ppanggolin. Les fichiers txt contenant les noms des séquences cloud / persistent / shell doivent se trouver dans un dossier "partitions", et le fasta "all_genes.fna" obtenu par la commande write avec l'option --all_genes doit se trouver dans un autre dossier intitulé "all_genes"

import os
import requests

from Bio import SeqIO
from Bio import SeqFeature as sf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


##VARIABLES
#Définir ici les chemins vers chaque fichier utilisé

cloudfile = "partitions/cloud.txt"                              #fichier PPanGGOLiN contenant la partition cloud
persistentfile = "partitions/persistent.txt"                    #fichier PPanGGOLiN contenant la partition persistent
shellfile = "partitions/shell.txt"                              #fichier PPanGGOLiN contenant la partition shell
deepecfile = "deepec_output/DeepEC_Result.txt"                  #fichier contenant les annotations de DeepEC
allgenesfile = "all_genes/all_genes.fna"                        #fichier de sortie de la commande PPanGGOLiN write --allgenes
infoprot= "all_prot_families/representative_gene_families.faa"  #fichier de sortie de la commande PPanGGOLiN write --all_prot_families

#Définir ici le nom de l'espèce étudiée
speciesname = "Methanosarcina"



##FONCTION

def create_taxonomic_data(species_name):
    """
    Query the EBI with the species name to create a dictionary containing taxon id,
    taxonomy and some other informations.
    """
    species_informations = {}
    species_name_url = species_name.replace(' ', '%20')

    url = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/' + species_name_url
    response = requests.get(url)
    temp_species_informations = response.json()[0]

    for temp_species_information in temp_species_informations:
        if temp_species_information == 'lineage':
            species_informations['taxonomy'] = temp_species_informations[temp_species_information].split('; ')[:-1]
        elif temp_species_information == 'division':
            species_informations['data_file_division'] = temp_species_informations[temp_species_information]
        elif temp_species_information == 'taxId':
            species_informations['db_xref'] = 'taxon:' + str(temp_species_informations[temp_species_information])
        else:
            species_informations[temp_species_information] = temp_species_informations[temp_species_information]

    compatible_species_name = species_name.replace('/', '_')
    species_informations['description'] = compatible_species_name + ' genome'
    species_informations['organism'] = compatible_species_name
    species_informations['keywords'] = [compatible_species_name]

    return species_informations



def Deepec_genbank(cloud_file, persistent_file, shell_file, deepec_file, all_genes_file, info_prot, species_name)

    #Ouvrir chaque fichier txt contenant les listes de séquences d'intérêt

    cloud = open(cloud_file, "r")
    persistent = open(persistent_file, "r")
    shell = open(shell_file, "r")

    #Importer la liste de numeros d'EC sous forme d'une liste de tuples

    with open(deepec_file) as f:
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
    
    #On importe les séquences protéiques correspondant aux traductions dans un dictionaire
        
    dictionaire_trad = {}

    for i in SeqIO.parse('all_prot_families/representative_gene_families.faa', 'fasta'):
        dictionaire_trad[i.id] = i.seq

    for i in range(0, len(read_info_prot), 2):
        nom_gene = (read_info_prot[i].replace(">","")).replace("\n", "")
        dictionaire_trad[nom_gene] = read_info_prot[i+1].replace("\n", "")
        
        

    for seq_record in SeqIO.parse(all_genes_file, "fasta"):
        if (seq_record.id + "\n" in cloudlines):
            cloudstring+=str(seq_record.seq)
            found_ec = [ec for (identifiant, ec) in ec_list if identifiant==seq_record.id]
            cloud_eclist.append([cloud_len + 1, cloud_len + len(str(seq_record.seq)), found_ec, dictionaire_trad[seq_record.id], seq_record.id])
            cloud_len += len(str(seq_record.seq))
           
        elif (seq_record.id + "\n" in persistentlines):
            persistentstring+=str(seq_record.seq)
            found_ec = [ec for (identifiant, ec) in ec_list if identifiant==seq_record.id]
            persistent_eclist.append([persistent_len + 1, persistent_len + len(str(seq_record.seq)), found_ec, dictionaire_trad[seq_record.id], seq_record.id])
            persistent_len += len(str(seq_record.seq))
    
        elif (seq_record.id + "\n" in shelllines):
            shellstring+=str(seq_record.seq)
            found_ec = [ec for (identifiant, ec) in ec_list if identifiant==seq_record.id]
            shell_eclist.append([shell_len + 1, shell_len + len(str(seq_record.seq)), found_ec, dictionaire_trad[seq_record.id], seq_record.id])
            shell_len += len(str(seq_record.seq))

       
    cloudlist = Seq(cloudstring, generic_dna)
    persistentlist= Seq(persistentstring, generic_dna)
    shelllist = Seq(shellstring, generic_dna)

       
    #On crée les objets SeqRecord correspondants

    cloud_record = SeqRecord(cloudlist, id="Methanosarcina_Cloud_Pseudocontig", name="cloud", description="Pseudo_contig for the methanosarcina cloud genome")
    persistent_record = SeqRecord(persistentlist, id="Methanosarcina_Persistent_Pseudocontig", name="persistent", description="Pseudo_contig for the methanosarcina persistent genome")
    shell_record = SeqRecord(shelllist, id="Methanosarcina_Shell_Pseudocontig", name="shell", description="Pseudo_contig for the methanosarcina shell genome")

    #On crée les Features contenant les EC correspondant à chaque partition, puis on les ajoute aux SeqRecords

    num_cloud = 0
    num_persistent = 0
    num_shell = 0
            
    for i in cloud_eclist:
        cloud_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1), type = 'CDS')
        if i[2] != []:
            cloud_feature_cds.qualifiers['EC_number'] = [(ec.replace('ec:', '')).replace('\n','') for ec in i[2]]
        cloud_feature_cds.qualifiers['old_locus_tag'] = i[4]
        
        tag_part1=i[4][19:27]
    
    
        cloud_feature_cds.qualifiers['locus_tag'] = "cloud_" + tag_part1 + str(num_cloud)
        num_cloud += 1
    
        cloud_feature_cds.qualifiers['translation'] = i[3]
        cloud_record.features.append(cloud_feature_cds)

    for i in persistent_eclist:
        persistent_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1) , type = 'CDS')
        if i[2] != []:
            persistent_feature_cds.qualifiers['EC_number'] = [(ec.replace('ec:', '')).replace('\n','') for ec in i[2]]
        persistent_feature_cds.qualifiers['old_locus_tag'] = i[4]
    
        tag_part1=i[4][19:27]

    
        persistent_feature_cds.qualifiers['locus_tag'] = "persistent_" + tag_part1 + str(num_persistent)
        num_persistent += 1

        persistent_feature_cds.qualifiers['translation'] = i[3]
        persistent_record.features.append(persistent_feature_cds)
        
    for i in shell_eclist:
        shell_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1), type = 'CDS')
        if i[2] != []:
            shell_feature_cds.qualifiers['EC_number'] = [(ec.replace('ec:', '')).replace('\n','') for ec in i[2]]
        shell_feature_cds.qualifiers['old_locus_tag'] = i[4]
        
        tag_part1=i[4][19:27]
    
        shell_feature_cds.qualifiers['locus_tag'] = "shell_" + tag_part1 + str(num_shell)
        num_shell += 1
    
        shell_feature_cds.qualifiers['translation'] = i[3]
        shell_record.features.append(shell_feature_cds)
        
   #On crée les features information


    species_information = create_taxonomic_data(species_name)

    new_feature_source_cloud = sf.SeqFeature(sf.FeatureLocation(1-1, len(cloudstring)),type="source")
    new_feature_source_shell = sf.SeqFeature(sf.FeatureLocation(1-1, len(shellstring)),type="source")
    new_feature_source_persistent = sf.SeqFeature(sf.FeatureLocation(1-1, len(persistentstring)),type="source")


    new_feature_source_cloud.qualifiers['db_xref'] = species_information['db_xref']
    new_feature_source_shell.qualifiers['db_xref'] = species_information['db_xref']
    new_feature_source_persistent.qualifiers['db_xref'] = species_information['db_xref']


    cloud_record.features.append(new_feature_source_cloud)
    shell_record.features.append(new_feature_source_shell)
    persistent_record.features.append(new_feature_source_persistent)
        


    #On crée les fichiers GenBank dans un dossier GenBank_data

    path = os.path.join(os.getcwd(), "GenBank_data_Deepec")
    os.mkdir(path)
    os.chdir(path)
        
    SeqIO.write(cloud_record, "cloud_genbank.gbk", "genbank")
    SeqIO.write(persistent_record, "persistent_genbank.gbk", "genbank")
    SeqIO.write(shell_record, "shell_genbank.gbk", "genbank")


    
## SCRIPT PRINCIPAL



Deepec_genbank(cloudfile, persistentfile, shellfile, deepecfile, allgenesfile, infoprot, speciesname)
