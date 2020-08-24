#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:16:43 2020

@author: baptiste
"""


import os
import requests

from Bio import SeqIO
from Bio import SeqFeature as sf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

##VARIABLES
#Renseigner ici les chemins vers les fichiers décrits


cloudfile = "fasta_gene_annotation/cloud_seq_methanosarcina_complete/cloud_seq_methanosarcina_complete.emapper.annotations"                         #fichier eggnog annotant la partition cloud
persistentfile = "fasta_gene_annotation/persistent_seq_methanosarcina_complete/persistent_seq_methanosarcina_complete.emapper.annotations"          #fichier eggnog annotant la partition persistent
shellfile = "fasta_gene_annotation/shell_seq_methanosarcina_complete/shell_seq_methanosarcina_complete.emapper.annotations"                         #fichier eggnog annotant la partition shell
allgenesfile = "all_genes/all_genes.fna"                                                                                                            #fichier de sortie de la commande write --all genes de PPanGGOLiN
infoprot= "all_prot_families/representative_gene_families.faa"                                                                                      #fichier de sortie de la commande PPanGGOLiN write --all_prot_families

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



def Eggnog_Genbank(cloud_file, persistent_file, shell_file, all_genes_file, info_prot, species_name):

    #On classe les séquences listées dans le fichier fasta en 3 strings, selon qu'on retrouve leur identifiant dans le fichier cloud, persistent ou shell

    cloudstring = ''
    persistentstring = ''
    shellstring = ''

    #On crée une liste de tuples à partir des résultats d'eggnog contenant les identifiants des gènes et les EC qui leurs sont associés

    with open(persistent_file) as f:
        persistent_ecidlist_avec_headers = [tuple(map(str, i.split('\t'))) for i in f]
        persistent_ecidlist = []
        for k in persistent_ecidlist_avec_headers[4:-3]:
            persistent_ecidlist.append([k[0],k[7]])

    with open(cloud_file) as f:
        cloud_ecidlist_avec_headers = [tuple(map(str, i.split('\t'))) for i in f]
        cloud_ecidlist = []
        for k in cloud_ecidlist_avec_headers[4:-3]:
            cloud_ecidlist.append([k[0],k[7]])

    with open(shell_file) as f:
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
    
    #On importe les séquences protéiques correspondant aux traductions dans un dictionaire
        
    dictionaire_trad = {}

    for i in SeqIO.parse(allgenesfile, 'fasta'):

        dictionaire_trad[i.id] = i.seq
        

    for seq_record in SeqIO.parse(all_genes_file, "fasta"):
        trouve = False
        for idec in cloud_ecidlist:
            if (seq_record.id in idec):
                cloudstring+=str(seq_record.seq)
                cloud_eclist.append([cloud_len, cloud_len + len(str(seq_record.seq)), idec[1], dictionaire_trad[idec[0]], idec[0]])
                cloud_len += len(str(seq_record.seq))
                trouve = True
            
        if not trouve:
            for idec in shell_ecidlist:
                if (seq_record.id in idec):
                    shellstring+=str(seq_record.seq)
                    shell_eclist.append([shell_len, shell_len + len(str(seq_record.seq)), idec[1], dictionaire_trad[idec[0]], idec[0]])
                    shell_len += len(str(seq_record.seq))
                    trouve = True
                
            if not trouve:
                for idec in persistent_ecidlist:
                    if (seq_record.id in idec):
                        persistentstring+=str(seq_record.seq)
                        persistent_eclist.append([persistent_len, persistent_len + len(str(seq_record.seq)), idec[1], dictionaire_trad[idec[0]], idec[0]])
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

            
    num_cloud = 0
    num_persistent = 0
    num_shell = 0
            
    for i in cloud_eclist:
        cloud_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1), type = 'CDS')
        if i[2] != '':
            eclist = i[2].split(',')
            cloud_feature_cds.qualifiers['EC_number'] = eclist
        cloud_feature_cds.qualifiers['old_locus_tag'] = i[4]
        
        tag_part1=i[4][19:27]
    
        cloud_feature_cds.qualifiers['locus_tag'] = "cloud_" + tag_part1 + str(num_cloud)
        num_cloud +=1
    
        cloud_feature_cds.qualifiers['translation'] = i[3]
        cloud_record.features.append(cloud_feature_cds)

    for i in persistent_eclist:
        persistent_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1) , type = 'CDS')
        if i[2] != '':
            eclist = i[2].split(',')
            persistent_feature_cds.qualifiers['EC_number'] = eclist
        persistent_feature_cds.qualifiers['old_locus_tag'] = i[4]
    
        tag_part1=i[4][19:27]
    
        persistent_feature_cds.qualifiers['locus_tag'] = "persistent_" + tag_part1 + str(num_persistent)
        num_persistent +=1
    
        persistent_feature_cds.qualifiers['translation'] = i[3]
        persistent_record.features.append(persistent_feature_cds)
        
    for i in shell_eclist:
        shell_feature_cds = sf.SeqFeature(sf.FeatureLocation(i[0], i[1], +1), type = 'CDS')
        if i[2] != '':
            eclist = i[2].split(',')
            shell_feature_cds.qualifiers['EC_number'] = eclist
        shell_feature_cds.qualifiers['old_locus_tag'] = i[4]
        
        tag_part1=i[4][19:27]
    
    
        shell_feature_cds.qualifiers['locus_tag'] = "shell_" + tag_part1 + str(num_shell)
        num_shell +=1
    
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

    path = os.path.join(os.getcwd(), "GenBank_data_eggnog")
    os.mkdir(path)
    os.chdir(path)
        
    SeqIO.write(cloud_record, "cloud_genbank.gbk", "genbank")
    SeqIO.write(persistent_record, "persistent_genbank.gbk", "genbank")
    SeqIO.write(shell_record, "shell_genbank.gbk", "genbank")

    
## SCRIPT

Eggnog_Genbank(cloudfile, persistentfile, shellfile, allgenesfile, infoprot, speciesname)

