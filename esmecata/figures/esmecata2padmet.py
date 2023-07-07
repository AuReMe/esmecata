import csv
import os.path
import shutil
import sys

import mpwt
from padmet.utils.connection.pgdb_to_padmet import from_pgdb_to_padmet
from padmet.utils.connection.sbmlGenerator import padmet_to_sbml
from padmet.utils.exploration.compare_padmet import compare_padmet
from padmet.classes.padmetRef import PadmetRef


def get_obs_to_keep(tax_diff_file):
    obs_to_keep = dict()
    with open(tax_diff_file, 'r') as f:
        lines = list(csv.reader(f, delimiter='\t'))
        for l in lines[1:]:
            obs = l[-1].split(';')[0]
            name = l[1].lower()
            forbidden_char = ['*', '.', '"', '/', '\\', '[', ']', ':', ';', '|', ' ']
            for char in forbidden_char:
                name = name.replace(char, '_')
            tax_id = l[5]
            value = (name, tax_id)
            if value not in obs_to_keep.values():
                obs_to_keep[obs] = value
    return obs_to_keep


def move_files(run_name, otk):
    dir_pathologic_source = os.path.join(run_name, '2_annotation', 'pathologic')
    network_path = os.path.join(run_name, 'Networks')
    if not os.path.exists(network_path):
        os.mkdir(network_path)
    dir_pathologic_destination = os.path.join(network_path, 'input')
    if not os.path.exists(dir_pathologic_destination):
        os.mkdir(dir_pathologic_destination)

    for obs, value in otk.items():
        name = value[0]
        dir_obs_destination = os.path.join(dir_pathologic_destination, name)
        if not os.path.exists(dir_obs_destination):
            os.mkdir(dir_obs_destination)
        file_pathologic_source = os.path.join(dir_pathologic_source, obs,
                                              os.listdir(os.path.join(dir_pathologic_source, obs))[0])
        file_pathologic_destination = os.path.join(dir_obs_destination, os.listdir(os.path.join(
            dir_pathologic_source, obs))[0].replace(obs, name))
        if not os.path.exists(file_pathologic_destination):
            shutil.copy(file_pathologic_source, file_pathologic_destination)

    taxon_destination = os.path.join(dir_pathologic_destination, 'taxon_id.tsv')
    if not os.path.exists(taxon_destination):
        with open(taxon_destination, 'w') as f:
            f.write('species\ttaxon_id')
            for obs, value in otk.items():
                f.write(f'\n{value[0]}\t{value[1]}')

    return network_path


def generate_input(run_name):
    tax_diff_file = os.path.join(run_name, 'Krona', f'{run_name}_taxonomy_diff.tsv')
    otk = get_obs_to_keep(tax_diff_file)
    move_files(run_name, otk)


def networks_reconstruction(run_name, nb_cpu):
    input_networks = os.path.join(run_name, 'Networks', 'input')
    output_networks = os.path.join(run_name, 'Networks', 'PGDBs')
    if not os.path.exists(output_networks):
        os.mkdir(output_networks)
    taxon_file = os.path.join(input_networks, 'taxon_id.tsv')
    log_path = os.path.join(run_name, 'Networks', 'logs')
    if not os.path.exists(log_path):
        os.mkdir(log_path)

    if os.listdir(output_networks) == list():
        mpwt.multiprocess_pwt(input_folder=input_networks,
                              output_folder=output_networks,
                              patho_inference=True,
                              flat_creation=True,
                              dat_extraction=True,
                              number_cpu=nb_cpu,
                              patho_log=log_path,
                              taxon_file=taxon_file,
                              verbose=True)


def create_padmet(run_name, database_path):
    pgdb_path = os.path.join(run_name, 'Networks', 'PGDBs')
    padmet_path = os.path.join(run_name, 'Networks', 'PADMETs')
    if not os.path.exists(padmet_path):
        os.mkdir(padmet_path)
    for obs in os.listdir(pgdb_path):
        input_pgdb = os.path.join(pgdb_path, obs)
        output_padmet = os.path.join(padmet_path, f'{obs}.padmet')
        if not os.path.exists(output_padmet):
            from_pgdb_to_padmet(input_pgdb, padmetRef_file=database_path, source="genome", extract_gene=True,
                                no_orphan=True, verbose=True, output_file=output_padmet)


def create_sbml(run_name):
    padmet_path = os.path.join(run_name, 'Networks', 'PADMETs')
    sbml_path = os.path.join(run_name, 'Networks', 'SBMLs')
    if not os.path.exists(sbml_path):
        os.mkdir(sbml_path)
    for obs in os.listdir(padmet_path):
        input_padmet = os.path.join(padmet_path, obs)
        output_sbml = os.path.join(sbml_path, obs.replace('.padmet', '.sbml'))
        if not os.path.exists(output_sbml):
            padmet_to_sbml(padmet=input_padmet, output=output_sbml, verbose=True)


def generate_tables(run_name, database_path):
    padmet_path = os.path.join(run_name, 'Networks', 'PADMETs')
    table_path = os.path.join(run_name, 'Networks', 'Tables')
    if not os.path.exists(table_path):
        os.mkdir(table_path)
    padmet_ref = PadmetRef(database_path)
    if not len(os.listdir(table_path)) == 4:
        compare_padmet(padmet_path=padmet_path, output=table_path, padmetRef=padmet_ref, verbose=True)


def main(run_name, database_path, nb_cpu=1):
    generate_input(run_name)
    networks_reconstruction(run_name, nb_cpu)
    create_padmet(run_name, database_path)
    create_sbml(run_name)
    generate_tables(run_name, database_path)


if __name__ == "__main__":
    RUN_NAME = sys.argv[1]
    DB = sys.argv[2]
    CPU = int(sys.argv[3])
    main(RUN_NAME, DB, CPU)
