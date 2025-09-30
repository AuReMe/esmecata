import csv
import os
import subprocess
import shutil

from esmecata.precomputed.create_database import create_database_from_esmecata_workflow_run
from esmecata.core.precomputed import precomputed_parse_affiliation
from esmecata.report.workflow_create_report import run_create_workflow_report

RESULTS = {
    'Cluster_1': {'proteomes': 2, 'protein_clusters': 603, 'GOs': 872, 'ECs': 314}
}

def test_workflow_online():
    buchnera_workflow_file = os.path.join('test_data', 'buchnera_workflow.tsv')

    subprocess.call(['esmecata', 'workflow_uniprot', '-i', buchnera_workflow_file, '-o', 'test_output', '-p', '1', '--minimal-nb-proteomes', '2', '--bioservices'])

    output_stat_file = os.path.join('test_output', 'stat_number_workflow.tsv')

    predicted_results = {}
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.reader(stat_file_read, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            predicted_results[line[0]] = {}
            predicted_results[line[0]]['proteomes'] = int(line[5])
            predicted_results[line[0]]['protein_clusters'] = int(line[7])
            predicted_results[line[0]]['GOs'] = int(line[8])
            predicted_results[line[0]]['ECs'] = int(line[9])

    for observation_name in predicted_results:
        for data in predicted_results[observation_name]:
            assert predicted_results[observation_name][data] == RESULTS[observation_name][data]


def test_create_database_from_workflow_online():
    buchnera_workflow_file = os.path.join('test_data', 'buchnera_workflow.tsv')
    esmecata_result_folder = 'test_output'
    esmecata_database_output_folder = 'output_database'
    output_folder = 'output_folder'
    report_output_folder = 'output_report_folder'

    create_database_from_esmecata_workflow_run(esmecata_result_folder, esmecata_database_output_folder, '0.1.0')

    esmecata_database = os.path.join(esmecata_database_output_folder, 'esmecata_database.zip')
    precomputed_parse_affiliation(buchnera_workflow_file, esmecata_database, output_folder)

    run_create_workflow_report(buchnera_workflow_file, output_folder, report_output_folder)

    html_report_file = os.path.join(report_output_folder, 'esmecata_summary.html')
    # Check creation of HTML report file.
    assert os.path.exists(html_report_file)

    shutil.rmtree(output_folder)
    shutil.rmtree(report_output_folder)
    shutil.rmtree(esmecata_database_output_folder)


def test_create_database_from_workflow_cli_online():
    buchnera_workflow_file = os.path.join('test_data', 'buchnera_workflow.tsv')
    esmecata_result_folder = 'test_output'
    esmecata_database_output_folder = 'output_database'
    output_folder = 'output_folder'

    subprocess.call(['esmecata_create_db', 'from_workflow', '-i', esmecata_result_folder, '-o', esmecata_database_output_folder, '--db-version', '0.1.0'])

    esmecata_database = os.path.join(esmecata_database_output_folder, 'esmecata_database.zip')
    precomputed_parse_affiliation(buchnera_workflow_file, esmecata_database, output_folder)

    report_output_folder = 'output_report_folder'
    run_create_workflow_report(buchnera_workflow_file, output_folder, report_output_folder)

    html_report_file = os.path.join(report_output_folder, 'esmecata_summary.html')
    # Check creation of HTML report file.
    assert os.path.exists(html_report_file)

    shutil.rmtree(output_folder)
    shutil.rmtree(report_output_folder)
    shutil.rmtree(esmecata_database_output_folder)
    shutil.rmtree(esmecata_result_folder)


if __name__ == "__main__":
    test_workflow_online()
