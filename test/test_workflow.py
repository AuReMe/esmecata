import csv
import os
import subprocess
import shutil

from esmecata.precomputed.create_database import create_database_from_esmecata_workflow_run
from esmecata.core.precomputed import precomputed_parse_affiliation
from esmecata.report.workflow_create_report import run_create_workflow_report

RESULTS = {
    'Cluster_1': {'proteomes': 2, 'protein_clusters': 603, 'GOs': 856, 'ECs': 313}
}

def test_workflow_online():
    subprocess.call(['esmecata', 'workflow_uniprot', '-i', 'buchnera_workflow.tsv', '-o', 'test_output', '-p', '1', '--minimal-nb-proteomes', '2', '--bioservices'])

    output_stat_file = os.path.join('test_output', 'stat_number_workflow.tsv')

    expected_results = {}
    with open(output_stat_file, 'r') as stat_file_read:
        csvreader = csv.reader(stat_file_read, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            expected_results[line[0]] = {}
            expected_results[line[0]]['proteomes'] = int(line[5])
            expected_results[line[0]]['protein_clusters'] = int(line[7])
            expected_results[line[0]]['GOs'] = int(line[8])
            expected_results[line[0]]['ECs'] = int(line[9])

    for observation_name in expected_results:
        for data in expected_results[observation_name]:
            assert expected_results[observation_name][data] == RESULTS[observation_name][data]


def test_create_database_from_workflow_online():
    esmecata_result_folder = 'test_output'
    esmecata_database_output_folder = 'output_database'
    create_database_from_esmecata_workflow_run(esmecata_result_folder, esmecata_database_output_folder)

    esmecata_database = os.path.join(esmecata_database_output_folder, 'esmecata_database.zip')
    output_folder = 'output_folder'
    precomputed_parse_affiliation('buchnera_workflow.tsv', esmecata_database, output_folder)

    report_output_folder = 'output_report_folder'
    run_create_workflow_report('buchnera_workflow.tsv', output_folder, report_output_folder)

    html_report_file = os.path.join(report_output_folder, 'esmecata_summary.html')
    # Check creation of HTML report file.
    assert os.path.exists(html_report_file)


def test_create_database_from_workflow_cli_online():
    esmecata_result_folder = 'test_output'
    esmecata_database_output_folder = 'output_database'
    subprocess.call(['esmecata_create_db', 'from_workflow', '-i', esmecata_result_folder, '-o', esmecata_database_output_folder])

    create_database_from_esmecata_workflow_run(esmecata_result_folder, esmecata_database_output_folder)

    esmecata_database = os.path.join(esmecata_database_output_folder, 'esmecata_database.zip')
    output_folder = 'output_folder'
    precomputed_parse_affiliation('buchnera_workflow.tsv', esmecata_database, output_folder)

    report_output_folder = 'output_report_folder'
    run_create_workflow_report('buchnera_workflow.tsv', output_folder, report_output_folder)

    html_report_file = os.path.join(report_output_folder, 'esmecata_summary.html')
    # Check creation of HTML report file.
    assert os.path.exists(html_report_file)

    shutil.rmtree('test_output')


if __name__ == "__main__":
    test_workflow_online()
