import csv
import os
import shutil
import subprocess

from esmecata.report.workflow_create_report import run_create_workflow_report
from esmecata.core.precomputed import precomputed_parse_affiliation


def test_create_report_from_precomputed_offline():
    output_folder = 'output_folder'
    buchnera_workflow_file = os.path.join('test_data', 'buchnera_workflow.tsv')
    buchnera_database_file = os.path.join('test_data', 'buchnera_database.zip')
    precomputed_parse_affiliation(buchnera_workflow_file, buchnera_database_file, output_folder)

    report_output_folder = 'output_report_folder'
    run_create_workflow_report(buchnera_workflow_file, output_folder, report_output_folder)

    html_report_file = os.path.join(report_output_folder, 'esmecata_summary.html')
    # Check creation of HTML report file.
    assert os.path.exists(html_report_file)


def test_create_report_from_precomputed_cli_offline():
    output_folder = 'output_folder'
    buchnera_workflow_file = os.path.join('test_data', 'buchnera_workflow.tsv')
    buchnera_database_file = os.path.join('test_data', 'buchnera_database.zip')

    subprocess.call(['esmecata', 'precomputed', '-i', buchnera_workflow_file, '-d', buchnera_database_file, '-o', output_folder])

    report_output_folder = 'output_report_folder'
    subprocess.call(['esmecata_report', 'create_report', '-i', buchnera_workflow_file, '-f', output_folder, '-o', report_output_folder])

    html_report_file = os.path.join(report_output_folder, 'esmecata_summary.html')
    # Check creation of HTML report file.
    assert os.path.exists(html_report_file)

if __name__ == "__main__":
    test_create_report_from_precomputed_offline()
    test_create_report_from_precomputed_cli_offline()
