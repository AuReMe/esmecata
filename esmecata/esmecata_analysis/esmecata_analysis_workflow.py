import os

from esmecata.esmecata_analysis.analysis import create_dataset_annotation_file
from esmecata.esmecata_analysis.ec_sunburst_per_model import create_sunburst_ec
from esmecata.esmecata_analysis.viz_esmecata_datapane import create_datapane

def run_analysis_pipeline(input_file, input_folder, output_folder):
    annotation_referene_folder_path = os.path.join(input_folder, '2_annotation', 'annotation_reference')
    dataset_annotation_ec_file_path = os.path.join(output_folder, 'dataset_annotation_ec.tsv')
    create_dataset_annotation_file(annotation_referene_folder_path, dataset_annotation_ec_file_path, content='EC')

    dataset_annotation_go_file_path = os.path.join(output_folder, 'dataset_annotation_go.tsv')
    create_dataset_annotation_file(annotation_referene_folder_path, dataset_annotation_go_file_path, content='GO')
    #create_sunburst_ec(input_folder, output_folder)
    create_datapane(input_file, input_folder, output_folder)

