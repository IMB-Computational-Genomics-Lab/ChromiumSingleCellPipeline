#!/usr/bin/env python2.7
import os
import argparse
import xmltodict
import yaml
import string
import pandas as pd

def createAggrSheet(samplesheet_df, output_dir):
    aggr_df = samplesheet_df[["Sample_ID"]]
    molecule_h5_paths = output_dir + aggr_df["Sample_ID"] + "/outs/molecule_info.h5"
    molecule_h5_paths.name = "molecule_h5"
    aggr_df["molecule_h5"] = molecule_h5_paths
    aggr_df.columns = ["library_id", "molecule_h5"]
    return aggr_df

def createSimpleSheet(samplesheet_df):
    simple_df = samplesheet_df[["Lane", "Sample", "Index"]]
    return simple_df

def loadSampleSheet(samplesheet_dir):
    samplesheet_df = pd.read_csv(samplesheet_dir, sep = ",")
    return samplesheet_df

def getFlowcellID(xml_file):
    with open(xml_file) as xml_stream:
        xml_dict = xmltodict.parse(xml_stream.read())
        return xml_dict["RunInfo"]["Run"]["Flowcell"]

def run_aggr(project_config, pipeline_config, output_path):
    # Load correct pipeline from project_config
    pipeline = project_config["pipeline"]["name"]
    version = "version " + project_config["pipeline"]["version"]
    pipeline_path = pipeline_config["pipelines"][pipeline][version]

    # Get information from project yaml
    project_name = project_config["project_name"]
    samplesheet_dir = project_config["input_dir"]["samplesheet_dir"]
    output_dir = project_config["output_dir"]

    # Create aggregation samplesheet
    samplesheet_df = loadSampleSheet(samplesheet_dir)
    aggr_df = createAggrSheet(samplesheet_df, output_dir)

    # Write out aggregation csv
    aggr_path = os.path.join(output_path, project_name + "_Aggr.csv")
    aggr_df.to_csv(aggr_path, index = False)

    # Collect variables into dict
    variable_dict = {"pipeline_path": pipeline_path,
    "version": version,
    "aggr_filepath": aggr_path,
    "project_name": project_name,
    "output_dir": output_dir}
    # Generate Aggr PBS file
    output_filename = os.path.join(output_path, project_name + "_Aggr.pbs")
    pbs_template = open("templates/aggr_template.pbs", "r").read()
    template_output = pbs_template.format(**variable_dict)
    open(output_filename, "w").write(template_output)

def run_count(project_config, pipeline_config, output_path):
    # Load correct pipeline from project_config
    pipeline = project_config["pipeline"]["name"]
    version = "version " + project_config["pipeline"]["version"]

    # Retrieve arguments from YAML files
    pipeline_path = pipeline_config["pipelines"][pipeline][version]

    # Retrieve variables from project yaml
    bcl_paths = project_config["input_dir"]["bcl_dir"]
    fastq_basedir = project_config["input_dir"]["fastq_dir"]
    ref_dir = project_config["input_dir"]["reference_dir"]
    output_dir = project_config["output_dir"]
    samplesheet_dir = project_config["input_dir"]["samplesheet_dir"]

    # Generate paths from each flowcell
    fastq_paths = []
    for flowcell in bcl_paths:
        flowcell_id = getFlowcellID(os.path.join(flowcell, "RunInfo.xml"))
        fastq_dir = os.path.join(fastq_basedir, flowcell_id, "outs", "fastq_path")
        fastq_paths.append(fastq_dir)

    # This makes the assumption that specified flowcells are for ALL samples
    if len(fastq_paths) > 1:
        fastq_paths = ("").join(fastq_paths)
    else:
        fastq_paths = fastq_paths[0]

    # Load samplesheet
    samplesheet_df = loadSampleSheet(samplesheet_dir)
    pbs_template = open("templates/count_template.pbs", "r").read()

    # Generate a PBS script for each sample - greater control over samples
    for index, row in samplesheet_df.iterrows():
        # Variables from the spreadsheet
        sample = row["Sample"]
        sample_id = row["Sample_ID"]
        ncells = row["Cells"]

        # Collect it into a nice dictionary
        variable_dict = {"fastq_dir": fastq_paths,
        "samplesheet_dir": samplesheet_dir,
        "output_dir": output_dir,
        "sample": sample,
        "sample_id": sample_id,
        "ncells": ncells,
        "pipeline_path": pipeline_path,
        "reference_dir": ref_dir}

        # Add variables to the supplied template
        sample_pbs = pbs_template.format(**variable_dict)

        # Write it out to disk
        output_filename = os.path.join(output_path, sample_id + ".pbs")
        open(output_filename, "w").write(sample_pbs)

def run_mkfastq(project_config, pipeline_config, output_path):
    # Load correct pipeline from project_config
    project_name = project_config["project_name"]
    pipeline = project_config["pipeline"]["name"]
    version = "version " + project_config["pipeline"]["version"]

    # Get pipeline path from pipeline_config
    pipeline_path = pipeline_config["pipelines"][pipeline][version]

    # Get variables stored in project_config
    input_dict = project_config["input_dir"]
    input_dirs = input_dict["bcl_dir"]
    samplesheet_dir = input_dict["samplesheet_dir"]
    fastq_dir = input_dict["fastq_dir"]

    # Open the samplesheet and simplify the sheet so Cell Ranger won't reject it
    samplesheet_df = loadSampleSheet(samplesheet_dir)
    simple_df = createSimpleSheet(samplesheet_df)
    simplesheet_path = os.path.join(output_path, project_name + "_SimpleSheet.csv")
    simple_df.to_csv(simplesheet_path, index = False)

    # Generate a PBS script for each specified BCL directory
    for input_dir in input_dirs:
        pbs_template = open("templates/mkfastq_template.pbs", "r").read()
        variable_dict = {"input_dir": input_dir,
        "samplesheet_dir": simplesheet_path,
        "fastq_dir": fastq_dir,
        "pipeline_path": pipeline_path}
        template_output = pbs_template.format(**variable_dict)
        flowcell_id = getFlowcellID(os.path.join(input_dir, "RunInfo.xml"))
        output_filename = os.path.join(output_path, flowcell_id + ".pbs")
        open(output_filename, "w").write(template_output)

# Wrapper for opening and parsing YAML files
def readYAML(yaml_filename):
    with open(yaml_filename, "r") as yaml_stream:
        yaml_dict = yaml.safe_load(yaml_stream)
    return yaml_dict

# Argument parser for input files
def parseArgs():
    parser = argparse.ArgumentParser(prog = "scRNAseqPipe.py", description = "Generates PBS scripts for the running of scRNAseq pipelines on Delta.")
    parser.add_argument("--project", type = str, help = "YAML file containing project information.", required = True)
    parser.add_argument("--config", type = str, help = "YAML file containing pipeline information.", required = True)
    args = parser.parse_args()
    return (args.project, args.config)

if __name__ == "__main__":
    # Parse filenames
    project_yaml, config_yaml = parseArgs()

    # Retrieve data from YAML files
    project_config = readYAML(project_yaml)
    pipeline_config = readYAML(config_yaml)

    # Get absolute path of project_yaml, then use as foundation for outputs
    output_path = os.path.dirname(os.path.abspath(project_yaml))

    # Organise functions
    function_dict = {"cellranger": {
    "mkfastq": run_mkfastq(project_config, pipeline_config, output_path),
    "count": run_count(project_config, pipeline_config, output_path),
    "aggr": run_aggr(project_config, pipeline_config, output_path)
    }}

    # Get pipeline information from project_yaml
    pipeline = project_config["pipeline"]["name"]
    stages = project_config["pipeline"]["stages"]

    # Load pipeline as specified by project
    pipeline_functions = function_dict[pipeline]

    # Run function for each stage
    for stage, status in stages.items():
        if status:
            pipeline_functions[stage]
