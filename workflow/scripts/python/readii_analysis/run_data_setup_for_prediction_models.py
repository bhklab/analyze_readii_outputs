import argparse
import os
from pathlib import Path

# from processing import 
from readii_analysis.data.helpers import (
    makeProcessedDataFolders, 
    loadImageDatasetConfig, 
    loadFileToDataFrame,
    subsetDataframe
)

from readii_analysis.data.labelling import (
    timeOutcomeColumnSetup,
    eventOutcomeColumnSetup
)

##### ARGUMENT INPUT #####
parser = argparse.ArgumentParser(description="Run data setup for prediction models.")
parser.add_argument("--dataset_name", type=str, help="Name of the dataset to run data setup for.")
parser.add_argument("--extraction_methods", type=str, nargs="+", help="List of extraction methods to run data setup for.")
args = parser.parse_args()

# Input arguments
DATASET_NAME = args.dataset_name
EXTRACTION_METHODS = args.extraction_methods

# Set variables
RAW_DATA_PATH = Path("../../../rawdata/")
PROC_DATA_PATH = Path("../../../procdata/")
RESULTS_DATA_PATH = Path("../../../results/")
CONFIG_DIR_PATH = Path("../../config/")

# Load config file
config = loadImageDatasetConfig(DATASET_NAME, CONFIG_DIR_PATH)

# Make output directories for this pipeline
makeProcessedDataFolders(dataset_name=DATASET_NAME,
                         proc_data_path=PROC_DATA_PATH,
                         data_sources=EXTRACTION_METHODS,
                         data_types=['clinical', 'features'],
                         train_test_split=config["train_test_split"]["split"])

# Load clinical data
clinical_data = loadFileToDataFrame(os.path.join(RAW_DATA_PATH, DATASET_NAME, "clinical", f"{DATASET_NAME}.csv"))
print(f"Clinical data loaded with {len(clinical_data)} patients.\n")

# Clean clinical data
exclusion_clinical_variables = config["exclusion_variables"]
if exclusion_clinical_variables:
    print("Will exclude clinical variables:", exclusion_clinical_variables)
    # Drop rows with values in the exclusion variables
    clinical_data = subsetDataframe(clinical_data, excludeDict=exclusion_clinical_variables)
    print("Clinical data updated, now has", len(clinical_data), "patients.\n")
else:
    print("No exclusion variables found in config file.\n")


# Outcome Variable setup
clinical_data = timeOutcomeColumnSetup(clinical_data, 
                                       outcome_column_label=config["outcome_variables"]["time_label"], 
                                       standard_column_label="survival_time_in_years",
                                       convert_to_years=config["outcome_variables"]["convert_to_years"])

clinical_data = eventOutcomeColumnSetup(clinical_data,
                                        outcome_column_label=config["outcome_variables"]["event_label"],
                                        standard_column_label="survival_event_binary",
                                        event_column_value_mapping=config["outcome_variables"]["event_value_mapping"])







