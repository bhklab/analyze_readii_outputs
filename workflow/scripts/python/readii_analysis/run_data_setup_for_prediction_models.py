import argparse
import os
from pathlib import Path

# from processing import 
from readii_analysis.data.helpers import makeProcessedDataFolders, loadImageDatasetConfig, loadFileToDataFrame

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
print(f"Clinical data loaded with {len(clinical_data)} patients.")



