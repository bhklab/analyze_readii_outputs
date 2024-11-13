#%%
import sys; sys.path.insert(0, '..')

# from processing import 
from readii_analysis.data.helpers import makeProcessedDataFolders, loadImageDatasetConfig

from pathlib import Path


# %%
# Input arguments
DATASET_NAME = "RADCURE"
EXTRACTION_METHODS = ["radiomics", "deep_learning"]

# Set variables
RAW_DATA_PATH = Path("../../../rawdata/")
PROC_DATA_PATH = Path("../../../procdata/")
RESULTS_DATA_PATH = Path("../../../results/")
CONFIG_DIR_PATH = Path("../../config/")


config = loadImageDatasetConfig(DATASET_NAME, CONFIG_DIR_PATH)

# Set up processed data folders
makeProcessedDataFolders(dataset_name = DATASET_NAME,
                         proc_data_path=PROC_DATA_PATH,
                         data_sources=EXTRACTION_METHODS,
                         data_types=['clinical', 'features'],
                         train_test_split=config["train_test_split"]["split"])


