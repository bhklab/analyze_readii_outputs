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

from readii_analysis.data.processing import (
    imageTypesFeatureProcessing
)



def run_data_setup_for_prediction_models(DATASET_NAME:str, EXTRACTION_METHOD:str, RAW_FEATURE_DIR_NAME:str):
    print(f"Running data setup for {DATASET_NAME} {EXTRACTION_METHOD} feature data.")

    # Set variables
    RAW_DATA_PATH = Path("../../../rawdata/")
    PROC_DATA_PATH = Path("../../../procdata/")
    CONFIG_DIR_PATH = Path("../../config/")

    # Load config file
    config = loadImageDatasetConfig(DATASET_NAME, CONFIG_DIR_PATH)

    # MAKE OUTPUT DIRECTORIES
    # Make output directories for this pipeline
    makeProcessedDataFolders(dataset_name=DATASET_NAME,
                            proc_data_path=PROC_DATA_PATH,
                            data_sources=EXTRACTION_METHOD,
                            data_types=['clinical', 'features/features_with_metadata', 'features/labelled_features_only'],
                            train_test_split=config["train_test_split"]["split"])

    #CLINICAL DATA PROCESSING
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

    # print(config["outcome_variables"]["event_value_mapping"])

    # print(clinical_data[config["outcome_variables"]["event_label"]].str.lower().unique())
    clinical_data = eventOutcomeColumnSetup(clinical_data,
                                            outcome_column_label=config["outcome_variables"]["event_label"],
                                            standard_column_label="survival_event_binary",
                                            event_column_value_mapping=config["outcome_variables"]["event_value_mapping"])

    # save out cleaned clinical data
    clinical_data.to_csv(os.path.join(PROC_DATA_PATH, DATASET_NAME, "clinical", f"cleaned_clinical_{DATASET_NAME}.csv"))

    # TODO: decide if setting patient ID as index here, or in intersection call 
    # TODO: save out clinical data at this point

    #%% IMAGE FEATURE PROCESSING
    # Construct path to the directory containing the raw image feature files
    feature_dir_path = os.path.join(RAW_DATA_PATH, DATASET_NAME, RAW_FEATURE_DIR_NAME)

    image_types_found = imageTypesFeatureProcessing(raw_data_dir=feature_dir_path,
                                feature_type=EXTRACTION_METHOD,
                                proc_data_path=PROC_DATA_PATH,
                                clinical_data=clinical_data,
                                dataset_name=DATASET_NAME,
                                outcome_labels=["survival_time_in_years", "survival_event_binary"],
                                train_test_split_settings=config['train_test_split'],
                                )
    
    # TODO: save image_types_found to config or to a file


    print(f"{DATASET_NAME} {EXTRACTION_METHOD} feature data has been set up for prediction models.")


if __name__ == "__main__":
    ##### ARGUMENT INPUT #####
    parser = argparse.ArgumentParser(description="Run data setup for prediction models.")
    parser.add_argument("--dataset_name", type=str, help="Name of the dataset to run data setup for.")
    parser.add_argument("--extraction_method", type=str, help="Extraction methods to run data setup for.")
    parser.add_argument("--raw_feature_dir", type=str, help="Path to the directory containing the raw feature files.")
    args = parser.parse_args()

    # Input arguments
    DATASET_NAME = args.dataset_name
    EXTRACTION_METHOD = args.extraction_method
    RAW_FEATURE_DIR_NAME = args.raw_feature_dir

    run_data_setup_for_prediction_models(DATASET_NAME, EXTRACTION_METHOD, RAW_FEATURE_DIR_NAME)






