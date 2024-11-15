import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse

from readii_analysis.data.helpers import ( 
    loadFeatureFilesFromImageTypes, 
    loadImageDatasetConfig, 
    savePlotFigure, 
    makeProcessedDataFolders)

from readii_analysis.analyze.correlation_functions import ( 
    getFeatureCorrelations, 
    plotCorrelationHeatmap, 
    plotCorrelationDistribution)


def run_correlation_analysis(dataset_name:str, extraction_method:str, extracted_feature_dir:str):
    print(f"Running correlation analysis for {DATASET_NAME} {EXTRACTION_METHOD} feature data.")

    PROC_DATA_PATH = "../../../procdata/"
    RESULTS_DATA_PATH = f"../../../results/"
    CONFIG_DIR_PATH = "../../config/"

    # Load config file
    config = loadImageDatasetConfig(dataset_name, CONFIG_DIR_PATH)

    # Results data folder creation
    makeProcessedDataFolders(dataset_name=dataset_name,
                            proc_data_path=RESULTS_DATA_PATH,
                            data_sources=extraction_method,
                            data_types=["correlation_distribution_plots", "correlation_heatmap_plots"],
                            train_test_split=config["train_test_split"]["split"])
    

    extracted_image_feature_dir = os.path.join(PROC_DATA_PATH, dataset_name, extraction_method, extracted_feature_dir)
    image_feature_sets = loadFeatureFilesFromImageTypes(extracted_feature_dir = extracted_image_feature_dir,
                                                        image_types=config["image_types"],
                                                        drop_labels=True)

    pass



if __name__ == "__main__":
    ##### ARGUMENT INPUT #####
    parser = argparse.ArgumentParser(description="Run data setup for prediction models.")
    parser.add_argument("--dataset_name", type=str, help="Name of the dataset to run data setup for.")
    parser.add_argument("--extraction_method", type=str, help="Extraction methods to run data setup for.")
    parser.add_argument("--extracted_feature_dir", type=str, help="Path to the directory containing the raw feature files.")
    args = parser.parse_args()

    # Input arguments
    DATASET_NAME = args.dataset_name
    EXTRACTION_METHOD = args.extraction_method
    EXTRACTED_FEATURE_DIR_NAME = args.extracted_feature_dir

    # print(DATASET_NAME, EXTRACTION_METHOD, EXTRACTED_FEATURE_DIR_NAME)
    run_correlation_analysis(DATASET_NAME, EXTRACTION_METHOD, EXTRACTED_FEATURE_DIR_NAME)