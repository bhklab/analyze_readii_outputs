import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
from typing import Optional

from readii_analysis.data.helpers import ( 
    loadFeatureFilesFromImageTypes, 
    loadImageDatasetConfig, 
    savePlotFigure, 
    makeProcessedDataFolders)

from readii_analysis.analyze.correlation_functions import ( 
    getFeatureCorrelations, 
    makeBothSelfCorrelationPlots,
    makeBothCrossCorrelationPlots)


def run_correlation_analysis(dataset_name:str, extraction_method:str, extracted_feature_dir:str,
                             correlation_method:Optional[str] = "pearson", corr_color_map:Optional[str] = "nipy_spectral",
                             self_dist_num_bins:Optional[int] = 450, self_dist_y_upper_bound:Optional[int] = None,
                             cross_dist_num_bins:Optional[int] = 450, cross_dist_y_upper_bound:Optional[int] = None):
    """ Function to run correlation analysis for a dataset.

    Parameters
    ----------
    dataset_name : str
        Name of the dataset to run correlation analysis for.
    extraction_method : str
        Name of the extraction method to use for the correlation analysis. Must be either "radiomic" or "deep_learning".
    extracted_feature_dir : str
        Path to the directory containing the extracted feature csv files.
    correlation_method : str, optional
        Correlation method to use for the correlation analysis. The default is "pearson".
    corr_color_map : str, optional
        Color map to use for the correlation plots. The default is "nipy_spectral".
    self_dist_num_bins : int, optional
        Number of bins to use for the self-correlation distribution plots. The default is 450.
    self_dist_y_upper_bound : int, optional
        Upper bound for the y-axis of the self-correlation distribution plots. The default is None.
    cross_dist_num_bins : int, optional
        Number of bins to use for the cross-correlation distribution plots. The default is 450.
    cross_dist_y_upper_bound : int, optional
        Upper bound for the y-axis of the cross-correlation distribution plots. The default is None.

    Returns
    -------
    None
    """

    print(f"Running correlation analysis for {dataset_name} {extraction_method} feature data.")

    PROC_DATA_PATH = "../../../procdata/"
    RESULTS_DATA_PATH = f"../../../results/"
    CONFIG_DIR_PATH = "../../config/"

    # Load config file
    config = loadImageDatasetConfig(dataset_name, CONFIG_DIR_PATH)

    # Results data folder creation
    makeProcessedDataFolders(dataset_name=dataset_name,
                            proc_data_path=RESULTS_DATA_PATH,
                            data_sources=extraction_method,
                            data_types=["correlation_distribution_plots", "correlation_heatmap_plots"])
    heatmap_dir_path=os.path.join(RESULTS_DATA_PATH, dataset_name, extraction_method.lower(), "correlation_heatmap_plots")
    dist_dir_path=os.path.join(RESULTS_DATA_PATH, dataset_name, extraction_method.lower(), "correlation_distribution_plots")

    
    # Load in all processed feature sets - should have one file for each image type in the extracted_image_feature_dir
    extracted_image_feature_dir = os.path.join(PROC_DATA_PATH, dataset_name, extraction_method, extracted_feature_dir)
    image_feature_sets = loadFeatureFilesFromImageTypes(extracted_feature_dir = extracted_image_feature_dir,
                                                        image_types=config["image_types"],
                                                        drop_labels=True,
                                                        labels_to_drop=["patient_ID","survival_time_in_years","survival_event_binary"])


    image_type_list = list(image_feature_sets.keys())
    # print("Feature sets available for analysis:")
    # for image_type in image_type_list:
    #     print("  ->", image_type)

    # Get list of negative control image types 
    negative_control_list = image_type_list
    negative_control_list.remove("original")

    make_original_plots = True
    for negative_control in negative_control_list:
        print(f"Calculating correlations for original vs. {negative_control}")
        
        # Going to perform Pearson correlation between feature sets from the original and each negative control image
        # y-axis is original (so vertical = original from now on)
        vertical_features = image_feature_sets["original"]
        vertical_feature_count = len(vertical_features.columns)

        # x-axis is the negative control (so horizontal = negative control from now on)
        horizontal_features = image_feature_sets[negative_control]
        horizontal_feature_count = len(horizontal_features.columns)

        # Perform the correlation - will return a dataframe that has length and width of 2x the number of features
        feature_correlation_matrix = getFeatureCorrelations(vertical_features = vertical_features,
                                                            horizontal_features = horizontal_features,
                                                            method = correlation_method,
                                                            vertical_feature_name="original",
                                                            horizontal_feature_name=negative_control)

       ################### SELF-CORRELATION PLOTS ###################
        # Only save out the original vs original plots once
        if make_original_plots:
            print("Making original self-correlation plots")
            make_original_plots = False

            original_self_corr_plot, original_self_corr_dist_plot = makeBothSelfCorrelationPlots(correlation_matrix = feature_correlation_matrix,
                                     axis = "vertical",
                                     num_axis_features = vertical_feature_count,
                                     feature_name = "original",
                                     corr_cmap = corr_color_map,
                                     dist_num_bins = self_dist_num_bins,
                                     dist_y_upper_bound = self_dist_y_upper_bound,
                                     correlation_method = correlation_method,
                                     extraction_method = extraction_method,
                                     dataset_name = dataset_name)
            plt.close('all')
            
            # Save out the correlation heatmap
            savePlotFigure(original_self_corr_plot,
                           plot_name=f"{dataset_name}_{correlation_method.lower()}_corr_original_v_original_{extraction_method.lower()}_plot.png",
                           output_dir_path=heatmap_dir_path)
            
            # Save out the distribution plot
            savePlotFigure(original_self_corr_dist_plot,
                           plot_name=f"{dataset_name}_{correlation_method.lower()}_corr_dist_original_v_original_{extraction_method.lower()}_plot.png",
                           output_dir_path=dist_dir_path)
        # End make original plots section

        print(f"Making {negative_control} self-correlation plots...")
        # Plot the correlation heatmap for the negative control vs negative control
        negative_control_self_corr_plot, negative_control_self_corr_dist_plot = makeBothSelfCorrelationPlots(correlation_matrix = feature_correlation_matrix,
                                                                                                    axis = "horizontal",
                                                                                                    num_axis_features = horizontal_feature_count,
                                                                                                    feature_name = negative_control,
                                                                                                    corr_cmap = corr_color_map,
                                                                                                    dist_num_bins = self_dist_num_bins,
                                                                                                    dist_y_upper_bound = self_dist_y_upper_bound,
                                                                                                    correlation_method = correlation_method,
                                                                                                    extraction_method = extraction_method,
                                                                                                    dataset_name = dataset_name)
        plt.close('all')

        # Save out the correlation heatmap
        savePlotFigure(negative_control_self_corr_plot,
                       plot_name=f"{dataset_name}_{correlation_method.lower()}_corr_{negative_control}_v_{negative_control}_{extraction_method.lower()}_plot.png",
                       output_dir_path=heatmap_dir_path)
        
        # Save out the distribution plot
        savePlotFigure(negative_control_self_corr_dist_plot,
                       plot_name=f"{dataset_name}_{correlation_method.lower()}_corr_dist_{negative_control}_v_{negative_control}_{extraction_method.lower()}_plot.png",
                       output_dir_path=dist_dir_path)
        

        ################### CROSS-CORRELATION PLOTS ###################
        print(f"Making original vs. {negative_control} cross-correlation plots...")

        orig_vs_neg_control_cross_corr_plot, orig_vs_neg_control_cross_corr_dist_plot = makeBothCrossCorrelationPlots(correlation_matrix = feature_correlation_matrix,
                                                                              num_vertical_features = vertical_feature_count,
                                                                              vertical_feature_name = "original",
                                                                              horizontal_feature_name =  negative_control,
                                                                              corr_cmap = corr_color_map,
                                                                              dist_num_bins = cross_dist_num_bins,
                                                                              dist_y_upper_bound = cross_dist_y_upper_bound,
                                                                              correlation_method = correlation_method,
                                                                              extraction_method = extraction_method,
                                                                              dataset_name = dataset_name)
        plt.close('all')
        
        
        savePlotFigure(orig_vs_neg_control_cross_corr_plot,
                       plot_name=f"{dataset_name}_{correlation_method.lower()}_corr_original_v_{negative_control}_{extraction_method.lower()}_plot.png",
                       output_dir_path=heatmap_dir_path)
        
        savePlotFigure(orig_vs_neg_control_cross_corr_dist_plot,
                       plot_name=f"{dataset_name}_{correlation_method.lower()}_corr_dist_original_v_{negative_control}_{extraction_method.lower()}_plot.png",
                       output_dir_path=dist_dir_path)

    print("Done!")





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