import os
from pandas import DataFrame
import pandas as pd

from .helpers import getPatientIdentifierLabel, getOnlyPyradiomicsFeatures, loadFileToDataFrame, splitDataSetup

from typing import Optional, Union
from pathlib import Path

def setPatientIdAsIndex(dataframe_to_index:DataFrame,
                        patient_id_col:str = None):
    """ Function to set the patient ID column as the index of a dataframe.

    Parameters
    ----------
    dataframe : DataFrame
        Dataframe to set the patient ID column as the index.
    patient_id_col : str, optional
        Name of the patient ID column to use as the index. If not provided, will find the patient ID column with getPatientIdentifierLabel.
    """
    if not patient_id_col:
        patient_id_col = getPatientIdentifierLabel(dataframe_to_index)
        
    pat_indexed_dataframe = dataframe_to_index.set_index(patient_id_col)
    return pat_indexed_dataframe
    


def getPatientIntersectionDataframes(dataframe_A:DataFrame, 
                              dataframe_B:DataFrame,
                              need_pat_index_A:bool = True,
                              need_pat_index_B:bool = True): 
    """ Function to get the subset of two dataframes based on the intersection of their indices. Intersection will be based on the index of dataframe A.

    Parameters
    ----------
    dataframe_A : DataFrame
        Dataframe A to get the intersection of based on the index.
    dataframe_B : DataFrame
        Dataframe B to get the intersection of based on the index.
    need_pat_index_A : bool, optional
        Whether to run setPatientIdAsIndex on dataframe A. If False, assumes the patient ID column is already set as the index. The default is True.
    need_pat_index_B : bool, optional
        Whether to run setPatientIdAsIndex on dataframe B. If False, assumes the patient ID column is already set as the index. The default is True.

    Returns
    -------
    intersection_index_dataframeA : DataFrame
        Dataframe containing the rows of dataframe A that are in the intersection of the indices of dataframe A and dataframe B.
    intersection_index_dataframeB : DataFrame
        Dataframe containing the rows of dataframe B that are in the intersection of the indices of dataframe A and dataframe B.
    """
    
    # Set the patient ID column as the index for dataframe A if needed
    if need_pat_index_A:
        dataframe_A = setPatientIdAsIndex(dataframe_A)

    # Set the patient ID column as the index for dataframe B if needed
    if need_pat_index_B:
        dataframe_B = setPatientIdAsIndex(dataframe_B)

    # Get patients in common between dataframe A and dataframe B
    intersection_index = dataframe_A.index.intersection(dataframe_B.index)

    # Select common patient rows from each dataframe
    intersection_index_dataframeA = dataframe_A.loc[intersection_index]
    intersection_index_dataframeB = dataframe_B.loc[intersection_index]

    return intersection_index_dataframeA, intersection_index_dataframeB



def addOutcomeLabels(feature_data_to_label:DataFrame,
                      clinical_data:DataFrame,
                      outcome_labels:list = ["survival_time_in_years", "survival_event_binary"]):
    """ Function to add survival labels to a feature dataframe based on a clinical dataframe.

    Parameters
    ----------
    feature_data_to_label : DataFrame
        Dataframe containing the feature data to add survival labels to.
    clinical_data : DataFrame
        Dataframe containing the clinical data to use for survival labels.
    outcome_labels : list, optional
        List of outcome labels to extract from the clinical dataframe. The default is ["survival_time_in_years", "survival_event_binary"].
    """
    # Get the survival time and event columns as a dataframe
    outcome_label_columns = clinical_data[outcome_labels]

    # Join the outcome label dataframe to the feature data dataframe
    outcome_labelled_feature_data = outcome_label_columns.join(feature_data_to_label)
    return outcome_labelled_feature_data



def getImageTypesFromDirectory(raw_data_dir:str,
                  feature_file_prefix:str = "",
                  feature_file_suffix:str = ".csv"):
    """ Function to get a list of image types from a directory containing image feature files.

    Parameters
    ----------
    raw_data_dir : str
        Path to the directory containing the image feature files.
    feature_file_prefix : str, optional
        Prefix to remove from the feature file name. The default is "".
    feature_file_suffix : str, optional
        Suffix to remove from the feature file name. The default is ".csv".
    
    Returns
    -------
    list
        List of image types from the image feature files.
    """

    return sorted([file.removeprefix(feature_file_prefix).removesuffix(feature_file_suffix) for file in os.listdir(raw_data_dir)])



def trainTestSplitSetup(clinical_data:DataFrame,
                        feature_data:DataFrame,
                        split_variable:str,
                        proc_data_path:str,
                        impute_value:Optional[str] = "",
                        clinical_out_file_suffix:Optional[str] = "",
                        feature_out_file_suffix:Optional[str] = "",
                        ):
    """ Function to split clinical and feature sets into training and test sets based on a label in the clinical data.
        The clinical and feature data will be saved out to the processed data folder.

    Parameters
    ----------
    clinical_data : DataFrame
        Dataframe containing the clinical data to use for train test splits. Must include a column name matching split_variable. Index must match the index of the feature data.
    feature_data : DataFrame
        Dataframe containing the feature data to use for train test splits. Index must match the index of the clinical data.
    split_variable : str
        Name of the column to split on in the clinical data.
    proc_data_path : str
        Path to the processed data folder that all folders will be made under.
    impute_value : str, optional
        Value to impute for missing values in the split variable. The default is "".
    clinical_out_file_suffix : str, optional
        Suffix to add to the clinical output file name. The default will be clinical.csv.
    feature_out_file_suffix : str, optional
        Suffix to add to the feature output file name. The default will be features.csv.

    Returns
    -------
    split_clinical : dict
        Dictionary containing the clinical data for the training and test sets.
    split_features : dict
        Dictionary containing the feature data for the training and test sets.
    """
    # Handle feature and clinical output file names
    if not clinical_out_file_suffix:
        clinical_out_file_suffix = "clinical.csv"
    # Check that clinical_out_file_suffix ends in ".csv" and add if it doesn't
    if not clinical_out_file_suffix.endswith(".csv"):
        clinical_out_file_suffix += ".csv"

    # Check that feature_out_file_suffix is not blank
    if not feature_out_file_suffix:
        feature_out_file_suffix = "features.csv"
    # Check that feature_out_file_suffix ends in ".csv" and add if it doesn't
    if not feature_out_file_suffix.endswith(".csv"):
        feature_out_file_suffix += ".csv"


    # Set up train test output path
    train_test_procdata_path = os.path.join(proc_data_path, "train_test_split")

    # Split the data into training and test sets, imputing rows with missing values
    split_clinical, split_features = splitDataSetup(clinical_data, feature_data, 
                                                    splitVariables = split_variable, 
                                                    imputeValue = impute_value)
    
    # Save out training and test clinical data
    split_clinical['training'].to_csv(os.path.join(train_test_procdata_path, f"clinical/train_{clinical_out_file_suffix}"))
    split_clinical['test'].to_csv(os.path.join(train_test_procdata_path, f"clinical/test_{clinical_out_file_suffix}"))

    # Save out training and test radiomic feature data
    split_features['training'].to_csv(os.path.join(train_test_procdata_path, f"train_features/train_{feature_out_file_suffix}"))
    split_features['test'].to_csv(os.path.join(train_test_procdata_path, f"test_features/test_{feature_out_file_suffix}"))

    print(f"Training and test splits saved to {train_test_procdata_path}")

    return split_clinical, split_features



def featureProcessingForPrediction(raw_image_data:DataFrame,
                                 clinical_data:DataFrame,
                                 outcome_labels:list = ["survival_time_in_years", "survival_event_binary"],
                                            ):
    """ Function to process a set of image features for prediction. Will return the intersected clinical and image data, and the image features with outcome labels added.

    Parameters
    ----------
    raw_image_data : DataFrame
        Dataframe containing the raw image features to process.
    clinical_data : DataFrame
        Dataframe containing the clinical data to use for survival labels.
    outcome_labels : list, optional
        List of outcome labels to extract from the clinical dataframe. The default is ["survival_time_in_years", "survival_event_binary"].

    Returns
    -------
    common_clinical_data : DataFrame
        Dataframe containing the clinical data for patients with both clinical and image features.
    common_image_data : DataFrame
        Dataframe containing the image data for patients with both clinical and image features.
    outcome_labelled_image_features : DataFrame
        Dataframe containing the image features from the intersected clinical and image data with outcome labels added.
    """
    
    # Get patient ID column name
    patient_identifier = getPatientIdentifierLabel(raw_image_data)

    # Check for duplicated rows and remove the second copy
    image_data = raw_image_data.drop_duplicates(subset=[patient_identifier, "series_UID","image_modality","seg_modality","seg_ref_image", "roi"])

    # Set patient ID as index for image data
    image_data.set_index(patient_identifier, inplace=True)

    # Filter the clinical and image features to only include patients with imaging and clinical data based on image features index
    # e.g. patients with only clinical data will not be included
    # Index of returned dataframes will be the patient IDs 
    common_image_data, common_clinical_data = getPatientIntersectionDataframes(image_data, clinical_data, need_pat_index_A=False, need_pat_index_B=True)

    print(f"Common patient count: {len(common_clinical_data)}")
    print(f"Number of segmentations: {len(common_image_data)}")

    # Get just the radiomic feature columns from the dataframe, remove any metadata/diagnostics columns
    image_features = getOnlyPyradiomicsFeatures(common_image_data)
    print(f"Number of radiomic features: {(image_features.shape[1])}")

    # Add survival labels to the image features
    outcome_labelled_image_features = addOutcomeLabels(image_features, common_clinical_data, outcome_labels)
    
    return common_clinical_data, common_image_data, outcome_labelled_image_features 



def imageTypesFeatureProcessing(raw_data_dir:str,
                       feature_type:str, # radiomic or deep_learning
                       proc_data_path:str,
                       clinical_data:DataFrame,
                       dataset_name:str,
                       outcome_labels:Optional[list] = ["survival_time_in_years", "survival_event_binary"],
                       train_test_split_settings:Optional[dict] = {"split": False, "split_variable": {}, "impute": ""},
                       ):
    """ Function to process a set of image features for prediction.

    Parameters
    ----------
    raw_data_dir : str
        Path to the directory containing the raw image feature files.
    feature_type : str
        Type of image feature to process. Something like "radiomic" or "deep_learning".
        The string will be used to construct the file names of the processed data.
    proc_data_path : str
        Path to the processed data folder that all folders will be made under.
    clinical_data : DataFrame
        Dataframe containing the clinical data to use for survival labels.
    dataset_name : str
        Name of the dataset to process.
    outcome_labels : list, optional
        List of outcome labels to extract from the clinical dataframe. The default is ["survival_time_in_years", "survival_event_binary"].
    train_test_split_settings : dict, optional
        Dictionary containing settings for training and test splits. The default is {"split": False, "split_variable": {}, "impute": ""}.

    Returns
    -------
    None

    """
    
    image_feature_file_list = sorted(os.listdir(raw_data_dir))

    feature_procdata_path = os.path.join(proc_data_path, dataset_name, feature_type)

    for feature_file in image_feature_file_list:
        # Get the image type from the feature file name
        image_type = feature_file.removeprefix(f"{feature_type}features_").removesuffix(f"_{dataset_name}.csv")
        print(f"Processing {feature_type} features for {image_type}")
        
        # Load the feature data
        feature_data = loadFileToDataFrame(os.path.join(raw_data_dir, feature_file))
    
        common_clinical_data, common_image_data, outcome_labelled_image_features = featureProcessingForPrediction(feature_data, clinical_data, outcome_labels)
        print(f"{image_type} {feature_type} feature data has been intersected with clinical data and labelled with outcome labels.")

        # Save processed data
        common_clinical_data.to_csv(os.path.join(feature_procdata_path, f"clinical/merged_clinical_{dataset_name}.csv"))
        common_image_data.to_csv(os.path.join(feature_procdata_path, f"features/merged_{feature_type}features_{image_type}_{dataset_name}.csv"))
        outcome_labelled_image_features.to_csv(os.path.join(feature_procdata_path, f"features/labelled_{feature_type}features_only_{image_type}_{dataset_name}.csv"))
        print(f"{image_type} {feature_type} feature data has been saved to {feature_procdata_path}/features")
    
        if train_test_split_settings:
            print("Splitting data into training and test sets.")

            clinical_out_file_suffix = f"merged_clinical_{dataset_name}.csv"
            feature_out_file_suffix = f"labelled_{feature_type}features_only_{image_type}_{dataset_name}.csv"
            
            # Split the clinical and labelled image features into training and test sets
            _, _ = trainTestSplitSetup(clinical_data = common_clinical_data,
                                       feature_data = outcome_labelled_image_features,
                                       split_variable = train_test_split_settings["split_variable"],
                                       proc_data_path = feature_procdata_path,
                                       impute_value = train_test_split_settings["impute"],
                                       clinical_out_file_suffix=clinical_out_file_suffix,
                                       feature_out_file_suffix=feature_out_file_suffix
                                       )
        # end train test split

        print("------------------------------------------------------------")
        print()
    # end image type loop
    return


def run_data_setup_for_prediction_models(dataset_name:str,
                                         config_file:Union[str, Path],
                                         ):
    
    pass