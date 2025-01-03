import itertools
import numpy as np 
import os
import pandas as pd 
import re
import yaml

from pandas import DataFrame
from seaborn import heatmap

from typing import Optional, Dict, Union

def loadImageDatasetConfig(dataset_name:str,
                           config_dir_path:Optional[str]="../../config") -> dict:
    # Make full path to config file
    config_file_path = os.path.join(config_dir_path, f"{dataset_name}.yaml")

    # Check if config file exists
    if os.path.exists(config_file_path):
        # Load the config file
        config = yaml.safe_load(open(config_file_path, "r"))
        return config
    else:
        print(f"Config file {config_file_path} does not exist.")
        return None


def loadFileToDataFrame(file_path:str) -> pd.DataFrame:
    """Load clinical data from a csv or xlsx file into a pandas dataframe.

    Args:
        clinical_data_path (str): Path to the clinical data file.

    Returns:
        pd.DataFrame: A pandas dataframe containing the clinical data.
    """
     # Get the file extension
    _, file_extension = os.path.splitext(file_path)
    
    try:
        # Check if the file is an Excel file
        if file_extension == '.xlsx':
            df = pd.read_excel(file_path)
        # Check if the file is a CSV file
        elif file_extension == '.csv':
            df = pd.read_csv(file_path)
        else:
            raise ValueError("Unsupported file format. Please provide a .csv or .xlsx file.")
        
        return df
    
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    


def loadFeatureFilesFromImageTypes(extracted_feature_dir:str,
                                    image_types:Optional[list]=['original'], 
                                    drop_labels:Optional[bool]=True, 
                                    labels_to_drop:Optional[list]=["patient_ID","survival_time_in_years","survival_event_binary"])->Dict[str,pd.DataFrame]:
    """Function to load in all the extracted imaging feature sets from a directory and return them as a dictionary of dataframes.

    Parameters
    ----------
    extracted_feature_dir : str
        Path to the directory containing the extracted feature csv files
    image_types : list, optional
        List of image types to load in. The default is ['original'].
    drop_labels : bool, optional
        Whether to drop the labels from the dataframes. Use when loading labelled data from data_setup_for_modeling.ipynb. The default is True.
    labels_to_drop : list, optional
        List of labels to drop from the dataframes. The default is ["patient_ID","survival_time_in_years","survival_event_binary"] based on code
        in data_setup_for_modeling.ipynb.

    Returns
    -------
    feature_sets : dict
        Dictionary of dataframes containing the extracted radiomics features.
    """
    # Initialize dictionary to store the feature sets
    feature_sets = {}

    feature_file_list = os.listdir(extracted_feature_dir)

    # Loop through all the files in the directory
    for image_type in image_types:
        try:
            # Extract the image type feature csv file from the feature directory
            # This should return a list of length 1, so we can just take the first element
            image_type_feature_file = [file for file in feature_file_list if (image_type in file) and (file.endswith(".csv"))][0]
            # Remove the image type file from the list of feature files
            feature_file_list.remove(image_type_feature_file)
        except Exception as e:
            print(f"{e}\n No {image_type} feature csv files found in {extracted_feature_dir}")
            # Skip to the next image type
            continue


        # Get the full path to the feature file
        feature_file_path = os.path.join(extracted_feature_dir, image_type_feature_file)
            
        # Load the feature data into a pandas dataframe
        raw_feature_data = loadFileToDataFrame(feature_file_path)

        try:
            # Drop the labels from the dataframe if specified
            if drop_labels:
                # Data is now only extracted features
                raw_feature_data.drop(labels_to_drop, axis=1, inplace=True)
        except Exception as e:
            print(f"{feature_file_path} does not have the labels {labels_to_drop} to drop.")
            # Skip to the next image type
            continue

        # Save the dataframe to the feature_sets dictionary
        feature_sets[image_type] = raw_feature_data
    
    return feature_sets


def getPatientIdentifierLabel(dataframeToSearch:pd.DataFrame) -> str:
    """Function to find a column in a dataframe that contains some form of patient ID or case ID (case-insensitive). 
       If multiple found, will return the first match.

    Parameters
    ----------
    dataframeToSearch : DataFrame
        Dataframe to look for a patient ID column in.

    Returns
    -------
    str
        Label for patient identifier column from the dataframe.
    """

    # regex to get patient identifier column name in the dataframes
    # catches case-insensitive variations of patient_id, patid, pat id, case_id, case id, caseid, id
    regexSearchTerm = re.compile(pattern= r'(pat)?(ient)?(case)?(\s|.)?(id|#)', flags=re.IGNORECASE)

    # Get any columns from the dataframe based on the regex
    patIdentifier = dataframeToSearch.filter(regex=regexSearchTerm).columns.to_list()

    if len(patIdentifier) > 1:
        print(f"Multiple patient identifier labels found. Using {patIdentifier[0]}.")
    
    elif len(patIdentifier) == 0:
        raise ValueError("Dataframe doesn't have a recognizeable patient ID column. Must contain patient or case ID.")

    return patIdentifier[0]


def subsetDataframe(dataframe:pd.DataFrame, 
                    includeDict:Optional[dict] = None,
                    excludeDict:Optional[dict] = None) -> pd.DataFrame:
    """
    Get rows of pandas DataFrame based on row values in the columns labelled as keys of the includeDict and not in the keys of the excludeDict.
    Include variables will be processed first, then exclude variables, in the order they are provided in the corresponding dictionaries.

    Parameters
    ----------
    dataframeToSubset : pd.DataFrame
        Dataframe to subset.
    includeDict : dict
        Dictionary of column names and values to include in the subset. ex. {"column_name": ["value1", "value2"]}
    excludeDict : dict
        Dictionary of column names and values to exclude from the subset. ex. {"column_name": ["value1", "value2"]}

    Returns
    -------
    pd.DataFrame
        Subset of the input dataframe.

    """
    try:
        if (includeDict is None) and (excludeDict is None):
            raise ValueError("Must provide one of includeDict or excludeDict.")
    
        if includeDict is not None:
            for key in includeDict.keys():
                if key in ["Index", "index"]:
                    dataframe = dataframe[dataframe.index.isin(includeDict[key])]
                else:
                    dataframe = dataframe[dataframe[key].isin(includeDict[key])]

        if excludeDict is not None:
            for key in excludeDict.keys():
                if key in ["Index", "index"]:
                    dataframe = dataframe[~dataframe.index.isin(excludeDict[key])]
                else:
                    dataframe = dataframe[~dataframe[key].isin(excludeDict[key])]
        
        return dataframe

    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    

def getOnlyPyradiomicsFeatures(dfPyradiomicsFeatures:DataFrame):
    """ Function to get out just the features from a Pyradiomics output that includes metadata/diagnostics columns before the features.
        Assumes the features start after the last metadata/diagnostics column.
    Parameters
    ----------
    dfPyradiomicsFeatures : DataFrame
        Dataframe of Pyradiomics features with diagnostics and other columns before the features
    
    Returns
    -------
    featsOnlyRadiomics : DataFrame
        Dataframe with just the radiomic features
    
    """
    # Find all the columns that begin with diagnostics
    diagnosticRadiomics = dfPyradiomicsFeatures.filter(regex=r"diagnostics_*")

    if not diagnosticRadiomics.empty:
        # Get the last diagnostics column name - the features begin in the next column
        lastDiagnosticColumn = diagnosticRadiomics.columns[-1]
        # Drop all the columns before the features start
        featsOnlyRadiomics = dropUpToFeature(dfPyradiomicsFeatures, lastDiagnosticColumn, keep_feature_name_column=False)

    else:
        originalRadiomics = dfPyradiomicsFeatures.filter(regex=r'^original_*')
        if not originalRadiomics.empty:
            # Get the first original feature column name - the features begin in this column
            firstOriginalFeature = originalRadiomics.columns[0]
            # Drop all the columns before the features start
            featsOnlyRadiomics = dropUpToFeature(dfPyradiomicsFeatures, firstOriginalFeature, keep_feature_name=True)
        else:
            raise ValueError("PyRadiomics file doesn't contain any diagnostics or original feature columns, so can't find beginning of features. Use dropUpToFeature and specify the last non-feature or first PyRadiomic feature column name to get only PyRadiomics features.")

    return featsOnlyRadiomics


def dropUpToFeature(dataframe:DataFrame,
                    feature_name:str,
                    keep_feature_name_column:Optional[bool] = False
                    ):
    """ Function to drop all columns up to and possibly including the specified feature.

    Parameters
    ----------
    dataframe : DataFrame
        Dataframe to drop columns from.
    feature_name : str
        Name of the feature to drop up to.
    keep_feature_name_column : bool, optional
        Whether to keep the specified feature name column in the dataframe or drop it. The default is False.
        
    Returns
    -------
    dataframe : DataFrame
        Dataframe with all columns up to and including the specified feature dropped.
    """
    try:
        if keep_feature_name_column:
            # Get the column names up to but not including the specified feature
            column_names = dataframe.columns.to_list()[:dataframe.columns.get_loc(feature_name)]
        else:
            # Get the column names up to and including the specified feature
            column_names = dataframe.columns.to_list()[:dataframe.columns.get_loc(feature_name)+1]

        # Drop all columns up to and including the specified feature
        dataframe_dropped_columns = dataframe.drop(columns=column_names)

        return dataframe_dropped_columns
    
    except KeyError:
        print(f"Feature {feature_name} was not found as a column in dataframe. No columns dropped.")
        return dataframe
    
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def splitDataSetup(dfClinical:DataFrame,
                   dfFeatures:DataFrame,
                   splitVariables:dict,
                   imputeValue:Optional[str] = None,
                   ):
    """Function to split clinical and feature data tables into subgroups (e.g. train and test) based on clinical variable.
       

    Parameters
    ----------
    dfClinical : DataFrame
        Dataframe containing clinical data. Must include a column with labels from splitVariables.
        The index of this dataframe must correspond with dfFeatures' index.

    dfFeatures : DataFrame
        Dataframe containing feature data. This table will be split based on the patient index of the clinical table.
        The index of this dataframe must correspond with dfClinical's index.

    splitVariables : dict
        Dictionary where keys are labels of columns in dfClinical to base the split upon and values is a list of values to select
        in each data split.

    Returns
    -------
    splitClinicalDataframes : dict
        Dictionary of split up clinical dataframes. Keys indicate the value used to select the subset, value is the Dataframe
    
    splitFeatureDataframes : dict
        Dictionary of split up feature dataframes. Keys indicate the value used to select the clinical subset, value is the Dataframe

    Examples
    --------
    >>> splitDataSetup(clinicalData, featureData, {'data_split': ['train', 'test']})
    """

    # Initialize dictionaries to contain the split up dataframes
    splitClinicalDataframes = {}
    splitFeatureDataframes = {}

    # Get the clinical variable to split by and the list of values to select in each subgroup
    for variable, values in splitVariables.items():
        # Impute rows in the splitVariable that don't have one of the specified values
        desired_values_string = "|".join(values)
        # Create a regex to match any row that doesn't have one of the desired values
        impute_regex = r"^(?!" + desired_values_string + ").*"
        replacement_regex = re.compile(impute_regex, re.IGNORECASE)
        updatedSplitVariable = dfClinical[variable].replace(regex=replacement_regex, value=imputeValue)

        variable = f"{variable}_imputed"
        # Update the variable column with the imputed values
        dfClinical[variable] = updatedSplitVariable
        print(f"Made copy of split variable with imputed columns: {variable}")

        print(f"Getting split for {variable}")

        for value in values:
            # Get the clinical subset for rows with value in the variable column
            splitClinical = subsetDataframe(dataframe = dfClinical,
                                            includeDict = {variable: [value]})
            # Save this subset dataframe into the output dictionary
            splitClinicalDataframes[value] = splitClinical

            # Get the feature subset based on the indices from splitClinical
            splitFeature = subsetDataframe(dataframe = dfFeatures,
                                           includeDict={"index": splitClinical.index})
            # Save this subset dataframe into the output dictionary
            splitFeatureDataframes[value] = splitFeature
    
    return splitClinicalDataframes, splitFeatureDataframes


def savePlotFigure(sns_plot:heatmap,
                   plot_name:str,
                   output_dir_path:Optional[str]="",
                   dataset_name:Optional[str]="",):
    """Function to save out a seaborn plot to a png file.

    Parameters
    ----------
    sns_plot : seaborn.heatmap
        Seaborn plot to save out.
    plot_name : str
        What to name the plot on save. Ex. "RADCURE_original_vs_shuffled_correlation_plot.png"
    output_dir_path : str, optional
        Path to the directory to save the plot to. The default is "../../results".
    dataset_name : str, optional
        Name of the dataset to save the plot for. The default is "".
    """

    # Set default output directory path
    if not output_dir_path:
        output_dir_path = os.path.join("../../results", dataset_name, "plot_figures")

    # Setup output path
    output_path = os.path.join(output_dir_path, plot_name)

    # Make directory if it doesn't exist, but don't fail if it already exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Save out the plot
    sns_plot.get_figure().savefig(output_path, bbox_inches='tight')
    # print(f"Saved out plot to {output_path}")

    return



def makeProcessedDataFolders(dataset_name:str,
                             proc_data_path:str,
                             data_sources:Optional[Union[str, list]] = [""],
                             data_types:Optional[Union[str,list]] = [""],
                             train_test_split:bool = False,
                             ):
    """ Function to make the processed data folders for a dataset.

    Parameters
    ----------
    dataset_name : str
        Name of the dataset to make the processed data folders for.
    proc_data_path : str
        Path to the processed data folder that all folders will be made under.
    data_sources : list
        List of data sources to make the processed data folders for.
    data_types : list
        List of data types to make the processed data folders for.
    train_test_split : bool, optional
        Whether to make the train and test processed data folders. The default is False.
    """
    # Set up data sources and data types as lists if they are not already
    if type(data_sources) is str:
        data_sources = [data_sources]
    if type(data_types) is str:
        data_types = [data_types]

    # Make clinical procdata folder
    path_to_proc_clinical = os.path.join(proc_data_path, dataset_name, 'clinical')
    if not os.path.exists(path_to_proc_clinical):
        os.makedirs(path_to_proc_clinical)
    
    # Make feature procdata folders
    # Will combine each data source and data type into a list of strings
    for combo in itertools.product([dataset_name], data_sources, data_types):
        path_to_create = os.path.join(proc_data_path, *combo)
        
        if not os.path.exists(path_to_create):
            os.makedirs(path_to_create)
    
    # Optionally make train and test output folders
    if train_test_split is True:
        split_data_types = ['clinical', 'train_features', 'test_features']
        for combo in itertools.product([dataset_name], data_sources, ['train_test_split'], split_data_types):
            path_to_create = os.path.join(proc_data_path, *combo)
            
            if not os.path.exists(path_to_create):
                os.makedirs(path_to_create)
    
    return