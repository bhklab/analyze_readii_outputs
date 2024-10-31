import numpy as np 
import os
import pandas as pd 
import re
import yaml

from pandas import Series, DataFrame
from seaborn import heatmap
from collections.abc import Sequence, Mapping

from typing import Optional, Union, Dict

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
        Dictionary of column names and values to include in the subset.
    excludeDict : dict
        Dictionary of column names and values to exclude from the subset.

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
        # Get the last diagnostics column index - the features begin in the next column
        lastDiagnosticIdx = dfPyradiomicsFeatures.columns.get_loc(diagnosticRadiomics.columns[-1])
        # Drop all the columns before the features start
        featsOnlyRadiomics = dfPyradiomicsFeatures.iloc[:, lastDiagnosticIdx+1:]

    else:
        originalRadiomics = dfPyradiomicsFeatures.filter(regex=r'^original_*')
        if not originalRadiomics.empty:
            # Get the first original feature column index - the features begin in this column
            firstOriginalIdx = dfPyradiomicsFeatures.columns.get_loc(originalRadiomics.columns[0])
            # Drop all the columns before the features start
            featsOnlyRadiomics = dfPyradiomicsFeatures.iloc[:, firstOriginalIdx:]
        else:
            raise ValueError("PyRadiomics file doesn't contain any diagnostics or original feature columns, so can't find beginning of features.")

    return featsOnlyRadiomics


def survivalTimeColumnSetup(clinical_dataframe:pd.DataFrame,
                            time_column_label:str,
                            output_time_column_label:Optional[str] = "survival_time_in_years",
                            convert_to_years:Optional[bool] = False,
                            divide_by:Optional[int] = 365,
                            dataset_name:Optional[str] = "config"):
    """ 
    Function to set up the survival time column in a clinical dataset for use in a predictive model (e.g. Cox PH)

    Parameters
    ----------
    clinical_dataframe:pd.DataFrame
    time_column_label:str
    convert_to_years:Optional[bool] = False
    divide_by:Optional[int] = 365
    dataset_name:Optional[str] = "config"
    """    

    # Confirm that time column is numeric
    if not np.issubdtype(clinical_dataframe[time_column_label].dtype, np.number):
        raise ValueError(f"Time column {time_column_label} is not numeric. Please confirm time label in {dataset_name}.yaml is the correct column or convert to numeric.")
    else:
        print(f"Time column {time_column_label} is numeric. Making copy with standardized column name.")
        if convert_to_years:
            print(f"Converting time column {time_column_label} from days to years.")
            clinical_dataframe["survival_time_in_years"] = clinical_dataframe[time_column_label] / divide_by
        else:
            clinical_dataframe["survival_time_in_years"] = clinical_dataframe[time_column_label]

    return clinical_dataframe


def survivalEventColumnSetup(clinical_dataframe:pd.DataFrame,
                            event_column_label:str,
                            output_event_column_label:Optional[str] = "survival_event_binary",
                            dataset_name:Optional[str] = "config"):
    
    # Determine what type of value is in event column
    event_variable_type = type(clinical_dataframe[event_column_label][0])
    if np.issubdtype(event_variable_type, np.number):
        print(f"Event column {event_column_label} is binary. Making copy with outpute event column name.")
        clinical_dataframe[output_event_column_label] = clinical_dataframe[event_column_label]

    elif np.issubdtype(event_variable_type, np.bool_):
        print(f"Event column {event_column_label} is boolean. Converting to binary and making copy with output event column name.")
        clinical_dataframe[output_event_column_label] = clinical_dataframe[event_column_label].astype(int)

    elif np.issubdtype(event_variable_type, np.str_):
        print(f"Event column {event_column_label} is string. Checking what values are present in the column.")

        event_column_values = clinical_dataframe[event_column_label].str.lower().unique()

        if len(event_column_values) != 2:
            raise ValueError(f"Event column {event_column_label} can only have two values. Please confirm event label in {dataset_name}.yaml is the correct column or update to have only two values.")
        
        # Check if alive and dead are present in the event column
        if 'alive' in event_column_values and 'dead' in event_column_values: 
            print(f"Converting to binary where 0 is alive and 1 is dead and making copy with output event column name.")

            clinical_dataframe[output_event_column_label] = clinical_dataframe[event_column_label].str.lower().replace({'alive': '0', 'dead': '1'}).astype(int)

        else:
            raise ValueError(f"Event column {event_column_label} doesn't contain any variation of 'alive' and 'dead'. Please confirm event label in {dataset_name}.yaml is the correct column.")
        
    return clinical_dataframe



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

        print("Getting split for ", variable)

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
    sns_plot.get_figure().savefig(output_path)
    print(f"Saved out plot to {output_path}")

    return

