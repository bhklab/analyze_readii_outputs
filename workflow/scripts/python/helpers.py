import numpy as np 
import os
import pandas as pd 
import re

from pandas import Series, DataFrame
from collections.abc import Sequence, Mapping

from typing import Optional, Union


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
        print("Multiple patient identifier labels found. Using the first one.")
    
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