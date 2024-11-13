import numpy as np
import pandas as pd

from pandas import DataFrame


def convertDaysToYears(dataframe_with_outcome:DataFrame,
                   time_column_label:str,
                   divide_by:int = 365):
    """ Function to create a copy of a time outcome column mesaured in days and convert it to years.

    Parameters
    ----------
    dataframe_with_outcome : DataFrame
        Dataframe containing the outcome column to convert.
    time_column_label : str
        Label for the time column to convert in the dataframe.
    divide_by : int, optional
        Value to divide the time column by. The default is 365.

    Returns
    -------
    dataframe_with_outcome : DataFrame
        Dataframe with a copy of the specified time column converted to years.
    """

    print(f"Converting time column {time_column_label} from days to years, using divide_by={divide_by}.")

    years_column_label = time_column_label + "_years"
    # Make a copy of the time column with the values converted from days to years and add suffic _years to the column name
    dataframe_with_outcome[years_column_label] = dataframe_with_outcome[time_column_label] / divide_by

    return dataframe_with_outcome


def timeOutcomeColumnSetup(dataframe_with_outcome:DataFrame,
                            outcome_column_label:str,
                            standard_column_label:str,
                            convert_to_years:bool = False,
                            ):
    """ Function to set up a time outcome column in a dataframe.

    Parameters
    ----------
    dataframe_with_outcome : DataFrame
        Dataframe containing the outcome column to convert.
    outcome_column_label : str
        Label for the outcome column to convert in the dataframe.
    standard_column_label : str
        Name of the column to save the standardized outcome column as.
    convert_to_years : bool, optional
        Whether to convert the time column to years. The default is False.

    Returns
    -------
    dataframe_with_outcome : DataFrame
        Dataframe with a copy of the specified outcome column converted to years.    
    """

    # Check if the outcome column is numeric
    if not np.issubdtype(dataframe_with_outcome[outcome_column_label].dtype, np.number):
        raise ValueError(f"{outcome_column_label} is not numeric. Please confirm outcome_column_label is the correct column or convert the column in the dataframe to numeric.")
    else:
        print(f"{outcome_column_label} is numeric. Making copy with standardized column name.")
        dataframe_with_standardized_outcome = dataframe_with_outcome.copy()

        if convert_to_years:
            print(f"Making copy of {outcome_column_label} with standardized column name and converting to years.")
            dataframe_with_standardized_outcome = convertDaysToYears(dataframe_with_standardized_outcome, outcome_column_label)

            # Rename the converted column with the standardized name
            dataframe_with_standardized_outcome.rename(columns={f"{outcome_column_label}_years": standard_column_label}, inplace=True)
        else:
            print(f"Making copy of {outcome_column_label} with standardized column name")
            # Make a copy of the outcome column with the standardized name
            dataframe_with_standardized_outcome[standard_column_label] = dataframe_with_standardized_outcome[outcome_column_label]
    
    return dataframe_with_outcome


def eventOutcomeColumnSetup(dataframe_with_outcome:DataFrame,
                            outcome_column_label:str,
                            standard_column_label:str,
                            event_column_values:dict):
    """ Function to set up an event outcome column in a dataframe.

    Parameters
    ----------
    dataframe_with_outcome : DataFrame
        Dataframe containing the outcome column to convert.
    outcome_column_label : str
        Label for the outcome column to convert in the dataframe.
    standard_column_label : str
        Name of the column to save the standardized outcome column as.
    event_column_values : dict
        Dictionary of event values to convert to numeric values. Keys are the event values, values are the numeric values.
    
    Returns
        
    """