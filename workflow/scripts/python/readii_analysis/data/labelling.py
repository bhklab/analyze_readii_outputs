import numpy as np
import pandas as pd

from pandas import DataFrame
from typing import Optional


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
            print(f"Making copy of {outcome_column_label} with standardized column name and converting to years.\n")
            dataframe_with_standardized_outcome = convertDaysToYears(dataframe_with_standardized_outcome, outcome_column_label)

            # Rename the converted column with the standardized name
            dataframe_with_standardized_outcome.rename(columns={f"{outcome_column_label}_years": standard_column_label}, inplace=True)
        else:
            print(f"Making copy of {outcome_column_label} with standardized column name.\n")
            # Make a copy of the outcome column with the standardized name
            dataframe_with_standardized_outcome[standard_column_label] = dataframe_with_standardized_outcome[outcome_column_label]
    
    return dataframe_with_standardized_outcome


def eventOutcomeColumnSetup(dataframe_with_outcome:DataFrame,
                            outcome_column_label:str,
                            standard_column_label:str,
                            event_column_value_mapping:Optional[dict] = {}):
    """ Function to set up an event outcome column in a dataframe.

    Parameters
    ----------
    dataframe_with_outcome : DataFrame
        Dataframe containing the outcome column to convert.
    outcome_column_label : str
        Label for the outcome column to convert in the dataframe.
    standard_column_label : str
        Name of the column to save the standardized outcome column as.
    event_column_value_mapping : dict, optional
        Dictionary of event values to convert to numeric values. Keys are the event values, values are the numeric values. If provided, all event values in the outcome column must be handled by the dictionary.
        If not provided, will attempt to convert based on the values in the outcome column. By default alive and dead will be converted to 0 and 1, respectively, if present in the outcome column.
        The default is an empty dictionary.
    
    Returns
    -------
    dataframe_with_standardized_outcome : DataFrame
        Dataframe with a copy of the specified outcome column converted to numeric values.
    """
    # Get the type of the existing event column
    event_variable_type = type(dataframe_with_outcome[outcome_column_label][0])

    # Make a copy of the dataframe to work on 
    dataframe_with_standardized_outcome = dataframe_with_outcome.copy()

    # Handle numeric event column
    if np.issubdtype(event_variable_type, np.number):
        print(f"{outcome_column_label} is numeric. Making copy with standardized column name.\n")
        dataframe_with_standardized_outcome[standard_column_label] = dataframe_with_outcome[outcome_column_label]

        return dataframe_with_standardized_outcome

    # Handle boolean event column
    elif np.issubdtype(event_variable_type, np.bool_):
        print(f"{outcome_column_label} is boolean. Making copy with standardized column name and converting to numeric.\n")
        dataframe_with_standardized_outcome[standard_column_label] = dataframe_with_outcome[outcome_column_label].astype(int)

        return dataframe_with_standardized_outcome

    # Handle string event column
    elif np.issubdtype(event_variable_type, np.str_):
        print(f"{outcome_column_label} is string.")

        # Check if user provided a dictionary handles all event values in the outcome column
        if event_column_value_mapping and not all(value == map for value, map in zip(dataframe_with_outcome[outcome_column_label].str.lower().unique(), event_column_value_mapping.keys())):
            raise ValueError(f"Not all event values in {outcome_column_label} are handled by the provided event_column_value_mapping dictionary.")
        
            # TODO: add handling for values not in the dictionary

        # Create event column value mapping if dictionary is not provided
        elif not event_column_value_mapping: 
            # Get the existing event values in the provided dataframe and and sort them
            existing_event_values = sorted(dataframe_with_outcome[outcome_column_label].str.lower().unique())

            # Set the conversion value for the first event value to 0
            other_event_num_value = 0

            # See if alive is present, and set the conversion value for alive to 0
            if "alive" in existing_event_values:
                event_column_value_mapping['alive'] = 0
                # Remove alive from the list of existing event values
                existing_event_values.remove("alive")
                # Update starting value for other event values
                other_event_num_value = 1

            # See if dead is present, and set the conversion value for dead to 1
            if "dead" in existing_event_values:
                event_column_value_mapping['dead'] = 1
                # Remove dead from the list of existing event values
                existing_event_values.remove("dead")
                # Update starting value for other event values
                other_event_num_value = 2

            # Set the conversion value for the other event values to the next available value
            for other_event_value in existing_event_values:
                event_column_value_mapping[other_event_value] = other_event_num_value
                other_event_num_value += 1
        # end creating event column value mapping
        
        print(f"Converting values with mapping of: {event_column_value_mapping}.")

        # get the existing event values, make them lowercase, replace the dictionary values with the dictionary keys, convert to numeric, and save to the standardized column copy
        dataframe_with_standardized_outcome[standard_column_label] = dataframe_with_outcome[outcome_column_label].str.lower().replace(event_column_value_mapping).astype(int)
        
        return dataframe_with_standardized_outcome
    # end string handling                  

    else:
        raise TypeError("Event column {outcome_column_label} is not a valid type. Must be a string, boolean, or numeric.")
