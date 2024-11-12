from pandas import DataFrame
import pandas as pd

def getIntersectionDataframes(dataframe_A:DataFrame, 
                              dataframe_B:DataFrame): 
    """ Function to get the subset of two dataframes based on the intersection of their indices.

    Parameters
    ----------
    dataframe_A : DataFrame
        Dataframe A to get the intersection of based on the index.
    dataframe_B : DataFrame
        Dataframe B to get the intersection of based on the index.

    Returns
    -------
    intersection_index_dataframeA : DataFrame
        Dataframe containing the rows of dataframe A that are in the intersection of the indices of dataframe A and dataframe B.
    intersection_index_dataframeB : DataFrame
        Dataframe containing the rows of dataframe B that are in the intersection of the indices of dataframe A and dataframe B.
    """
    
    intersection_index = dataframe_A.index.intersection(dataframe_B.index)
    intersection_index_dataframeA = dataframe_A.loc[intersection_index]
    intersection_index_dataframeB = dataframe_B.loc[intersection_index]

    return intersection_index_dataframeA, intersection_index_dataframeB
