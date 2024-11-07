import pandas as pd
from typing import Optional
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from matplotlib.colors import Colormap

def getFeatureCorrelations(vertical_features:pd.DataFrame,
                           horizontal_features:pd.DataFrame,
                           method:Optional[str] = "pearson",
                           vertical_feature_name:Optional[str] = "",
                           horizontal_feature_name:Optional[str] = ""):
    """ Function to calculate correlation between two sets of features.

    Parameters
    ----------
    vertical_features : pd.DataFrame
        Dataframe containing features to calculate correlations with.
    horizontal_features : pd.DataFrame
        Dataframe containing features to calculate correlations with.
    method : str
        Method to use for calculating correlations. Default is "pearson".
    vertical_feature_name : str
        Name of the vertical features to use as suffix in correlation dataframe. Default is blank "".
    horizontal_feature_name : str
        Name of the horizontal features to use as suffix in correlation dataframe. Default is blank "".
    
    Returns
    -------
    correlation_matrix : pd.DataFrame
        Dataframe containing correlation values.
    """
    # Join the features into one dataframe
    # Use inner join to keep only the rows that have a value in both vertical and horizontal features
    features_to_correlate = vertical_features.join(horizontal_features, 
                                                how='inner', 
                                                lsuffix=f"_{vertical_feature_name}", 
                                                rsuffix=f"_{horizontal_feature_name}") 

    # Calculate correlation between vertical features and horizontal features
    correlation_matrix = features_to_correlate.corr(method=method)

    return correlation_matrix


def plotCorrelationHeatmap(correlation_matrix_df:pd.DataFrame,
                           diagonal:Optional[bool] = False,
                           triangle:Optional[str] = "lower",
                           cmap:Optional[str] = "nipy_spectral",
                           xlabel:Optional[str] = "",
                           ylabel:Optional[str] = "",
                           title:Optional[str] = "",
                           subtitle:Optional[str] = "",
                           ):
    """Function to plot a correlation heatmap.

    Parameters
    ----------
    correlation_matrix_df : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    diagonal : bool, optional
        Whether to only plot half of the matrix. The default is False.
    triangle : str, optional
        Which triangle half of the matrixto plot. The default is "lower".
    xlabel : str, optional
        Label for the x-axis. The default is "".
    ylabel : str, optional
        Label for the y-axis. The default is "".
    title : str, optional
        Title for the plot. The default is "".
    subtitle : str, optional
        Subtitle for the plot. The default is "".

    Returns
    -------
    matplolib.pyplot.figure
    """
    if diagonal:
        if triangle == "lower":
            # Mask out the upper triangle half of the matrix
            mask = np.triu(correlation_matrix_df)
        elif triangle == "upper":
            # Mask out the lower triangle half of the matrix
            mask = np.tril(correlation_matrix_df)
        else:
            raise ValueError("If diagonal is True, triangle must be either 'lower' or 'upper'.")
    else:
        mask = None
    
    if not title:
        title = "Correlation Heatmap"


    corr_fig, corr_ax = plt.subplots()


    # Plot the correlation matrix
    corr_ax = sns.heatmap(correlation_matrix_df,
                         mask = mask,
                         cmap=cmap,
                         vmin=-1.0,
                         vmax=1.0)
    
    # Remove the individual feature names from the axes
    corr_ax.set_xticklabels(labels=[])
    corr_ax.set_yticklabels(labels=[])

    # Set axis labels
    corr_ax.set_xlabel(xlabel)
    corr_ax.set_ylabel(ylabel)
    
    plt.title(title, fontsize=14, y = 1.13)
    plt.suptitle(subtitle, fontsize=12)
    
    return corr_fig
