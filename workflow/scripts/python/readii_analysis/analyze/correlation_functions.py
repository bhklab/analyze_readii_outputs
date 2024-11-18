import pandas as pd
from typing import Optional
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.linalg import issymmetric

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
    corr_fig : matplotlib.pyplot.figure
        Figure object containing a Seaborn heatmap.
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
    
    plt.title(subtitle, fontsize=12)
    plt.suptitle(title, fontsize=14)
    
    return corr_fig


def plotCorrelationDistribution(correlation_matrix:pd.DataFrame,
                                num_bins:Optional[int] = 100,
                                xlabel:Optional[str] = "Correlations",
                                ylabel:Optional[str] = "Frequency",
                                y_lower_bound:Optional[int] = 0,
                                y_upper_bound:Optional[int] = None,
                                title:Optional[str] = "Distribution of Correlations for Features",
                                subtitle:Optional[str] = "",
                                ):
    """ Function to plot a distribution of correlation values for a correlation matrix.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    num_bins : int, optional
        Number of bins to use for the distribution plot. The default is 100.
    xlabel : str, optional
        Label for the x-axis. The default is "Correlations".
    ylabel : str, optional
        Label for the y-axis. The default is "Frequency".
    y_lower_bound : int, optional
        Lower bound for the y-axis of the distribution plot. The default is 0.
    y_upper_bound : int, optional
        Upper bound for the y-axis of the distribution plot. The default is None.
    title : str, optional
        Title for the plot. The default is "Distribution of Correlations for Features".
    subtitle : str, optional
        Subtitle for the plot. The default is "".

    Returns
    -------
    dist_fig : plt.Figure
        Figure object containing the histogram of correlation values.
    bin_values : np.ndarray or list of arrays
        Numpy array containing the values in each bin for the histogram.
    bin_edges : np.ndarray
        Numpy array containing the bin edges for the histogram.
    """
    
    # Convert to numpy to use histogram function
    feature_correlation_arr = correlation_matrix.to_numpy()

    # Check if matrix is symmetric
    if issymmetric(feature_correlation_arr):
        print("Correlation matrix is symmetric.")
        # Get only the bottom left triangle of the correlation matrix since the matrix is symmetric 
        lower_half_idx = np.mask_indices(feature_correlation_arr.shape[0], np.tril)
        # This is a 1D array for binning and plotting
        correlation_vals = feature_correlation_arr[lower_half_idx]
    else:
        # Flatten the matrix to a 1D array for binning and plotting
        correlation_vals = feature_correlation_arr.flatten()

    dist_fig, dist_ax = plt.subplots()
    bin_values, bin_edges, _ = dist_ax.hist(correlation_vals, bins=num_bins)
    dist_ax.set_xlabel(xlabel)
    dist_ax.set_ylabel(ylabel)
    dist_ax.set_ybound(y_lower_bound, y_upper_bound)
    plt.suptitle(title, fontsize=14)
    plt.title(subtitle, fontsize=10)

    return dist_fig, bin_values, bin_edges


def getVerticalSelfCorrelations(correlation_matrix:pd.DataFrame,
                                num_vertical_features:int):
    """ Function to get the vertical (y-axis) self correlations from a correlation matrix. Gets the top left quadrant of the correlation matrix.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to get the vertical self correlations from.
    num_vertical_features : int
        Number of vertical features in the correlation matrix.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the vertical self correlations from the correlation matrix.    
    """
    if num_vertical_features > correlation_matrix.shape[0]:
        raise ValueError(f"Number of vertical features ({num_vertical_features}) is greater than the number of rows in the correlation matrix ({correlation_matrix.shape[0]}).")
    
    if num_vertical_features > correlation_matrix.shape[1]:
        raise ValueError(f"Number of vertical features ({num_vertical_features}) is greater than the number of columns in the correlation matrix ({correlation_matrix.shape[1]}).")

    # Get the correlation matrix for vertical vs vertical - this is the top left corner of the matrix
    return correlation_matrix.iloc[0:num_vertical_features, 0:num_vertical_features]



def getHorizontalSelfCorrelations(correlation_matrix:pd.DataFrame,
                                  num_horizontal_features:int):
    """ Function to get the horizontal (x-axis) self correlations from a correlation matrix. Gets the bottom right quadrant of the correlation matrix.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to get the horizontal self correlations from.
    num_horizontal_features : int
        Number of horizontal features in the correlation matrix.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the horizontal self correlations from the correlation matrix.
    """
    
    if num_horizontal_features > correlation_matrix.shape[0]:
        raise ValueError(f"Number of horizontal features ({num_horizontal_features}) is greater than the number of rows in the correlation matrix ({correlation_matrix.shape[0]}).")
    
    if num_horizontal_features > correlation_matrix.shape[1]:
        raise ValueError(f"Number of horizontal features ({num_horizontal_features}) is greater than the number of columns in the correlation matrix ({correlation_matrix.shape[1]}).")

    # Get the index of the start of the horizontal correlations
    start_of_horizontal_correlations = len(correlation_matrix.columns) - num_horizontal_features

    # Get the correlation matrix for horizontal vs horizontal - this is the bottom right corner of the matrix
    return correlation_matrix.iloc[start_of_horizontal_correlations:, start_of_horizontal_correlations:]


def getCrossCorrelationMatrix(correlation_matrix:pd.DataFrame,
                              num_vertical_features:int):
    """ Function to get the cross correlation matrix subsection for a correlation matrix. Gets the top right quadrant of the correlation matrix so vertical and horizontal features are correctly labeled.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to get the cross correlation matrix subsection from.
    num_vertical_features : int
        Number of vertical features in the correlation matrix.
    
    Returns
    -------
    pd.DataFrame
        Dataframe containing the cross correlations from the correlation matrix.
    """

    if num_vertical_features > correlation_matrix.shape[0]:
        raise ValueError(f"Number of vertical features ({num_vertical_features}) is greater than the number of rows in the correlation matrix ({correlation_matrix.shape[0]}).")
    
    if num_vertical_features > correlation_matrix.shape[1]:
        raise ValueError(f"Number of vertical features ({num_vertical_features}) is greater than the number of columns in the correlation matrix ({correlation_matrix.shape[1]}).")
    
    return correlation_matrix.iloc[0:num_vertical_features, num_vertical_features:]



def plotSelfCorrelationHeatMaps(correlation_matrix:pd.DataFrame,
                                axis:str,
                                num_axis_features:int,
                                feature_name:str,
                                correlation_method:Optional[str] = "",
                                extraction_method:Optional[str] = "",
                                dataset_name:Optional[str] = "",
                                cmap:Optional[str] = "nipy_spectral",
                                ):
    """ Function to plot a correlation heatmap for the vertical (y-axis) and horizontal (x-axis) self correlations.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    axis : str
        Axis to plot the self correlations for. Must be either "vertical" or "horizontal". The default is "vertical".
    num_axis_features : int
        Number of features in the axis to plot the self correlations for. This is used to get the self correlations from the correlation matrix.
    feature_name : str
        Name of the feature to use for the plot title and subtitle.
    correlation_method : str, optional
        Name of the correlation method to use for the plot title and subtitle. The default is "".
    extraction_method : str, optional
        Name of the extraction method to use for the plot title and subtitle. The default is "".
    dataset_name : str, optional
        Name of the dataset to use for the plot title and subtitle. The default is "".
    cmap : str, optional
        Name of the matplotlib colormap to use for the heatmap. The default is "nipy_spectral".

    Returns
    -------
    self_plot : matplotlib.pyplot.figure
        Figure object containing a Seaborn heatmap of the vertical or horizontalself correlations from the correlation matrix.
    """

    if axis == "vertical":
        # Get the correlation matrix for vertical vs vertical
        # This is the top left corner of the matrix
        self_correlations = getVerticalSelfCorrelations(correlation_matrix, num_axis_features)
    elif axis == "horizontal":
        # Get the correlation matrix for horizontal vs horizontal
        # This is the bottom right corner of the matrix
        self_correlations = getHorizontalSelfCorrelations(correlation_matrix, num_axis_features)
    else:
        raise ValueError(f"Axis must be either 'vertical' or 'horizontal'. Provided axis is {axis}.")

    # Create correlation heatmap for vertical vs vertical
    self_plot = plotCorrelationHeatmap(self_correlations,
                                       diagonal = True,
                                       triangle = "lower",
                                       cmap = cmap,
                                       xlabel = feature_name,
                                       ylabel = feature_name,
                                       title = f"{correlation_method.capitalize()} Self Correlations for {dataset_name} {extraction_method.capitalize()} Features",
                                       subtitle = f"{feature_name} vs. {feature_name}")

    return self_plot



def plotCrossCorrelationHeatmap(correlation_matrix:pd.DataFrame,
                                num_vertical_features:int,
                                vertical_feature_name:str,
                                horizontal_feature_name:str,
                                correlation_method:Optional[str] = "",
                                extraction_method:Optional[str] = "",
                                dataset_name:Optional[str] = "",
                                cmap:Optional[str] = "nipy_spectral"):
    """ Function to plot heatmap for a the cross-correlation section of a correlation matrix. Will be the top right quadrant of the correlation matrix so vertical and horizontal features are correctly labeled.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    num_vertical_features : int
        Number of vertical (y-axis) features in the correlation matrix.
        The number of vertical features must be less than the number of rows in the correlation matrix.
    vertical_feature_name : str
        Name of the vertical feature to use for the plot title and subtitle.
    horizontal_feature_name : str
        Name of the horizontal feature to use for the plot title and subtitle.
    correlation_method : str, optional
        Name of the correlation method to use for the plot title and subtitle. The default is "".
    extraction_method : str, optional
        Name of the extraction method to use for the plot title and subtitle. The default is "".
    dataset_name : str, optional
        Name of the dataset to use for the plot title and subtitle. The default is "".
    cmap : str, optional
        Name of the matplotlib colormap to use for the heatmap. The default is "nipy_spectral".
    
    Returns
    -------
    cross_corr_plot : matplotlib.pyplot.figure
        Figure object containing a Seaborn heatmap of the cross correlations from the correlation matrix.
    """
    
    # Get the cross correlation matrix from the main correlation matrix
    cross_corr_matrix = getCrossCorrelationMatrix(correlation_matrix, num_vertical_features)

    # Create heatmap for the cross correlation matrix
    cross_corr_plot = plotCorrelationHeatmap(cross_corr_matrix,
                                             diagonal = False,
                                             cmap=cmap,
                                             xlabel = vertical_feature_name,
                                             ylabel = horizontal_feature_name,
                                             title = f"{correlation_method.capitalize()} Cross Correlations for {dataset_name} {extraction_method.capitalize()} Features",
                                             subtitle = f"{vertical_feature_name} vs. {horizontal_feature_name}")

    return cross_corr_plot



def plotSelfCorrelationDistributionPlots(correlation_matrix:pd.DataFrame,
                                         axis:str,
                                         num_axis_features:int,
                                         feature_name:str,
                                         num_bins: Optional[int] = 450,
                                         y_upper_bound:Optional[int] = None,
                                         correlation_method:Optional[str] = "",
                                         extraction_method:Optional[str] = "",
                                         dataset_name:Optional[str] = "",
                                         ):
    """ Function to plot a distribution of self correlation values for a correlation matrix.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    axis : str
        Axis to plot the self correlations for. Must be either "vertical" or "horizontal".
    num_axis_features : int
        Number of features in the axis to plot the self correlations for. This is used to get the self correlations from the correlation matrix.
    feature_name : str
        Name of the feature to use for the plot title and subtitle.
    num_bins : int, optional
        Number of bins to use for the distribution plot. The default is 450.
    y_upper_bound : int, optional
        Upper bound for the y-axis of the distribution plot. The default is None.
    correlation_method : str, optional
        Name of the correlation method to use for the plot title and subtitle. The default is "".
    extraction_method : str, optional
        Name of the extraction method to use for the plot title and subtitle. The default is "".
    dataset_name : str, optional
        Name of the dataset to use for the plot title and subtitle. The default is "".

    Returns
    -------
    self_corr_dist_fig : plt.Figure
        Figure object containing the histogram of self correlation values.
    """
    
    if axis == "vertical":
        # Get the correlation matrix for vertical vs vertical
        # This is the top left corner of the matrix
        self_correlations = getVerticalSelfCorrelations(correlation_matrix, num_axis_features)
    elif axis == "horizontal":
        # Get the correlation matrix for horizontal vs horizontal
        # This is the bottom right corner of the matrix
        self_correlations = getHorizontalSelfCorrelations(correlation_matrix, num_axis_features)
    else:
        raise ValueError(f"Axis must be either 'vertical' or 'horizontal'. Provided axis is {axis}.")
    
    # Plot the distribution of correlation values for the self correlations
    self_corr_dist_fig, _, _ = plotCorrelationDistribution(self_correlations,
                                                               num_bins = num_bins,
                                                               xlabel = f"{correlation_method.capitalize()} Correlation",
                                                               ylabel = "Frequency",
                                                               y_upper_bound=y_upper_bound,
                                                               title = f"Distribution of {correlation_method.capitalize()} Self Correlations for {dataset_name} {extraction_method.capitalize()} Features",
                                                               subtitle = f"{feature_name} vs. {feature_name}"
                                                               )
                                                                                                     
    return self_corr_dist_fig


def plotCrossCorrelationDistributionPlots(correlation_matrix:pd.DataFrame,
                                          num_vertical_features:int,
                                          vertical_feature_name:str,
                                          horizontal_feature_name:str,
                                          num_bins: Optional[int] = 450,
                                          y_upper_bound:Optional[int] = None,
                                          correlation_method:Optional[str] = "",
                                          extraction_method:Optional[str] = "",
                                          dataset_name:Optional[str] = "",
                                          ):
    """ Function to plot a distribution of cross correlation values for a correlation matrix. Will be the top right quadrant of the correlation matrix so vertical and horizontal features are correctly labeled.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    num_vertical_features : int
        Number of vertical (y-axis) features in the correlation matrix.
        The number of vertical features must be less than the number of rows and columns in the correlation matrix.
    vertical_feature_name : str
        Name of the vertical feature to use for the plot title and subtitle.
    horizontal_feature_name : str
        Name of the horizontal feature to use for the plot title and subtitle.
    num_bins : int, optional
        Number of bins to use for the distribution plot. The default is 450.
    y_upper_bound : int, optional
        Upper bound for the y-axis of the distribution plot. The default is None.
    correlation_method : str, optional
        Name of the correlation method to use for the plot title and subtitle. The default is "".
    extraction_method : str, optional
        Name of the extraction method to use for the plot title and subtitle. The default is "".
    dataset_name : str, optional
        Name of the dataset to use for the plot title and subtitle. The default is "".

    Returns
    -------
    cross_corr_dist_fig : plt.Figure
        Figure object containing the histogram of cross correlation values.
    """
    
    # Get the cross correlation matrix from the main correlation matrix
    cross_corr_matrix = getCrossCorrelationMatrix(correlation_matrix, num_vertical_features)

    # Create heatmap for the cross correlation matrix
    cross_corr_dist_fig, _, _ = plotCorrelationDistribution(cross_corr_matrix,
                                                      num_bins = num_bins,
                                                      xlabel = f"{correlation_method.capitalize()} Correlation",
                                                      ylabel = "Frequency",
                                                      y_upper_bound = y_upper_bound,
                                                      title = f"Distribution of {correlation_method.capitalize()} Cross Correlations for {dataset_name} {extraction_method.capitalize()} Features",
                                                      subtitle = f"{vertical_feature_name} vs. {horizontal_feature_name}"
                                                      )

    return cross_corr_dist_fig



def makeBothSelfCorrelationPlots(correlation_matrix:pd.DataFrame,
                             axis:str,
                             num_axis_features:int,
                             feature_name:str,
                             corr_cmap:Optional[str] = "nipy_spectral",
                             dist_num_bins: Optional[int] = 450,
                             dist_y_upper_bound:Optional[int] = None,
                             correlation_method:Optional[str] = "",
                             extraction_method:Optional[str] = "",
                             dataset_name:Optional[str] = ""
                             ):
    """ Function to make both the correlation heatmap and distribution plots for a correlation matrix.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    axis : str
        Axis to plot the self correlations for. Must be either "vertical" or "horizontal". The default is "vertical".
    num_axis_features : int
        Number of features in the axis to plot the self correlations for. This is used to get the self correlations from the correlation matrix.
    feature_name : str
        Name of the feature to use for the plot title and subtitle.
    corr_cmap : str, optional
        Name of the matplotlib colormap to use for the heatmap. The default is "nipy_spectral".
    dist_num_bins : int, optional
        Number of bins to use for the distribution plot. The default is 450.
    dist_y_upper_bound : int, optional
        Upper bound for the y-axis of the distribution plot. The default is None.
    correlation_method : str, optional
        Name of the correlation method to use for the plot title and subtitle. The default is "".
    extraction_method : str, optional
        Name of the extraction method to use for the plot title and subtitle. The default is "".
    dataset_name : str, optional
        Name of the dataset to use for the plot title and subtitle. The default is "".

    Returns
    -------
    self_corr_plot : matplotlib.pyplot.figure
        Figure object containing a Seaborn heatmap of the self correlations from the correlation matrix.
    self_corr_dist_plot : matplotlib.pyplot.figure
        Figure object containing a histogram of the distribution of self correlations from the correlation matrix.
    """
    # Plot the correlation heatmap for the self correlations
    self_corr_plot = plotSelfCorrelationHeatMaps(correlation_matrix = correlation_matrix,
                                                 axis = axis,
                                                 num_axis_features = num_axis_features,
                                                 feature_name = feature_name,
                                                 cmap = corr_cmap,
                                                 correlation_method = correlation_method,
                                                 extraction_method = extraction_method,
                                                 dataset_name = dataset_name)

    # Plot the distribution of correlation values for the self correlations
    self_corr_dist_plot = plotSelfCorrelationDistributionPlots(correlation_matrix = correlation_matrix,
                                                               axis = axis,
                                                               num_axis_features = num_axis_features,
                                                               feature_name = feature_name,
                                                               num_bins = dist_num_bins,
                                                               y_upper_bound=dist_y_upper_bound,
                                                               correlation_method = correlation_method,
                                                               extraction_method = extraction_method,
                                                               dataset_name = dataset_name)

    return self_corr_plot, self_corr_dist_plot


def makeBothCrossCorrelationPlots(correlation_matrix:pd.DataFrame,
                                  num_vertical_features:int,
                                  vertical_feature_name:str,
                                  horizontal_feature_name:str,
                                  corr_cmap:Optional[str] = "nipy_spectral",
                                  dist_num_bins: Optional[int] = 450,
                                  dist_y_upper_bound:Optional[int] = None,
                                  correlation_method:Optional[str] = "",
                                  extraction_method:Optional[str] = "",
                                  dataset_name:Optional[str] = "",
                                  ):
    """ Function to make both the self correlation heatmap and distribution plots for a correlation matrix.

    Parameters
    ----------
    correlation_matrix : pd.DataFrame
        Dataframe containing the correlation matrix to plot.
    num_vertical_features : int
        Number of vertical (y-axis) features in the correlation matrix.
        The number of vertical features must be less than the number of rows and columns in the correlation matrix.
    vertical_feature_name : str
        Name of the vertical feature to use for the plot title and subtitle.
    horizontal_feature_name : str
        Name of the horizontal feature to use for the plot title and subtitle.
    corr_cmap : str, optional
        Name of the matplotlib colormap to use for the heatmap. The default is "nipy_spectral".
    dist_num_bins : int, optional
        Number of bins to use for the distribution plot. The default is 450.
    dist_y_upper_bound : int, optional
        Upper bound for the y-axis of the distribution plot. The default is None.
    correlation_method : str, optional
        Name of the correlation method to use for the plot title and subtitle. The default is "".
    extraction_method : str, optional
        Name of the extraction method to use for the plot title and subtitle. The default is "".
    dataset_name : str, optional
        Name of the dataset to use for the plot title and subtitle. The default is "".

    Returns
    -------
    cross_corr_plot : matplotlib.pyplot.figure
        Figure object containing a Seaborn heatmap of the cross correlations from the correlation matrix.
    cross_corr_dist_plot : matplotlib.pyplot.figure
        Figure object containing a histogram of the distribution of cross correlations from the correlation matrix.
    """
    # Plot the correlation heatmap for the cross correlations
    cross_corr_plot = plotCrossCorrelationHeatmap(correlation_matrix = correlation_matrix,
                                                  num_vertical_features = num_vertical_features,
                                                  vertical_feature_name = vertical_feature_name,
                                                  horizontal_feature_name = horizontal_feature_name,
                                                  cmap = corr_cmap,
                                                  correlation_method = correlation_method,
                                                  extraction_method = extraction_method,
                                                  dataset_name = dataset_name)
    
    # Plot the distribution of correlation values for the cross correlations
    cross_corr_dist_plot = plotCrossCorrelationDistributionPlots(correlation_matrix = correlation_matrix,
                                                                 num_vertical_features = num_vertical_features,
                                                                 vertical_feature_name = vertical_feature_name,
                                                                 horizontal_feature_name = horizontal_feature_name,
                                                                 num_bins = dist_num_bins,
                                                                 y_upper_bound = dist_y_upper_bound,
                                                                 correlation_method = correlation_method,
                                                                 extraction_method = extraction_method,
                                                                 dataset_name = dataset_name)

    return cross_corr_plot, cross_corr_dist_plot