#' Function to drop columns by their label from a feature data frame.
#' 
#' @param labelled_feature_data A data frame containing the feature data to drop columns from.
#' @param labels_to_drop A vector of labels to drop from the feature data. Must be a subset of the column names in labelled_feature_data.
#' 
#' @return A data frame containing the feature data with the specified columns dropped.
dropLabelsFromFeatureData <- function(labelled_feature_data, labels_to_drop) {

    # Get the data column names of the labels to drop
    columns_to_drop <- tryCatch({ names(labelled_feature_data) %in% labels_to_drop
    }, error = function(e) {
        stop(paste("dropLabelsFromFeatureData: labels_to_drop must be a subset of the labels in labelled_feature_data."))
    })

    # Remove those columns from the feature data
    feature_data_only <- labelled_feature_data[!columns_to_drop]

    return(feature_data_only)
}