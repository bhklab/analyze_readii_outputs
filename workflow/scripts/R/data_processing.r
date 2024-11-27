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

#' Function to make a report for a trained CPH model. If no inputs are provided, will provide an empty version
#' of the report to use for initialization. 
#' 
#' @param model_feature_weights A named list of model weights for the trained CPH model.
#' @param performance_results A named list of performance results for the trained CPH model (e.g. output from concordance.index).
#' 
#' @return A named list containing the model weights, and concordance index, confidence intervals, and p-value from the performance results.
makeCPHModelReport <- function(model_feature_weights = list(), 
                               performance_results = list()) {

    # Initialize empty CPH model report - this will be what is returned if no inputs are provided
    model_data <- c(features = list(), ci = 0, conf_lower = 0, conf_upper = 0, pval = 0)

    # Store the model feature weights if provided
    if (length(model_feature_weights) > 0) { model_data$features = model_feature_weights }    
    
    # Store a subset of the performance results if provided
    if (length(performance_results) > 0) {
        model_data$ci = performance_results$c.index,
        model_data$conf_lower = performance_results$lower,
        model_data$conf_upper = performance_results$upper,
        model_data$pval = performance_results$p.value
    }

    return(model_data)
}