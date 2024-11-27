
#### PSEUDO-CODE ####
# solution_length_options = [4, 10, 50, 100]
# solution_count = 15
#
#
# for solution_length in solution_length_options:
#     best_c_index = 0
#     selected_features_and_weights = []

#     run mrmr ensemble 

#     for solution in solution_count:
#         fit cph model to solution
#         evaluate cph model 
#         if better than best_c_index:
#             best_c_index = solution c_index
#             selected_features_and_weights = solution features and weights

# install.packages("BiocManager", repos = "http://cran.us.r-project.org")
# install.packages("mRMRe", repos = "http://cran.us.r-project.org")
# install.packages("checkmate", repos = "http://cran.us.r-project.org")
# BiocManager::install("survcomp")

source("workflow/scripts/R/data_processing.r")


#' Function to run mRMRe classic method on a set of features
#' 
#' @param data A data.frame containing the features to run mRMR on.
#' @param n_features Number of features to use for mRMR. Default is 30.
#' 
#' @return A vector of indices for the features selected by mRMR.
runMRMRClassic <- function(feature_data, 
                           n_features = 30) {
    # setup data for the mRMRe function
    mrmr_feature_data <- mRMR.data(data=data)

    # run MRMR with classic method 
    mrmr_results <- mRMR.classic(data = mrmr_feature_data,   
                                 target_indices = c(1),
                                 feature_count = n_features)

    # Extract a vector contain n_feature feature indices into feature_data
    fmcib_indices <- solutions(fmcib_results)[[1]]

    return(fmcib_indices)
}


#' Function to run mRMRe bootstrap ensemble method on a set of features
#' 
#' @param feature_data A data.frame containing the features to run mRMR on.
#' @param n_features Number of features to use for mRMR. Default is 30.
runMRMRBootstrap <- function(feature_data, n_features, n_solutions) {
    # setup data for the mRMRe function
    mrmr_feature_data <- mRMR.data(data = feature_data)

    # run mRMR with bootstrap ensemble method
    mrmr_results <- mRMR.ensemble(data = mrmr_feature_data,
                                    target_indices = c(1),
                                    solution_count = n_solutions,
                                    feature_count = n_features)
    
    # Extract a vector containing n_solutions lists of n_features feature indices into feature_data
    # rows = selected feature indices, columns = solution_number
    mrmr_solution_feature_indices <- solutions(mrmr_results)[[1]]

    return(mrmr_solution_feature_indices)
}


trainMRMRCoxModel <- function(labelled_train_data,
                              n_features,
                              n_solutions,
                              mrmr_method='bootstrap',
                              surv_time_label="survival_time_in_years",
                              surv_event_label="survival_event_binary") {

    # Get the feature data for the training set with outcome labels removed
    train_feature_data <- dropLabelsFromFeatureData(labelled_train_data, labels_to_drop=c("patient_ID", surv_time_label, surv_event_label))
    
    # Initialize best solution holder variables
    best_c_index = 0
    best_features_and_weights = list()

    # Perform mRMR feature selection using the specified method
    if (mrmr_method == 'bootstrap') {
        # Will return a n_features x n_solutions matrix with feature indices into train_data as values
        all_mrmr_solutions_matrix <- runMRMRBootstrap(train_feature_data,
                                        n_features,
                                        n_solutions)
    } else { # Have this here so future methods can be added (e.g. exhaustive mRMRe)
        print('Error: Invalid mRMR method')
        return(NULL)
    }

    # Loop through solutions to fit CPH models and determine best one
    for (solution_idx in 1:n_solutions) {
        print(paste("Solution", solution_idx))
        # Get the solution, will be a vector of length n_features
        solution_feature_indices <- all_mrmr_solutions_matrix[,solution_idx]

        # # Get the feature names for the solution to pass to the Cox model
        solution_feature_names <- names(train_feature_data)[solution_feature_indices]

        # Train the Cox model with the solution features
        solution_cph_coefficients <- trainCoxModel(labelled_train_data,
                                                    surv_time_label = surv_time_label,
                                                    surv_event_label = surv_event_label,
                                                    model_feature_list = solution_feature_names)
    
        
    }
}
