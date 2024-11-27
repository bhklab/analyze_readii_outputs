
# install.packages("BiocManager", repos = "http://cran.us.r-project.org")
# install.packages("mRMRe", repos = "http://cran.us.r-project.org")
# install.packages("checkmate", repos = "http://cran.us.r-project.org")
# install.packages("caret", repos = "http://cran.us.r-project.org")
# BiocManager::install("survcomp")

source("workflow/scripts/R/io.r")
source("workflow/scripts/R/data_processing.r")
source("workflow/scripts/R/mrmr_functions.r")

library(survival)
library(survcomp)
library(tools)
library(caret)
library(checkmate)
library(mRMRe)

# #' Function to train a CoxPH model to make a signature based on a set of features.
# #' Will return the trained weights for the model.
# #'
# #' @param train_data Data.frame containing the training features with outcome labels included.
# #' @param val_data Data.frame containing the validation features with outcome labels included.
# #' @param surv_time_label Label for the time column in the training features file.
# #' @param surv_event_label Label for the event column in the training features file.
# #' @param model_feature_list List of feature names to use for the Cox model.
# #' @param start Starting number of features to use for MRMR model training. Default is 2.
# #' @param end Ending number of features to use for MRMR model training. Default is 50.
# #' @param by Increment for number of features to use for MRMR model training. Default is 2.
# #' 
# #' @return vector of trained weights.
# trainMRMRCoxModel <- function(train_data,
#                               val_data,
#                               surv_time_label,
#                               surv_event_label,
#                               model_feature_list,
#                               start=2,
#                               end=50,
#                               by=2){
#     best_ci <- 0
#     for (k_mrmr in seq(start, end, by=by)) {
#             print(paste("Starting k_mrmr:", k_mrmr))
#             mrmr_idx = runMRMR(train_data[,model_feature_list], n_features=k_mrmr)
#             mrmr_feature_list = names(train_data[,model_feature_list])[mrmr_idx]

#             coefs <- trainCoxModel(train_data,
#                                    surv_time_label,
#                                    surv_event_label,
#                                    mrmr_feature_list)   
#                         val_ci <- testCoxModel(val_data,
#                                    surv_time_label,
#                                    surv_event_label,
#                                    mrmr_feature_list,
#                                    coefs)$c.index

#             if (val_ci > best_ci) {
#                 best_ci <- val_ci
#                 best_k <- k_mrmr``
#                 best_feats <- mrmr_feature_list
#                 best_coefs <- coefs
#             }
#         }
#     return(setNames(as.list(best_coefs), best_feats))
# }

#' Function to train a CoxPH model to make a signature based on a set of features.
#' Will return the trained weights for the model.
#' 
#' @param train_labelled_features_file_path Path to the file containing the training features with outcome labels included.
#' @param surv_time_label Label for the time column in the training features file.
#' @param surv_event_label Label for the event column in the training features file.
#' @param model_feature_list List of feature names to use for the Cox model.
#' 
#' @return vector of trained weights with feature names.
trainCoxModel <- function(labelled_feature_data,
                          surv_time_label,
                          surv_event_label,
                          model_feature_list){ #nolint
     # Get just features selected for the model
    train_feature_data <- tryCatch({
        labelled_feature_data[, model_feature_list]
    }, error = function(e) {
        stop(paste("trainCoxModel:Model features not found in provided feature set:", model_feature_list))
    })

    # Get the time and event label columns from the feature data
    time_label <- tryCatch({
        labelled_feature_data[, surv_time_label]
    }, error = function(e) {
        stop(paste("trainCoxModel:Time column not found in provided feature set:", surv_time_label))
    })
    event_label <- tryCatch ({
        labelled_feature_data[, surv_event_label] 
    }, error = function(e) {
    stop(paste("trainCoxModel:Event column not found in provided feature set:", surv_event_label))
    })

    # Fit a CoxPH model to the training features
    model_fit <- coxph(Surv(time_label, event_label) ~ .,
                       x = TRUE,
                       y = TRUE,
                       method = "breslow",
                       data = train_feature_data)

    # Get weights from model and return them
    return(model_fit$coefficients)
}


#' Function to test a CPH model with weights on a set of features.
#' 
#' @param test_labelled_features_file_path Path to the file containing the test features with outcome labels included.
#' @param surv_time_label Label for the time column in the test features file.
#' @param surv_event_label Label for the event column in the test features file.
#' @param model_feature_list List of feature names to use for the Cox model.
#' @param model_feature_weights Vector of weights for the Cox model
#' 
#' @return vector of test results.
testCoxModel <- function(labelled_feature_data,
                         surv_time_label,
                         surv_event_label,
                         model_feature_list,
                         model_feature_weights){ #nolint
    
    # Get only features selected for the model
    test_feature_data <- tryCatch({
        labelled_feature_data[, model_feature_list]
    }, error = function(e) {
        stop(paste("testCoxModel:Model features not found in provided feature set:", model_feature_list, "\n"))
    })

    # Convert the features dataframe to a matrix
    test_feature_matrix <- data.matrix(test_feature_data)

    # Multiply the feature matrix by the weights - this is applying the Cox model to the test set
    feature_hazards <- test_feature_matrix %*% model_feature_weights

    # Get the time and event label columns from the feature data
    time_label <- tryCatch({
        labelled_feature_data[, surv_time_label]
    }, error = function(e) {
        stop(paste("testCoxModel:Time column not found in provided feature set:", surv_time_label))
    })
    event_label <- tryCatch ({
        labelled_feature_data[, surv_event_label] 
    }, error = function(e) {
    stop(paste("testCoxModel:Event column not found in provided feature set:", surv_event_label))
    })

    # Calculate concordance index for the test set
    performance_results <- concordance.index(x = feature_hazards,
                                             surv.time = time_label,
                                             surv.event = event_label,
                                             method = "noether",
                                             alpha = 0.5,
                                             alternative = "two.sided")

    return(performance_results)
}

#' Function to run MRMR feature selection and train a CoxPH model on the selected features.
#' Returns the feature weights and performance results for the best performing model based on the concordance index.
#' 
#' @param labelled_train_data Data.frame containing the training features with outcome labels included.
#' @param n_features Number of features to use for MRMR model training. Default is 30.
#' @param n_solutions Number of solutions to use for MRMR model training. Default is 1. Must be greater than 1 for the bootstrap method.
#' @param mrmr_method Name of the method to use for MRMR feature selection. Must be either "classic" or "bootstrap". The default is "classic".
#' @param surv_time_label Label for the time column in labelled_train_data.
#' @param surv_event_label Label for the event column in labelled_train_data.
#' 
#' @return A named list containing the model weights, concordance index, confidence intervals, and p-value for the best performing model.
trainMRMRCoxModel <- function(labelled_train_data,
                              n_features = 30,
                              n_solutions = 1,
                              mrmr_method='classic',
                              surv_time_label="survival_time_in_years",
                              surv_event_label="survival_event_binary") {

    # Get the feature data for the training set with outcome labels removed
    train_feature_data <- dropLabelsFromFeatureData(labelled_train_data, labels_to_drop=c("patient_ID", surv_time_label, surv_event_label))
    
    # Perform mRMR feature selection using the specified method
    if (mrmr_method == 'classic' | n_solutions == 1) {
        all_mrmr_solutions_matrix <- runMRMRClassic(train_feature_data, n_features)
    } else if (mrmr_method == 'bootstrap') {
        # Will return a n_features x n_solutions matrix with feature indices into train_data as values
        all_mrmr_solutions_matrix <- runMRMRBootstrap(train_feature_data,
                                        n_features,
                                        n_solutions)
    } else { # Have this here so future methods can be added (e.g. exhaustive mRMRe)
        print('Error: Invalid mRMR method. Must be either "classic" or "bootstrap".')
        return(NULL)
    }

    # Initialize best solution CPH model report with empty values
    best_solution_model_data = makeCPHModelReport()

    # Loop through solutions to fit CPH models and determine best one
    for (solution_idx in 1:n_solutions) {
        # Get the solution, will be a vector of length n_features
        solution_feature_indices <- all_mrmr_solutions_matrix[,solution_idx]

        # # Get the feature names for the solution to pass to the Cox model
        solution_feature_names <- names(train_feature_data)[solution_feature_indices]

        # Train the Cox model with the solution features
        # Returns a named list of model weights
        solution_model_feature_weights <- trainCoxModel(labelled_train_data,
                                                    surv_time_label = surv_time_label,
                                                    surv_event_label = surv_event_label,
                                                    model_feature_list = solution_feature_names)
    
        # Get the model performance with the solution weights
        solution_performance_results <- testCoxModel(labelled_train_data,
                                                    surv_time_label = surv_time_label,
                                                    surv_event_label = surv_event_label,
                                                    model_feature_list = solution_feature_names,
                                                    model_feature_weights = solution_model_feature_weights)
        # Get the concordance index for the solution
        solution_c_index <- solution_performance_results$c.index

        # Compare to best solution so far
        if (solution_c_index > best_solution_model_data$ci) {
            # Update best solution
            best_solution_model_data <- makeCPHModelReport(solution_model_feature_weights, solution_performance_results)
        }
    }

    return (best_solution_model_data)
}


trainKFoldMRMRCoxModel <- function(labelled_feature_data,
                                   k = 5,
                                   n_features = 30,
                                   surv_time_label="survival_time_in_years",
                                   surv_event_label="survival_event_binary" ) { #nolint
    # Initialize list to store the feature weights and validation performance results for each fold
    fold_results <- list()
    # Initialize holder for best fold index
    best_fold_idx <- 0
    # Initialize empty best fold model CPH report
    best_fold_model_data <- makeCPHModelReport()

    # Generate k-fold cross-validation folds with the patient IDs 
    folds <- caret::createFolds(labelled_feature_data$patient_ID, k = k)

    for (i in 1:k) {
        print(paste("Training fold:", i))

        # Initialize empty CPH model report for the current fold
        fold_results <- makeCPHModelReport()

        # Get the training data for the current fold
        train_fold <- labelled_feature_data[-folds[[i]],]
        # Get the validation data for the current fold
        val_fold <- labelled_feature_data[folds[[i]],]

        # Run classic MRMR feature selection and train a CoxPH model on the selected features
        trained_model_results <- trainMRMRCoxModel(train_fold,
                                            n_features = n_features,
                                            n_solutions = 1,
                                            mrmr_method = 'classic',
                                            surv_time_label = surv_time_label,
                                            surv_event_label = surv_event_label)
        # Get the trained weights for the model
        train_selected_weights <- trained_model_results$features

        # Run trained CPH model on the validation set
        validation_performance_results <- testCoxModel(val_fold,
                                surv_time_label = surv_time_label,
                                surv_event_label = surv_event_label,
                                model_feature_list = names(train_selected_weights),
                                model_feature_weights = train_selected_weights)

        # Organize the validation results into a named list
        validation_model_data <- makeCPHModelReport(train_selected_weights, validation_performance_results)
        # Store results for the current fold
        list[[i]] <- validation_model_data

        # Check if validation results are better than best fold results
        if (validation_model_data$ci > best_fold_model_data$ci) {
            # Update best fold model data
            best_fold_model_data <- validation_model_data
            # Update best fold index
            best_fold_idx <- i
        }
    } # end for loop

    k_fold_results <- list(best_fold_idx = best_fold_idx, best_fold_model_data = best_fold_model_data, fold_results = fold_results)

    return(k_fold_results)
}

