
# install.packages("BiocManager", repos = "http://cran.us.r-project.org")
# install.packages("mRMRe", repos = "http://cran.us.r-project.org")
# install.packages("checkmate", repos = "http://cran.us.r-project.org")
# BiocManager::install("survcomp")

source("workflow/scripts/R/io.r")
source("workflow/scripts/R/data_processing.r")
source("workflow/scripts/R/signature_handling.r")
source("workflow/scripts/R/mrmr_functions.r")

library(survival)
library(survcomp)
library(tools)
# library(caret)
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
#                 best_k <- k_mrmr
#                 best_feats <- mrmr_feature_list
#                 best_coefs <- coefs
#             }
#         }
#     return(setNames(as.list(best_coefs), best_feats))
# }

# trainKFoldCoxModel <- function(labelled_feature_data,
#                                surv_time_label,
#                                surv_event_label,
#                                model_feature_list,
#                                k){ #nolint
#     trained_weights <- list()
#     folds <- createFolds(labelled_feature_data$patient_ID, k = k)
#     for (i in 1:k) {
#         print(paste("Training fold:", i))
#         train_fold <- labelled_feature_data[folds[[i]],]
#         val_fold <- labelled_feature_data[-folds[[i]],]
#         trained_model <- trainMRMRCoxModel(train_fold,
#                                                   val_fold,
#                                                   surv_time_label = "survival_time_in_years",
#                                                   surv_event_label = "survival_event_binary",
#                                                   model_feature_list = model_feature_list)
#         trained_weights[[i]] <- trained_model
#         }
#     return(trained_weights)
# }

#' Function to train a CoxPH model to make a signature based on a set of features.
#' Will return the trained weights for the model.
#' 
#' @param train_labelled_features_file_path Path to the file containing the training features with outcome labels included.
#' @param surv_time_label Label for the time column in the training features file.
#' @param surv_event_label Label for the event column in the training features file.
#' @param model_feature_list List of feature names to use for the Cox model.
#' 
#' @return vector of trained weights.
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
        stop(paste("testCoxModel:Model features not found in provided feature set:", model_feature_list))
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


trainMRMRCoxModel <- function(labelled_train_data,
                              n_features,
                              n_solutions,
                              mrmr_method='bootstrap',
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
        print('Error: Invalid mRMR method')
        return(NULL)
    }

    # Initialize best solution holder variables
    best_c_index = 0
    best_features_and_weights = list()

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

