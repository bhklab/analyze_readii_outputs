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

