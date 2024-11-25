
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

install.packages("BiocManager", repos = "http://cran.us.r-project.org")
install.packages("mRMRe", repos = "http://cran.us.r-project.org")
install.packages("checkmate", repos = "http://cran.us.r-project.org")
BiocManager::install("survcomp")


runMRMRBootstrap <- function(data, n_features, n_solutions) {
    # setup data for the mRMRe function
    mrmr_data <- mRMR.data(data = data)

    # run mRMR with bootstrap ensemble method
    # will return n_solutions lists of n_features features
    mrmr_solutions <- mRMR.ensemble(data = mrmr_data,
                                    solution_count = n_solutions,
                                    feature_count = n_features,)
}
