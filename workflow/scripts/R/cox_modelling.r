library(survival)
library(survcomp)
# library(readxl)
library(tools)
# library(yaml)

##### VARIABLES #####
# Main directory - directory containing labelled feature csvs
# Negative controls - list of negative controls to run CPH on
# Model features - list of features to use for CPH model
# Model weights - matrix of weights to apply for trained CPH model


###### STEPS #######
# Pretrained CPH Model
# 1. Load csv file for testing patient set
# 2. Drop the patient ID and survival labels from the dataframe
# 3. Confirm that features in model features list match the radiomics
# 3. Convert the features dataframe to a data.matrix
# 4. Multiple the weight matrix by the features
# 5. Set up time and event label
# 6. Calculate concordance index 


# From scratch CPH Model training
# 1. Load csv file of training patient set 
# 2. Drop the patient ID and survival labels from the dataframe
# 3. Confirm that features in model features list match the radiomics
# 4. Fit cph model with the select features with the training set
# 5. Get weights for model
# 6. Return the trained CPH model

#' Function to load in the feature data file for CPH model training or testing
#'
#' @param data_file_path A string path to the file to load.
#' 
#' @return A data.table containing the loaded data.
loadDataFile <- function(data_file_path) { #nolint
    checkmate::assertFile(data_file_path, access = "r", extension = c("csv", "xlsx"))

    switch(tools::file_ext(data_file_path),
        "csv" = read.csv(data_file_path, header = TRUE, sep = ",", check.names = FALSE),
        "xlsx" = readxl::read_excel(data_file_path)
    ) 
}


#' Function to train a CoxPH model to make a signature based on a set of features.
#' Will return the trained weights for the model.
#' 
#' @param train_labelled_features_file_path Path to the file containing the training features with outcome labels included.
#' @param surv_time_label Label for the time column in the training features file.
#' @param surv_event_label Label for the event column in the training features file.
#' @param model_feature_list List of feature names to use for the Cox model.
#' 
#' @return vector of trained weights.
trainCoxModel <- function(train_labelled_features_file_path,
                          surv_time_label,
                          surv_event_label,
                          model_feature_list){ #nolint

    # Load the feature data as a dataframe
    labelled_feature_data <- loadDataFile(train_labelled_features_file_path)

     # Get just features selected for the model
    train_feature_data <- tryCatch({
        labelled_feature_data[, model_feature_list]
    }, error = function(e) {
        stop(paste("Model features not found in provided feature set:", model_feature_list))
    })

    # Get the time and event label columns from the feature data
    time_label <- tryCatch({
        labelled_feature_data[, surv_time_label]
    }, error = function(e) {
        stop(paste("Time column not found in provided feature set:", surv_time_label))
    })
    event_label <- tryCatch ({
        labelled_feature_data[, surv_event_label] 
    }, error = function(e) {
    stop(paste("Event column not found in provided feature set:", surv_event_label))
    })

    # Fit a CoxPH model to the training features
    model_fit <- coxph(Surv(time_label, event_label) ~ .,
                       x = TRUE,
                       y = TRUE,
                       method = "breslow",
                       data = train_features_data)

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
testCoxModel <- function(test_labelled_features_file_path,
                         surv_time_label,
                         surv_event_label,
                         model_feature_list,
                         model_feature_weights){ #nolint

    # Load the feature data as a dataframe
    labelled_feature_data <- loadDataFile(test_features_file_path)
    
    # Get only features selected for the model
    test_feature_data <- tryCatch({
        labelled_feature_data[, model_feature_list]
    }, error = function(e) {
        stop(paste("Model features not found in provided feature set:", model_feature_list))
    })

    # Convert the features dataframe to a matrix
    test_feature_matrix <- data.matrix(test_feature_data)

    # Multiply the feature matrix by the weights - this is applying the Cox model to the test set
    feature_hazards <- test_feature_matrix %*% model_feature_weights

    # Get the time and event label columns from the feature data
    time_label <- tryCatch({
        labelled_feature_data[, surv_time_label]
    }, error = function(e) {
        stop(paste("Time column not found in provided feature set:", surv_time_label))
    })
    event_label <- tryCatch ({
        labelled_feature_data[, surv_event_label] 
    }, error = function(e) {
    stop(paste("Event column not found in provided feature set:", surv_event_label))
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


#' Function to save out a CPH signature file with the trained weights. A signature file with the names of the features must already exist.
#' 
#' @param signature_name Name of the signature to load in signature names from yaml file. Will be saved out with the same name.
#' @param model_feature_weights Vector of trained weights for the signature. Must have the same number of values as the number of features in the signature.
#' @param output_dir Directory to save the signature file to. Default is "workflow/signatures/". Must end in a "/".
#' @param overwrite_signature Boolean to indicate whether to overwrite existing signature weights. Default is FALSE.
#' 
#' @return None
saveSignatureYAML <- function(signature_name, model_feature_weights, output_dir = "workflow/signatures/", overwrite_signature = FALSE){ #nolint 
    # Load in the signature file to get the feature names
    signature <- loadSignatureYAML(signature_name)

    # Check if signature weights already exist
    if (signature$weights[1] != 0) {
        # Check whether to overwrite the signature weights
        if (overwrite_signature == TRUE) {
            print("Signature weights are being overwritten.")
        } else {
            print("Signature weights already exist. Set overwriteSignature to TRUE to overwrite.")
            stop()
        }
    }

    # Make sure output directory ends in a "/"
    if (endsWith(output_dir, "/") == FALSE) { output_dir <- paste(output_dir, "/", sep = "") }

    # Convert model weights to list and set up names from existing signature file
    model_feature_list <- as.list(model_feature_weights)
    names(model_feature_list) <- signature$names

    # Name the signature so ouptut is correctly formatted
    final_signature <- list(signature = model_feature_list)

    # Setupt output file name
    output_file <- file(paste(output_dir, signature_name, ".yaml", sep = ""), "w")
    # Write out the signature
    write_yaml(final_signature, output_file)
    close(output_file)
}


#' Function to read in a CPH signature file and get the feature names and weights
#' 
#' @param signature_name Name of the signature to read in, should have a signature.yaml file in the signatures folder. Weights are optional in the file.
#' 
#' @return list of feature names and weights
loadSignatureYAML <- function(signature_name) { #nolint 
    # Paste together the path to the signature file
    signature_file_path <- paste("workflow/signatures/", signature_name, ".yaml", sep = "")
    # Confirm it is a yaml file
    checkmate::assert_file(signature_file_path, access = "r", extension = "yaml")
    # Load the signature file
    signature_config <- read_yaml(signature_file_path)
    # Names of the features in the signature
    sig_feature_names <- names(signature_config$signature)
    # Weights for the features in the signature
    sig_weights <- matrix(unlist(signature_config$signature))

    signature <- list(sig_feature_names, sig_weights)
    names(signature) <- c("names", "weights")

    return(signature)
}


#' Function to train a CoxPH model to make a signature based on a set of radiomic features.
#' Will return the trained weights for the model.
#' 
#' @param datasetConfigFilePath Path to the config file for the dataset
#' @param signatureName Name of the signature to train, should have a signature.yaml file in the signatures folder. If this has any weights, they will be overwritten.
#' @param outputDir Directory to save the trained weights signature file to. Default is "workflow/signatures/".
#' @param overwriteSignature Boolean to indicate whether to overwrite existing signature weights. Default is FALSE.
#' @param testSignature Boolean to indicate whether to run test data on the signature. The dataset provided must have a training/test split. Default is FALSE.
#' 
#' @return vector of trained weights.
createSignature <- function(dataset_config_file_path, signature_name, output_dir = "workflow/signatures/", overwrite_signature = FALSE, test_signature = FALSE) { #nolint
    dataset_config <- read_yaml(dataset_config_file_path)
    # Name of the dataset to run CPH on
    dataset_name <- dataset_config$dataset_name

    # Check if dataset has a training/test split needed for signature creation
    if (datasetConfig$train_test_split$split == FALSE) {
        print("Dataset must have a training subset to create a signature.")
        stop()
    }

    # Signature setup - get the signature features and weights
    signature <- loadSignatureYAML(signature_name)
    sig_feature_names <- signature$names
    sig_weights <- signature$weights

    if (sig_weights[1] != 0) {
        # Check whether to overwrite the signature weights
        if (overwrite_signature == TRUE) {
            print("Signature weights are being overwritten.")
        } else {
            print("Signature weights already exist. Set overwriteSignature to TRUE to overwrite.")
            stop()
        }
    }

    # Path to training radiomics features from the original image (not a negative control)
    train_feature_file_path <- paste("procdata/", dataset_name, "/radiomics/train_test_split/train_features/train_labelled_radiomicfeatures_only_original_", dataset_name, ".csv", sep = "")

    # Fit a CoxPH model to the training radiomics features
    trained_weights <- trainCoxModel(train_feature_file_path,
                                     survTimeLabel = "survival_time_in_years",
                                     survEventLabel = "survival_event_binary",
                                     modelFeatureList = sig_feature_names)

    # Save out the model weights for the signature as a yaml file
    saveSignatureYAML(signature_name = signature_name,
                      model_feature_weights = trained_weights,
                      output_dir = output_dir,
                      overwrite_signature = overwrite_signature)

    # Apply the signature to the test set if specified
    if (test_signature == TRUE) {
        applySignature(dataset_config_file_path = dataset_config_file_path,
                       signature_name = signature_name)
    }

    return(trained_weights)

}


#' Function to apply a trained CPH model to a test or validation set of radiomics features
#' 
#' @param dataset_config_file_path Path to the config file for the dataset
#' @param signature_name Name of the signature to apply, should have a signature.yaml file in the signatures folder
#' 
#' @return None
applySignature <- function(dataset_config_file_path, signature_name) { #nolint
    
    # Load in config file for the dataset
    checkmate::assert_file(dataset_config_file_path, access = "r", extension = "yaml")
    dataset_config <- read_yaml(dataset_config_file_path)

    # Name of the dataset to run CPH on
    dataset_name <- dataset_config$dataset_name

    # Signature setup - get the signature features and weights
    signature <- loadSignatureYAML(signature_name)
    sig_feature_names <- signature$names
    sig_weights <- signature$weights

    print(paste("DATASET: ", dataset_name))
    print(paste("SIGNATURE: ", signature_name))

    # Check if applying signature to test or validation set based on training/test split
    if (dataset_config$train_test_split$split == TRUE) {
        # Path to directory containing test radiomics features
        feature_dir_path <- paste0("procdata/", dataset_name, "/radiomics/train_test_split/test_features", sep="")
    } else {
        # Path to directory containing validation radiomics features
        feature_dir_path <- paste0("procdata/", dataset_name, "/radiomics/features", sep="")
    }

    # Get relative paths to all radiomic feature files containing the outcome labels and save in character vector
    # test set will be test_labelled_radiomicfeatures_* and validation set will be train_labelled_radiomicfeatures_*
    feature_files <- list.files(feature_dir_path, pattern="labelled.*\\.csv$", ignore.case=TRUE, full.names=TRUE)

    # Run the signature CPH model on each feature set
    cph_model_results <- lapply(feature_files, testCoxModel, 
                                surv_time_label = "survival_time_in_years",
                                surv_event_label = "survival_event_binary",
                                model_feature_list = sig_feature_names,
                                model_feature_weights = sig_weights)

    #TODO: update this to just be the file name, not the full path
    # Make the file path of the feature file the name of the list element
    names(cph_model_results) <- feature_files

    # Return the list of results
    return(cph_model_results)
}



##### VARIABLES #####
# Main directory - directory containing labelled feature csvs
# Negative controls - list of negative controls to run CPH on
# Model features - list of features to use for CPH model
# Model weights - matrix of weights to apply for trained CPH model
# dataset_config_path <- "workflow/config/Head-Neck-Radiomics-HN1.yaml"
# signature_name <- "aerts_original"
# output_dir <- "workflow/signatures"

# # trained_weights <- create_signature(dataset_config_path, signature_name, output_dir = output_dir, test = TRUE)
# test_results = applySignature(dataset_config_path, signature_name)
