source("workflow/scripts/R/cox_modelling.r")
source("workflow/scripts/R/io.r")


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
    tryCatch({
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
        yaml::write_yaml(final_signature, output_file)
        close(output_file)
        }, error = function(e) {
            print(paste("Signature file not found for signature:", signature_name))
            
            # Setup the final signature dictionary
            final_signature <- list(signature = model_feature_weights)

            # Setup output file name
            output_file <- file(paste(output_dir, signature_name, ".yaml", sep = ""), "w")
            yaml::write_yaml(final_signature, output_file)
            close(output_file)
    })
}


#' Function to read in a CPH signature file and get the feature names and weights
#' 
#' @param signature_name Name of the signature to read in, should have a signature.yaml file in the signatures folder. Weights are optional in the file.
#' 
#' @return list of feature names and weights
loadSignatureYAML <- function(signature_name) { #nolint 
    # Paste together the path to the signature file
    signature_file_path <- paste("workflow/signatures/", signature_name, ".yaml", sep = "")
    # Load the signature file
    signature_config <- loadYAMLFile(signature_file_path)
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
#' @param k Number of k-folds to use for MRMR model training. Default is 1.
#' @param outputDir Directory to save the trained weights signature file to. Default is "workflow/signatures/".
#' @param overwriteSignature Boolean to indicate whether to overwrite existing signature weights. Default is FALSE.
#' @param testSignature Boolean to indicate whether to run test data on the signature. The dataset provided must have a training/test split. Default is FALSE.
#' 
#' @return vector of trained weights.
createSignature <- function(dataset_config_file_path, signature_name, feature_type = "radiomics", k=1, output_dir = "workflow/signatures/", overwrite_signature = FALSE, test_signature = FALSE) { #nolint
    dataset_config <- loadYAMLFile(dataset_config_file_path)
    # Name of the dataset to run CPH on
    dataset_name <- dataset_config$dataset_name

    # Check if dataset has a training/test split needed for signature creation
    if (dataset_config$train_test_split$split == FALSE) {
        print("Dataset must have a training subset to create a signature.")
        stop()
    }

    tryCatch({
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
    }, error = function(e) {
        print(paste("Signature file not found for signature:", signature_name))
        print("Creating new signature.")
    })

    # Path to training radiomics features from the original image (not a negative control)
    feature_dir_path <- paste("procdata", dataset_name, feature_type, "train_test_split", "train_features", sep="/")
    train_feature_file_path <- list.files(feature_dir_path, pattern="train_labelled.*original.*csv", ignore.case=TRUE, full.names=TRUE) # re-written to use regex to support both radiomic and deep learning (fmcib) features
    print(paste("Original feature path:", train_feature_file_path))
    # Load the feature data as a dataframe
    labelled_feature_data <- loadDataFile(train_feature_file_path)

    if (!exists("sig_feature_names")) {
        if (feature_type == "radiomics") {
            sig_feature_names <- names(labelled_feature_data)
        } else if (feature_type == "deep_learning") {
            sig_feature_names <- names(labelled_feature_data)[grep("^pred", names(labelled_feature_data))]
        }
    }

    # Fit a CoxPH model to the training radiomics features
    if (k <= 1){    
        trained_weights <- c(signature_name = trainCoxModel(labelled_feature_data,
                                                            surv_time_label = "survival_time_in_years",
                                                            surv_event_label = "survival_event_binary",
                                                            model_feature_list = sig_feature_names))
    }
    else{
        print(paste("Training with", k, "folds"))
        trained_weights <- trainKFoldCoxModel(labelled_feature_data,
                                              surv_time_label = "survival_time_in_years",
                                              surv_event_label = "survival_event_binary",
                                              model_feature_list = sig_feature_names,
                                              k)
    }
    
    print(paste("Done training! Saving signatures...", names(trained_weights)))
    for (signature_key in names(trained_weights)) {
        # Save out the model weights for the signature as a yaml file
        saveSignatureYAML(signature_name = signature_key,
                          model_feature_weights = trained_weights[signature_key],
                          output_dir = output_dir,
                          overwrite_signature = overwrite_signature)
    }

    print("Done saving!")
    # Apply the signature to the test set if specified
    if (test_signature == TRUE) {
        applySignature(dataset_config_file_path = dataset_config_file_path, feature_type = feature_type,
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
applySignature <- function(dataset_config_file_path, signature_name, feature_type = "radiomics") { #nolint

    # Load in the config file for the dataset
    dataset_config <- loadYAMLFile(dataset_config_file_path)

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
        feature_dir_path <- paste0("procdata/", dataset_name, "/", feature_type, "/train_test_split/test_features", sep="")
    } else {
        # Path to directory containing validation radiomics features
        feature_dir_path <- paste0("procdata/", dataset_name, "/", feature_type, "/features", sep="")
    }

    # Get relative paths to all radiomic feature files containing the outcome labels and save in character vector
    # test set will be test_labelled_radiomicfeatures_* and validation set will be train_labelled_radiomicfeatures_*
    feature_files <- list.files(feature_dir_path, pattern="labelled.*\\.csv$", ignore.case=TRUE, full.names=TRUE)

    # Get all the feature files in the directory
    feature_datas <= lapply(feature_files, loadDataFile)

    # Run the signature CPH model on each feature set
    cph_model_results <- lapply(feature_datas, testCoxModel, 
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