library(survival)
library(survcomp)
library(readxl)
library(tools)
library(yaml)

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

load_data_file <- function(data_file_path) { #nolint
    data_file_type = file_ext(data_file_path)
    if (length(data_file_type) == 0) {
        # No file passed, return 0
        return(NULL)
    }
    if (data_file_type == "csv") {
        loaded_dataframe <- read.csv(data_file_path, header = TRUE, sep = ",", check.names = FALSE)
    } else if (data_file_type == "xlsx") {
        loaded_dataframe <- read_excel(data_file_path)
    } else {
        stop("Radiogenomic data file must be a .csv or .xlsx.")
    }

    return(loaded_dataframe)
}


train_cox_model <- function(csv_training_features,
                            surv_time_label,
                            surv_event_label,
                            model_feature_list){ #nolint

    train_radiomics_data <- load_data_file(csv_training_features)

    time_data <- train_radiomics_data[, surv_time_label]
    event_data <- train_radiomics_data[, surv_event_label]

    train_model_features <- as.data.frame(train_radiomics_data[, model_feature_list])

    model_fit <- coxph(Surv(time_data, event_data) ~ .,
                       x = TRUE,
                       y = TRUE,
                       method = "breslow",
                       data = train_model_features)

    # Get weights from model and return those two
    return(model_fit$coefficients)
}

test_cox_model <- function(csv_testing_features,
                           surv_time_label,
                           surv_event_label,
                           model_feature_list,
                           model_feature_weights){ #nolint

    test_radiomics_data <- load_data_file(csv_testing_features)

    # INSERT CHECK FOR RADIOMIC COLUMNS AND MODELFEATURELIST

    test_radiomic_features <- test_radiomics_data[, model_feature_list]

    test_rad_matrix <- data.matrix(test_radiomic_features)

    radiomic_hazards <- test_rad_matrix %*% model_feature_weights

    time_label <- test_radiomics_data[, surv_time_label]
    event_label <- test_radiomics_data[, surv_event_label]

    performance_results <- concordance.index(x = radiomic_hazards,
                                             surv.time = time_label,
                                             surv.event = event_label,
                                             method = "noether",
                                             alpha = 0.5,
                                             alternative = "two.sided")

    return(performance_results)
}


save_signature <- function(signature_name, model_feature_weights, output_dir){ #nolint 
    # Save out model weights for CPH model

    signature <- signature_yaml_setup(signature_name)

    # Convert model weights to list and set up names
    model_feature_list <- as.list(model_feature_weights)
    names(model_feature_list) <- signature$names

    # Name the signature so ouptut is correctly formatted
    final_signature <- list(signature = model_feature_list)

    # Setupt output file name
    output_file <- file(paste(output_dir, "/", signature_name, ".yaml", sep = ""), "w")
    # Write out the signature
    write_yaml(final_signature, output_file)
    close(output_file)
}


signature_yaml_setup <- function(signature_name) { #nolint 
    signature_config <- read_yaml(paste("workflow/signatures/", signature_name, ".yaml", sep = ""))
    # Names of the features in the signature
    sig_feature_names <- names(signature_config$signature)
    # Weights for the features in the signature
    sig_weights <- matrix(unlist(signature_config$signature))

    signature <- list(sig_feature_names, sig_weights)
    names(signature) <- c("names", "weights")

    return(signature)
}


setup_outcome_status <- function(dataset_config){
    time_label <- dataset_config$outcome_variables$time_label
    event_label <- dataset_config$outcome_status$event_label

    return(list("time_label" = time_label, "event_label" = event_label))
}



create_signature <- function(config_file_path, signature_name, output_dir, test = FALSE) { #nolint
    dataset_config <- read_yaml(config_file_path)
    # Name of the dataset to run CPH on
    dataset_name <- dataset_config$dataset_name

    # Signature setup - get the signature features and weights
    signature <- signature_yaml_setup(signature_name)
    sig_feature_names <- signature$names
    sig_weights <- signature$weights

    if (sig_weights[1] != 0) {
        print("Signature weights are being overwritten.")
    }

    # Path to directory containing radiomics features
    feature_dir_path <- paste("../../../procdata", dataset_name, "all_features", sep="/") #, "test/labelled_readii/", sep="/") 
    

    # Determine whether to load data as train test split or not
    if (dataset_config$train_test_split$split == TRUE) {
        train_dir_path <- paste(feature_dir_path, "/train/labelled_readii/training_labelled_radiomicfeatures_original_", dataset_name, ".csv", sep = "")
    } else {
        print("Dataset must have a training subset to create a signature.")
        stop()
    }

    # outcomeLabels <- setupOutcomeStatus(datasetConfig)

    trained_weights <- train_cox_model(csv_training_features = train_dir_path,
                                       surv_time_label = "survival_time_in_years",
                                       surv_event_label = "survival_event_binary",
                                       model_feature_list = sig_feature_names)

    save_signature(signature_name = signature_name,
                   model_feature_weights = trained_weights,
                   output_dir = output_dir)

    if (test == TRUE) {
        apply_signature(config_file_path = config_file_path,
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
apply_signature <- function(dataset_config_file_path, signature_name) { #nolint
    
    # Load in config file for the dataset
    checkmate::assert_file(dataset_config_file_path, access = "r", extension = "yaml")
    dataset_config <- read_yaml(dataset_config_file_path)

    # Name of the dataset to run CPH on
    dataset_name <- dataset_config$dataset_name

    # Signature setup - get the signature features and weights
    signature <- signature_yaml_setup(signature_name)
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
    cph_model_results <- lapply(feature_files, test_cox_model, 
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
dataset_config_path <- "workflow/config/Head-Neck-Radiomics-HN1.yaml"
signature_name <- "aerts_original"
output_dir <- "workflow/signatures"

# trained_weights <- create_signature(dataset_config_path, signature_name, output_dir = output_dir, test = TRUE)
test_results = apply_signature(dataset_config_path, signature_name)
