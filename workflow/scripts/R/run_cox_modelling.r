source("workflow/scripts/R/cox_modelling.r")

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


##### VARIABLES #####
# Main directory - directory containing labelled feature csvs
# Negative controls - list of negative controls to run CPH on
# Model features - list of features to use for CPH model
# Model weights - matrix of weights to apply for trained CPH model

#### SETUP ####
signature_name <- "aerts_RADCURE"
output_dir <- "workflow/signatures"

#### TRAINING ####
train_dataset_name <- "RADCURE"
train_dataset_config_path <- paste0("workflow/config/", train_dataset_name, ".yaml", sep="")

trained_weights <- createSignature(train_dataset_config_path, signature_name, output_dir = output_dir, overwrite_signature = TRUE, test_signature = TRUE)

#### TESTING ####
# test_dataset_name <- "HNSCC"
# test_dataset_config_path <- paste0("workflow/config/", test_dataset_name, ".yaml", sep="")
# test_results <- applySignature(test_dataset_config_path, signature_name)

# TODO: add print and saving out results options