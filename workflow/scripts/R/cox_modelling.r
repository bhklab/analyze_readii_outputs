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

loadDataFile <- function(dataFilePath) { #nolint
    dataFileType = file_ext(dataFilePath)
    if (length(dataFileType) == 0) {
        # No file passed, return 0
        return(NULL)
    }
    if (dataFileType == "csv") {
        loadedDataframe <- read.csv(dataFilePath, header = TRUE, sep = ",", check.names = FALSE)
    } else if (dataFileType == "xlsx") {
        loadedDataframe <- read_excel(dataFilePath)
    } else {
        stop("Radiogenomic data file must be a .csv or .xlsx.")
    }

    return(loadedDataframe)
}


trainCoxModel <- function(csvTrainingFeatures,
                          survTimeLabel,
                          survEventLabel,
                          modelFeatureList){ #nolint

    trainRadiomicsData <- loadDataFile(csvTrainingFeatures)

    timeData <- trainRadiomicsData[, survTimeLabel]
    eventData <- trainRadiomicsData[, survEventLabel]

    trainModelFeatures <- as.data.frame(trainRadiomicsData[, modelFeatureList])

    model.fit <- coxph(Surv(timeData, eventData) ~ .,
                       x = TRUE,
                       y = TRUE,
                       method = "breslow",
                       data = trainModelFeatures)

    # Get weights from model and return those two
    return(model.fit$coefficients)
}

testCoxModel <- function(csvTestingFeatures,
                         survTimeLabel,
                         survEventLabel,
                         modelFeatureList,
                         modelFeatureWeights){ #nolint

    testRadiomicsData <- loadDataFile(csvTestingFeatures)

    # INSERT CHECK FOR RADIOMIC COLUMNS AND MODELFEATURELIST

    testRadiomicFeatures <- testRadiomicsData[, modelFeatureList]

    testRadMatrix <- data.matrix(testRadiomicFeatures)

    radiomicHazards <- testRadMatrix %*% modelFeatureWeights

    timeLabel <- testRadiomicsData[, survTimeLabel]
    eventLabel <- testRadiomicsData[, survEventLabel]

    performanceResults <- concordance.index(x = radiomicHazards,
                                            surv.time = timeLabel,
                                            surv.event = eventLabel,
                                            method = "noether",
                                            alpha = 0.5,
                                            alternative = "two.sided")

    return(performanceResults)
}


saveSignature <- function(signatureName, modelFeatureWeights, outputDir){ #nolint 
    # Save out model weights for CPH model

    signature <- signatureYAMLSetup(signatureName)

    # Convert model weights to list and set up names
    modelFeatureList <- as.list(modelFeatureWeights)
    names(modelFeatureList) <- signature$names

    # Name the signature so ouptut is correctly formatted
    finalSignature <- list(signature = modelFeatureList)

    # Setupt output file name
    outputFile <- file(paste(outputDir, "/", signatureName, ".yaml", sep = ""), "w")
    # Write out the signature
    write_yaml(finalSignature, outputFile)
    close(outputFile)
}


signatureYAMLSetup <- function(signatureName) { #nolint 
    signatureConfig <- read_yaml(paste("workflow/signatures/", signatureName, ".yaml", sep = ""))
    # Names of the features in the signature
    sigFeatureNames <- names(signatureConfig$signature)
    # Weights for the features in the signature
    sigWeights <- matrix(unlist(signatureConfig$signature))

    signature <- list(sigFeatureNames, sigWeights)
    names(signature) <- c("names", "weights")

    return(signature)
}


setupOutcomeStatus <- function(datasetConfig){
    timeLabel <- datasetConfig$outcome_variables$time_label
    eventLabel <- datasetConfig$outcome_status$event_label

    return(list("timeLabel" = timeLabel, "eventLabel" = eventLabel))
}



createSignature <- function(configFilePath, signatureName, outputDir, test = FALSE) { #nolint
    datasetConfig <- read_yaml(configFilePath)
    # Name of the dataset to run CPH on
    datasetName <- datasetConfig$dataset_name

    # Signature setup - get the signature features and weights
    signature <- signatureYAMLSetup(signatureName)
    sigFeatureNames <- signature$names
    sigWeights <- signature$weights

    if (sigWeights[1] != 0) {
        print("Signature weights are being overwritten.")
    }

    # Path to directory containing radiomics features
    featureDirPath <- paste("../../../procdata", datasetName, "all_features", sep="/") #, "test/labelled_readii/", sep="/") 
    

    # Determine whether to load data as train test split or not
    if (datasetConfig$train_test_split$split == TRUE) {
        trainDirPath <- paste(featureDirPath, "/train/labelled_readii/training_labelled_radiomicfeatures_original_", datasetName, ".csv", sep = "")
    } else {
        print("Dataset must have a training subset to create a signature.")
        stop()
    }

    # outcomeLabels <- setupOutcomeStatus(datasetConfig)

    trainedWeights <- trainCoxModel(csvTrainingFeatures = trainDirPath,
                                    survTimeLabel = "survival_time_in_years",
                                    survEventLabel = "survival_event_binary",
                                    modelFeatureList = sigFeatureNames)

    saveSignature(signatureName = signatureName,
                  modelFeatureWeights = trainedWeights,
                  outputDir = outputDir)

    if (test == TRUE) {
        applySignature(configFilePath = configFilePath,
                       signatureName = signatureName)
    }

    return(trainedWeights)

}

#' Function to apply a trained CPH model to a test or validation set of radiomics features
#' 
#' @param datasetConfigFilePath Path to the config file for the dataset
#' @param signatureName Name of the signature to apply, should have a signature.yaml file in the signatures folder
#' 
#' @return None
applySignature <- function(datasetConfigFilePath, signatureName) { #nolint
    
    # Load in config file for the dataset
    checkmate::assertFile(datasetConfigFilePath, access = "r", extension = "yaml")
    datasetConfig <- read_yaml(datasetConfigFilePath)

    # Name of the dataset to run CPH on
    datasetName <- datasetConfig$dataset_name

    # Signature setup - get the signature features and weights
    signature <- signatureYAMLSetup(signatureName)
    sigFeatureNames <- signature$names
    sigWeights <- signature$weights

    print(paste("DATASET: ", datasetName))
    print(paste("SIGNATURE: ", signatureName))

    # Check if applying signature to test or validation set based on training/test split
    if (datasetConfig$train_test_split$split == TRUE) {
        # Path to directory containing test radiomics features
        featureDirPath <- paste0("procdata/", datasetName, "/radiomics/train_test_split/test_features", sep="")
    } else {
        # Path to directory containing validation radiomics features
        featureDirPath <- paste0("procdata/", datasetName, "/radiomics/features", sep="")
    }

    # Get relative paths to all radiomic feature files containing the outcome labels and save in character vector
    # test set will be test_labelled_radiomicfeatures_* and validation set will be train_labelled_radiomicfeatures_*
    featureFiles <- list.files(featureDirPath, pattern="labelled.*\\.csv$", ignore.case=TRUE, full.names=TRUE)

    # Run the signature CPH model on each feature set
    cphModelResults <- lapply(featureFiles, testCoxModel, 
                              survTimeLabel = "survival_time_in_years",
                              survEventLabel = "survival_event_binary",
                              modelFeatureList = sigFeatureNames,
                              modelFeatureWeights = sigWeights)

    #TODO: update this to just be the file name, not the full path
    # Make the file path of the feature file the name of the list element
    names(cphModelResults) <- featureFiles

    # Return the list of results
    return(cphModelResults)
}



##### VARIABLES #####
# Main directory - directory containing labelled feature csvs
# Negative controls - list of negative controls to run CPH on
# Model features - list of features to use for CPH model
# Model weights - matrix of weights to apply for trained CPH model
datasetConfigPath <- "workflow/config/Head-Neck-Radiomics-HN1.yaml"
signatureName <- "aerts_original"
outputDir <- "workflow/signatures"

# trainedWeights <- createSignature(datasetConfigPath, signatureName, outputDir = outputDir, test = TRUE)
testResults = applySignature(datasetConfigPath, signatureName)
