

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


#' Function to load in a YAML file with proper checks
#'
#' @param yaml_file_path A string path to the file to load.
#' 
#' @return A data.table containing the loaded data.
loadYAMLFile <- function(yaml_file_path) { #nolint
    checkmate::assertFile(yaml_file_path, access = "r", extension = "yaml")
    yaml::read_yaml(yaml_file_path)
}