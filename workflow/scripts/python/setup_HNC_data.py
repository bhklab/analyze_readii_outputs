from readii_analysis.run_data_setup_for_prediction_models import run_data_setup_for_prediction_models

# RADCURE
run_data_setup_for_prediction_models(DATASET_NAME="RADCURE", 
                                     EXTRACTION_METHOD="radiomic", 
                                     RAW_FEATURE_DIR_NAME="readii_outputs")
run_data_setup_for_prediction_models(DATASET_NAME="RADCURE", 
                                     EXTRACTION_METHOD="deep_learning", 
                                     RAW_FEATURE_DIR_NAME="fmcib_outputs")

print("\n_______________________________________________________________________________________________________________\n")

# Head-Neck-Radiomics-HN1
run_data_setup_for_prediction_models(DATASET_NAME="Head-Neck-Radiomics-HN1", 
                                     EXTRACTION_METHOD="radiomic", 
                                     RAW_FEATURE_DIR_NAME="readii_outputs")
run_data_setup_for_prediction_models(DATASET_NAME="Head-Neck-Radiomics-HN1", 
                                     EXTRACTION_METHOD="deep_learning",
                                     RAW_FEATURE_DIR_NAME="fmcib_outputs")

print("\n_______________________________________________________________________________________________________________\n")

# HNSCC
run_data_setup_for_prediction_models(DATASET_NAME="HNSCC", 
                                     EXTRACTION_METHOD="radiomic", 
                                     RAW_FEATURE_DIR_NAME="readii_outputs")
run_data_setup_for_prediction_models(DATASET_NAME="HNSCC", 
                                     EXTRACTION_METHOD="deep_learning", 
                                     RAW_FEATURE_DIR_NAME="fmcib_outputs")

