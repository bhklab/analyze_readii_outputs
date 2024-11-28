# analyze_readii_outputs
Code for analyzing outputs from the READII package or the readii-orcestra pipeline.


# Workflow Contents
* shell script for setting up the directory structure
* Python Jupyter Notebook for pre-processing the clinical and image features data and setting up the pre-existing radiomic signatures
* R notebook for performing feature selection and CPH modeling

# Data Organization

Raw data (currently what _would_ be the deconstructed output from ORCESTRA object)
```raw
└──  rawdata
    └──  {DATASETNAME}
        ├──  clinical
        ├──  fmcib_outputs
        └──  readii_outputs
        {DATASETNAME}_READII-RADIOMICS_MAE.RDS
```

Processed data = filtered clinical, radiomic, and deep learning features, possibly split into training and test sets
```raw
└──  procdata
    └──  {DATASETNAME}
        └──  clinical
             ├──  cleaned_filtered_clinical_{DATASETNAME}.csv
             └──  [OPTIONAL] train_test_labelled_clinical_{DATASETNAME}.csv
        └──  radiomics
             ├──  clinical
                  └──  merged_clinical_{DATASETNAME}.csv
             ├──  features
                  ├──  merged_radiomicfeatures_{image_type}_{DATASETNAME}.csv
                  └──  labelled_radiomicfeatures_only_{image_type}_{DATASETNAME}.csv
             └──  [OPTIONAL] train_test_split
                  ├──  clinical
                       └──  train_merged_clinical_{DATASETNAME}.csv
                       └──  test_merged_clinical_{DATASETNAME}.csv
                  ├──  train_features
                       └──  train_labelled_radiomicfeatures_only_{image_type}_{DATASETNAME}.csv
                  └──  test_features
                       └──  test_labelled_radiomicfeatures_only_{image_type}_{DATASETNAME}.csv
        └──  deep_learning
             ├──  clinical
                  └──  merged_clinical_{DATASETNAME}.csv
             ├──  features
                  └──  merged_fmcibfeatures_{image_type}_{DATASETNAME}.csv
                  └──  labelled_fmcibfeatures_only_{image_type}_{DATASETNAME}.csv
             └──  [OPTIONAL] train_test_split
                  ├──  clinical
                       └──  train_merged_clinical_{DATASETNAME}.csv
                       └──  test_merged_clinical_{DATASETNAME}.csv
                  ├──  train_features
                       └──  train_labelled_fmcibfeatures_only_{image_type}_{DATASETNAME}.csv
                  └──  test_features
                       └──  test_labelled_fmcibfeatures_only_{image_type}_{DATASETNAME}.csv           
```
# TODO:

- [ ] Implement logger
- [ ] Implement ORCESTRA download
- [ ] Implement MAE deconstructor
    - [ ] Unpack clinical data --> save to csv
    - [ ] Unpack radiomic features --> save each experiment to csv
    - [ ] Unpack deep learning features --> save each experiment to csv
    - [ ] Get list of experiments, specifically the negative controls
- [ ] Make this into a snakemake pipeline
- [ ] Implement config file creation if one is not present
- [ ] Move data_setup_for_modelling from scripts into notebooks
- [ ] Finish implementing survival time and event setup as functions
- [x] Supports MRMR training over k folds
- [ ] Doesn't support MRMR training over 1 fold
- [ ] Doesn't support regular training over k folds
- [ ] Doesn't support loading model weights across k folds