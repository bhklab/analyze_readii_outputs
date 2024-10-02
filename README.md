# analyze_readii_outputs
Code for analyzing outputs from the READII package or the readii-orcestra pipeline.


# Workflow Contents
* shell script for setting up the directory structure
* Python Jupyter Notebook for pre-processing the clinical and image features data and setting up the pre-existing radiomic signatures
* R notebook for performing feature selection and CPH modeling


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