from readii_analysis.run_correlation_analysis import run_correlation_analysis


run_correlation_analysis(dataset_name="Head-Neck-Radiomics-HN1",
                         extraction_method="radiomic",
                         extracted_feature_dir="features/labelled_features_only")

run_correlation_analysis(dataset_name="Head-Neck-Radiomics-HN1",
                         extraction_method="deep_learning",
                         extracted_feature_dir="features/labelled_features_only")


run_correlation_analysis(dataset_name = "HNSCC",
                         extraction_method = "radiomic",
                         extracted_feature_dir="features/labelled_features_only")

run_correlation_analysis(dataset_name = "HNSCC",
                         extraction_method = "deep_learning",
                         extracted_feature_dir="features/labelled_features_only")


run_correlation_analysis(dataset_name = "RADCURE",
                         extraction_method = "radiomic",
                         extracted_feature_dir="train_test_split/train_features")

run_correlation_analysis(dataset_name = "RADCURE",
                         extraction_method = "deep_learning",
                         extracted_feature_dir="train_test_split/train_features",
                         self_dist_y_upper_bound= 55000)



