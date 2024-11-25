from readii_analysis.run_correlation_analysis import run_correlation_analysis

CORRELATION_METHOD = "pearson"
CORR_COLOR_MAP = "nipy_spectral"

SELF_DIST_NUM_BINS = 450
SELF_DIST_Y_UPPER_BOUND = 55000

CROSS_DIST_NUM_BINS = 450
CROSS_DIST_Y_UPPER_BOUND = 160000


run_correlation_analysis(dataset_name="Head-Neck-Radiomics-HN1",
                         extraction_method="radiomic",
                         extracted_feature_dir="features/labelled_features_only",
                         correlation_method=CORRELATION_METHOD,
                         corr_color_map=CORR_COLOR_MAP,
                         self_dist_num_bins=SELF_DIST_NUM_BINS,
                         self_dist_y_upper_bound=SELF_DIST_Y_UPPER_BOUND,
                         cross_dist_num_bins=CROSS_DIST_NUM_BINS,
                         cross_dist_y_upper_bound=CROSS_DIST_Y_UPPER_BOUND
                         )

print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")

run_correlation_analysis(dataset_name="Head-Neck-Radiomics-HN1",
                         extraction_method="deep_learning",
                         extracted_feature_dir="features/labelled_features_only",
                         correlation_method=CORRELATION_METHOD,
                         corr_color_map=CORR_COLOR_MAP,
                         self_dist_num_bins=SELF_DIST_NUM_BINS,
                         self_dist_y_upper_bound=SELF_DIST_Y_UPPER_BOUND,
                         cross_dist_num_bins=CROSS_DIST_NUM_BINS,
                         cross_dist_y_upper_bound=CROSS_DIST_Y_UPPER_BOUND
                         )

print("\n_______________________________________________________________________________________________________________\n")

run_correlation_analysis(dataset_name = "HNSCC",
                         extraction_method = "radiomic",
                         extracted_feature_dir="features/labelled_features_only",
                         correlation_method=CORRELATION_METHOD,
                         corr_color_map=CORR_COLOR_MAP,
                         self_dist_num_bins=SELF_DIST_NUM_BINS,
                         self_dist_y_upper_bound=SELF_DIST_Y_UPPER_BOUND,
                         cross_dist_num_bins=CROSS_DIST_NUM_BINS,
                         cross_dist_y_upper_bound=CROSS_DIST_Y_UPPER_BOUND
                         )
print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")

run_correlation_analysis(dataset_name = "HNSCC",
                         extraction_method = "deep_learning",
                         extracted_feature_dir="features/labelled_features_only",
                         correlation_method=CORRELATION_METHOD,
                         corr_color_map=CORR_COLOR_MAP,
                         self_dist_num_bins=SELF_DIST_NUM_BINS,
                         self_dist_y_upper_bound=SELF_DIST_Y_UPPER_BOUND,
                         cross_dist_num_bins=CROSS_DIST_NUM_BINS,
                         cross_dist_y_upper_bound=CROSS_DIST_Y_UPPER_BOUND
                         )

print("\n_______________________________________________________________________________________________________________\n")


run_correlation_analysis(dataset_name = "RADCURE",
                         extraction_method = "radiomic",
                         extracted_feature_dir="train_test_split/train_features",
                         correlation_method=CORRELATION_METHOD,
                         corr_color_map=CORR_COLOR_MAP,
                         self_dist_num_bins=SELF_DIST_NUM_BINS,
                         self_dist_y_upper_bound=SELF_DIST_Y_UPPER_BOUND,
                         cross_dist_num_bins=CROSS_DIST_NUM_BINS,
                         cross_dist_y_upper_bound=CROSS_DIST_Y_UPPER_BOUND
                         )

print("\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n")

run_correlation_analysis(dataset_name = "RADCURE",
                         extraction_method = "deep_learning",
                         extracted_feature_dir="train_test_split/train_features",
                         correlation_method=CORRELATION_METHOD,
                         corr_color_map=CORR_COLOR_MAP,
                         self_dist_num_bins=SELF_DIST_NUM_BINS,
                         self_dist_y_upper_bound=SELF_DIST_Y_UPPER_BOUND,
                         cross_dist_num_bins=CROSS_DIST_NUM_BINS,
                         cross_dist_y_upper_bound=CROSS_DIST_Y_UPPER_BOUND
                         )



