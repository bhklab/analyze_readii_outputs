
import pandas as pd
import os


def select_GTVp_for_RADCURE_readii_outputs(RADCURE_readii_output_dir:str):
    """ Function to select out the GTVp radiomic features from the RADCURE readii outputs.
    """
    for feature_file in os.listdir(RADCURE_readii_output_dir):
        if feature_file.endswith(".csv"):
            print(f"Processing {feature_file}")
            feature_df = pd.read_csv(os.path.join(RADCURE_readii_output_dir, feature_file))

            # Select out the rows with the GTVp label in roi column
            GTVp_df = feature_df[feature_df["roi"] == "GTVp"]

            # Save out the GTVp dataframe
            GTVp_df.to_csv(os.path.join(RADCURE_readii_output_dir, feature_file), index=False)


if __name__ == "__main__":
    RADCURE_readii_output_dir = "../../../rawdata/RADCURE/readii_outputs/"

    select_GTVp_for_RADCURE_readii_outputs(RADCURE_readii_output_dir)
