import os
import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--fastqs", help="FASTQ directory", required=True)
    parser.add_argument("--rename", help="Mapping file", required=True)
    args = parser.parse_args()
    fastqs_dir = args.fastqs
    rename = args.rename
    if rename is not None:
        df = pd.read_csv(rename, sep="\t", header=None, names=["original_name", "new_name"])
        df = df.dropna()
        # strip path
        df["original_name"] = df["original_name"].str.split("/").str[-1]
        df["new_name"] = df["new_name"].str.split("/").str[-1]

        for i in range(len(df)):
            d = df.iloc[i]
            original_name = d["original_name"]
            new_name = d["new_name"]
            src = os.path.join(fastqs_dir, original_name)
            dest = os.path.join(fastqs_dir, new_name)
            if os.path.exists(src):
                os.rename(src, dest)
            else:
                print(original_name + " not found")
