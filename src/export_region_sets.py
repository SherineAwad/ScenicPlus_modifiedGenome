#!/usr/bin/env python3
import os
import pickle
import argparse
from pycisTopic.utils import region_names_to_coordinates

def export_region_sets(cistopic_pickle, out_dir):
    # Create folders
    os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok=True)

    # Load binarised CistopicObject
    with open(cistopic_pickle, "rb") as f:
        obj = pickle.load(f)

    print("[DEBUG] Available attributes in object:")
    for attr in dir(obj):
        if "region_bin" in attr or "markers" in attr:
            print(f"  {attr}")

    # --- Export Topics using Otsu threshold ---
    if hasattr(obj, "region_bin_topics_otsu"):
        print(f"[INFO] Exporting Topics_otsu from region_bin_topics_otsu")
        print(f"[DEBUG] region_bin_topics_otsu shape: {obj.region_bin_topics_otsu.shape}")

        # Correct: columns = topics, index = regions
        for topic in obj.region_bin_topics_otsu.columns:
            regions = obj.region_bin_topics_otsu.index[obj.region_bin_topics_otsu[topic] > 0]
            if len(regions) > 0:
                df = region_names_to_coordinates(regions)
                df.sort_values(["Chromosome", "Start", "End"]).to_csv(
                    os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
                    sep="\t", header=False, index=False
                )
                print(f"[INFO] Exported {len(regions)} regions for {topic}")
            else:
                print(f"[WARNING] No regions found for {topic}")

    # --- Export top-n Topics ---
    if hasattr(obj, "region_bin_topics_top"):
        print(f"[INFO] Exporting Topics_top_3k from region_bin_topics_top")
        print(f"[DEBUG] region_bin_topics_top shape: {obj.region_bin_topics_top.shape}")

        # Correct: columns = topics, index = regions
        for topic in obj.region_bin_topics_top.columns:
            regions = obj.region_bin_topics_top.index[obj.region_bin_topics_top[topic] > 0]
            if len(regions) > 0:
                df = region_names_to_coordinates(regions)
                df.sort_values(["Chromosome", "Start", "End"]).to_csv(
                    os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
                    sep="\t", header=False, index=False
                )
                print(f"[INFO] Exported {len(regions)} regions for {topic}")
            else:
                print(f"[WARNING] No regions found for {topic}")

    # --- Export DARs per cell type ---
    if hasattr(obj, "markers_dict"):
        print(f"[INFO] Exporting DARs from markers_dict")
        for cell_type, marker_data in obj.markers_dict.items():
            if not hasattr(marker_data, "index"):
                print(f"[WARNING] Skipping {cell_type} - not a DataFrame")
                continue

            regions = marker_data.index
            if len(regions) > 0:
                df = region_names_to_coordinates(regions)
                df.sort_values(["Chromosome", "Start", "End"]).to_csv(
                    os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
                    sep="\t", header=False, index=False
                )
                print(f"[INFO] Exported {len(regions)} DARs for {cell_type}")
            else:
                print(f"[WARNING] No regions found for {cell_type}")

    print(f"[INFO] Region sets exported to: {os.path.join(out_dir, 'region_sets')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export Topics and DARs region sets from a binarised CistopicObject")
    parser.add_argument("-i", "--input_pickle", required=True, help="Binarised CistopicObject pickle file")
    parser.add_argument("-o", "--out_dir", required=True, help="Output directory for BED files")
    args = parser.parse_args()

    export_region_sets(args.input_pickle, args.out_dir)

