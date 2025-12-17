#!/usr/bin/env python3
import argparse
import os
import pickle
import numpy as np
import pandas as pd
from pycisTopic.topic_binarization import binarize_topics
import matplotlib.pyplot as plt
from types import SimpleNamespace


def run_binarizations(input_pickle, output_dir, ntop, scale_factor):
    os.makedirs(output_dir, exist_ok=True)

    # --- Load object ---
    with open(input_pickle, "rb") as f:
        obj = pickle.load(f)
    print(f"[INFO] Loaded object from {input_pickle}")

    # FIX: Ensure cell_topic exists at object level for pycisTopic compatibility
    if hasattr(obj, 'selected_model') and isinstance(obj.selected_model, dict):
        if 'cell_topic' in obj.selected_model and not hasattr(obj, 'cell_topic'):
            print("[INFO] Adding cell_topic directly to object for pycisTopic compatibility")
            obj.cell_topic = obj.selected_model['cell_topic']

    # --- Convert selected_model dict to object for pycisTopic ---
    if isinstance(obj.selected_model, dict):
        print("[INFO] Converting selected_model dict to object for pycisTopic compatibility")

        # Create a SimpleNamespace object with same keys as attributes
        class ModelObject:
            def __init__(self, model_dict):
                # Copy all dictionary items to object attributes
                for key, value in model_dict.items():
                    setattr(self, key, value)

        # Store original dict for later
        original_model_dict = obj.selected_model.copy()

        # Create object version
        obj.selected_model = ModelObject(original_model_dict)

        # Also attach as attribute for compatibility
        obj._selected_model_dict = original_model_dict

    # --- Optional scaling ---
    if scale_factor is not None and scale_factor > 0:
        if hasattr(obj.selected_model, 'topic_region'):
            obj.selected_model.topic_region *= scale_factor
            print(f"[INFO] Scaled topic_region by {scale_factor}. New max: {obj.selected_model.topic_region.values.max()}")
        elif isinstance(obj.selected_model, dict) and 'topic_region' in obj.selected_model:
            obj.selected_model['topic_region'] *= scale_factor
            print(f"[INFO] Scaled topic_region by {scale_factor}. New max: {obj.selected_model['topic_region'].values.max()}")

    # --- Ensure selected model exists ---
    if not hasattr(obj, "selected_model") or obj.selected_model is None:
        raise ValueError("[ERROR] No selected_model found in object. Run MALLET first.")

    # Check cell_topic exists
    cell_topic_exists = (
        hasattr(obj.selected_model, 'cell_topic') and obj.selected_model.cell_topic is not None
    ) or (
        isinstance(obj.selected_model, dict) and 'cell_topic' in obj.selected_model and obj.selected_model['cell_topic'] is not None
    )

    if not cell_topic_exists:
        raise ValueError("[ERROR] selected_model missing cell_topic. Check MALLET run.")

    print("[INFO] LDA model detected. Proceeding to binarization...")

    # --- 1) Binarize regions (ntop) ---
    print(f"[INFO] Running region binarization (top {ntop})...")
    try:
        region_bin_topics_top_dict = binarize_topics(
            obj, method='ntop', ntop=ntop, plot=True, num_columns=5
        )
        plt.savefig(os.path.join(output_dir, f"region_bin_top{ntop}.png"), bbox_inches='tight')
        plt.close()
        print(f"[INFO] Region binarization (ntop) completed successfully")
    except Exception as e:
        print(f"[ERROR] Region binarization (ntop) failed: {e}")
        print("[INFO] Creating empty result for ntop binarization")
        region_bin_topics_top_dict = {}

    region_bin_topics_top_combined = _combine_topic_matrices(region_bin_topics_top_dict, obj.region_names)
    print(f"[INFO] Top-n combined matrix shape: {region_bin_topics_top_combined.shape}")

    # --- 2) Binarize regions (Otsu) ---
    print("[INFO] Running region binarization (Otsu)...")
    try:
        region_bin_topics_otsu_dict = binarize_topics(
            obj, method='otsu', plot=True, num_columns=5
        )
        plt.savefig(os.path.join(output_dir, "region_bin_otsu.png"), bbox_inches='tight')
        plt.close()
        print(f"[INFO] Region binarization (Otsu) completed successfully")
    except Exception as e:
        print(f"[WARNING] Region binarization (Otsu) failed: {e}")
        region_bin_topics_otsu_dict = {}

    region_bin_topics_otsu_combined = _combine_topic_matrices(region_bin_topics_otsu_dict, obj.region_names)
    print(f"[INFO] Otsu combined matrix shape: {region_bin_topics_otsu_combined.shape}")

    # --- 3) Binarize cells (Li) ---
    print("[INFO] Running cell binarization (Li)...")
    try:
        binarized_cell_topic_dict = binarize_topics(
            obj, target='cell', method='li', plot=True, num_columns=5, nbins=100
        )
        plt.savefig(os.path.join(output_dir, "cell_topic_li.png"), bbox_inches='tight')
        plt.close()
        print(f"[INFO] Cell binarization (Li) completed successfully")
    except Exception as e:
        print(f"[WARNING] Cell binarization (Li) failed: {e}")
        binarized_cell_topic_dict = {}

    binarized_cell_topic_combined = _combine_cell_matrices(binarized_cell_topic_dict, obj.cell_names)
    print(f"[INFO] Cell-topic combined matrix shape: {binarized_cell_topic_combined.shape}")

    # --- Restore original selected_model dict if we converted it ---
    if hasattr(obj, '_selected_model_dict'):
        obj.selected_model = obj._selected_model_dict
        delattr(obj, '_selected_model_dict')

    # --- Attach binarized outputs to object ---
    obj.region_bin_topics_top = region_bin_topics_top_combined
    obj.region_bin_topics_otsu = region_bin_topics_otsu_combined
    obj.binarized_cell_topic = binarized_cell_topic_combined

    # --- Save object and outputs ---
    binarized_file = os.path.join(output_dir, "cistopic_obj_binarized.pkl")
    with open(binarized_file, "wb") as f:
        pickle.dump(obj, f)
    print(f"[INFO] Saved binarized CistopicObject -> {binarized_file}")

    binarized_cell_topic_file = os.path.join(output_dir, "binarized_cell_topic.pkl")
    with open(binarized_cell_topic_file, "wb") as f:
        pickle.dump(binarized_cell_topic_combined, f)
    print(f"[INFO] Saved binarized cell_topic DataFrame -> {binarized_cell_topic_file}")

    # --- Sanity check ---
    print("[INFO] Binarization complete. Summary:")
    print(f"  - region_bin_topics_top shape: {region_bin_topics_top_combined.shape}")
    print(f"  - region_bin_topics_otsu shape: {region_bin_topics_otsu_combined.shape}")
    print(f"  - binarized_cell_topic shape: {binarized_cell_topic_combined.shape}")


# --- Combine topic–region matrices ---
def _combine_topic_matrices(topic_dict, all_region_names):
    if not topic_dict:
        # Return empty DataFrame with correct structure
        return pd.DataFrame(index=all_region_names)

    topic_names = list(topic_dict.keys())
    combined_df = pd.DataFrame(0, index=all_region_names, columns=topic_names, dtype=int)

    for topic_name, topic_matrix in topic_dict.items():
        if isinstance(topic_matrix, pd.DataFrame) and not topic_matrix.empty:
            # typical shape: regions × 1
            if topic_matrix.shape[1] > 0:
                selected_regions = topic_matrix.index[topic_matrix.iloc[:, 0] > 0]
                combined_df.loc[selected_regions, topic_name] = 1

    return combined_df


# --- Combine topic–cell matrices ---
def _combine_cell_matrices(topic_dict, all_cell_names):
    if not topic_dict:
        # Return empty DataFrame with correct structure
        return pd.DataFrame(index=all_cell_names)

    topic_names = list(topic_dict.keys())
    combined_df = pd.DataFrame(0, index=all_cell_names, columns=topic_names, dtype=int)

    for topic_name, topic_matrix in topic_dict.items():
        if isinstance(topic_matrix, pd.DataFrame) and not topic_matrix.empty:
            # typical shape: cells × 1
            if topic_matrix.shape[1] > 0:
                selected_cells = topic_matrix.index[topic_matrix.iloc[:, 0] > 0]
                combined_df.loc[selected_cells, topic_name] = 1

    return combined_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run binarization of CistopicObject topics")
    parser.add_argument("-i", "--input_pickle", required=True, help="Input CistopicObject pickle")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    parser.add_argument("-n", "--ntop", type=int, default=3000, help="Top N regions for ntop method")
    parser.add_argument("-s", "--scale_topic_region", type=float, default=None, help="Optional scale factor")
    args = parser.parse_args()

    run_binarizations(args.input_pickle, args.output_dir, args.ntop, args.scale_topic_region)
