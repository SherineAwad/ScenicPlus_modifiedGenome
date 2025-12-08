#!/usr/bin/env python3
import argparse
import pickle
import re

def inspect_cistopic_object(obj, idx=None):
    prefix = f"Object {idx}:" if idx is not None else ""
    print(f"\n===== Inspecting {prefix} =====")

    # Object type
    print("Object type:", type(obj))

    # 1. Fragment matrix info
    if hasattr(obj, "fragment_matrix") and obj.fragment_matrix is not None:
        fm = obj.fragment_matrix
        print("\nFragment matrix info:")
        print("Shape (cells x regions):", fm.shape)
        print("Number of cells:", len(getattr(obj, "cell_names", [])))
        print("Number of regions/features:", len(getattr(obj, "region_names", [])))
        if hasattr(fm, "nnz"):
            print("Number of non-zero entries:", fm.nnz)
            print("Sparsity:", 1 - fm.nnz / (fm.shape[0] * fm.shape[1]))

        # 2. Sample of region names
        if hasattr(obj, "region_names") and obj.region_names is not None:
            regions = list(obj.region_names)
            print("\nSample regions/features (first 50):", regions[:50])
            # Check for problematic characters
            pattern = re.compile(r'[^A-Za-z0-9._:-]')
            bad_names = [r for r in regions if pattern.search(r)]
            print(f"Number of regions with problematic characters: {len(bad_names)}")
            if len(bad_names) > 0:
                print("Sample problematic region names:", bad_names[:20])
        else:
            print("No 'region_names' found.")

        # 3. Sample of cell names
        if hasattr(obj, "cell_names") and obj.cell_names is not None:
            cells = list(obj.cell_names)
            print("\nSample cells (first 50):", cells[:50])
        else:
            print("No 'cell_names' found.")
    else:
        print("No 'fragment_matrix' found in this object.")

    # 4. Metadata / cell_data
    if hasattr(obj, "cell_data") and obj.cell_data is not None:
        print("\nCell metadata / cell_data columns:", list(obj.cell_data.columns))
        print("Sample metadata for first 5 cells:")
        print(obj.cell_data.head())
    else:
        print("\nNo 'cell_data' found.")

    # 5. Other relevant attributes
    extras = ["fragments", "fragment_coords", "fragment_names"]
    for attr in extras:
        if hasattr(obj, attr):
            val = getattr(obj, attr)
            print(f"\nExtra field '{attr}': type {type(val)}", end="")
            try:
                print(", sample (first 10):", list(val)[:10] if hasattr(val, "__len__") else val)
            except Exception as e:
                print(f", could not print sample: {e}")

def main(cistopic_pickle):
    # Load object
    with open(cistopic_pickle, "rb") as f:
        loaded = pickle.load(f)

    # Check if it's a list
    if isinstance(loaded, list):
        print(f"Loaded a list of {len(loaded)} CistopicObjects.")
        for i, obj in enumerate(loaded):
            inspect_cistopic_object(obj, i)
    else:
        inspect_cistopic_object(loaded)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect CistopicObject(s) for MALLET debugging")
    parser.add_argument("cistopic_obj_pickle", help="Path to CistopicObject pickle (or list of objects)")
    args = parser.parse_args()

    main(args.cistopic_obj_pickle)

