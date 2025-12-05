#!/usr/bin/env python3
"""
FIXED VERSION - Add scRNA metadata to cisTopicObject with proper barcode handling
"""

import pickle
import pandas as pd
import re
import scipy.sparse as sp
import argparse

def add_scrna_metadata_fixed(cistopic_pickle, scrna_csv, output_pickle):
    """
    FIXED VERSION: Add scRNA metadata with consistent barcode handling
    """
    print(f"🔹 Loading cisTopic object: {cistopic_pickle}")
    with open(cistopic_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    print(f"Original object: {cistopic_obj.cell_data.shape[0]} cells, {len(cistopic_obj.region_names)} regions")

    # --- Load and clean scRNA metadata ---
    scrna_meta = pd.read_csv(scrna_csv)
    if 'barcode' not in scrna_meta.columns or 'celltype' not in scrna_meta.columns:
        raise ValueError("❌ scRNA CSV must contain 'barcode' and 'celltype' columns")

    # Clean scRNA barcodes
    scrna_meta['barcode'] = scrna_meta['barcode'].astype(str).str.strip()
    scrna_meta['barcode_clean'] = scrna_meta['barcode'].str.replace('-1$', '', regex=True)
    scrna_meta.set_index('barcode', inplace=True)

    # --- Prepare ATAC barcodes ---
    original_barcodes = cistopic_obj.cell_data.index
    cistopic_obj.cell_data = cistopic_obj.cell_data.copy()

    # Add clean barcodes as a column, but keep original index
    cistopic_obj.cell_data['barcode_clean'] = original_barcodes.str.replace('-1$', '', regex=True)
    cistopic_obj.cell_data['orig_barcode'] = original_barcodes  # Keep original barcodes

    # --- COMPREHENSIVE BARCODE DEBUG ---
    print("\n" + "="*80)
    print("COMPREHENSIVE BARCODE DEBUG")
    print("="*80)

    print(f"ATAC total cells: {len(original_barcodes)}")
    print(f"scRNA total cells: {len(scrna_meta)}")

    # Check for duplicates in original index
    duplicate_mask = original_barcodes.duplicated(keep=False)
    duplicate_count = duplicate_mask.sum()
    if duplicate_count > 0:
        print(f"❌ FOUND DUPLICATE BARCODES IN ORIGINAL INDEX: {duplicate_count}")
        print("Sample duplicates:")
        # Convert to list and show first 10 unique duplicates
        duplicate_values = original_barcodes[duplicate_mask].unique()[:10]
        for i, barcode in enumerate(duplicate_values):
            print(f"  {i+1}: '{barcode}'")

        # Show how many times each duplicate appears
        duplicate_freq = original_barcodes.value_counts()
        duplicate_freq = duplicate_freq[duplicate_freq > 1]
        print(f"\nDuplicate frequency:")
        for barcode, count in duplicate_freq.head(5).items():
            print(f"  '{barcode}': {count} times")

        # FIX: Handle duplicates by making them unique while preserving format
        print("\n🔧 Making duplicate barcodes unique...")
        seen = {}
        new_barcodes = []
        for barcode in original_barcodes:
            if barcode not in seen:
                seen[barcode] = 1
                new_barcodes.append(barcode)
            else:
                seen[barcode] += 1
                new_barcodes.append(f"{barcode}_dup{seen[barcode]}")
        
        # Update the index with unique barcodes
        cistopic_obj.cell_data.index = new_barcodes
        cistopic_obj.cell_data['barcode'] = new_barcodes
        print(f"✅ Made barcodes unique while preserving format")

    # --- Join scRNA metadata using clean barcodes ---
    scrna_clean_to_celltype = dict(zip(scrna_meta['barcode_clean'], scrna_meta['celltype']))
    cistopic_obj.cell_data['celltype'] = cistopic_obj.cell_data['barcode_clean'].map(scrna_clean_to_celltype)

    # --- Report match statistics ---
    n_matched = cistopic_obj.cell_data['celltype'].notna().sum()
    total_cells = cistopic_obj.cell_data.shape[0]
    print(f"\n✅ Matched {n_matched} / {total_cells} cells ({n_matched / total_cells * 100:.2f}%)")

    # --- Filter matrices to only keep matched cells ---
    print("\nFiltering matrices to match cell data...")
    keep_cell_indices = cistopic_obj.cell_data['celltype'].notna()
    print(f"Keeping {keep_cell_indices.sum()} cells in matrices")

    # Filter cell_data
    cistopic_obj.cell_data = cistopic_obj.cell_data[keep_cell_indices]

    # Filter matrices
    if hasattr(cistopic_obj, 'fragment_matrix') and cistopic_obj.fragment_matrix is not None:
        if sp.issparse(cistopic_obj.fragment_matrix):
            cistopic_obj.fragment_matrix = cistopic_obj.fragment_matrix[keep_cell_indices, :]
        else:
            cistopic_obj.fragment_matrix = cistopic_obj.fragment_matrix[keep_cell_indices, :]
        print(f"✓ Filtered fragment_matrix to shape: {cistopic_obj.fragment_matrix.shape}")

    if hasattr(cistopic_obj, 'binary_matrix') and cistopic_obj.binary_matrix is not None:
        if sp.issparse(cistopic_obj.binary_matrix):
            cistopic_obj.binary_matrix = cistopic_obj.binary_matrix[keep_cell_indices, :]
        else:
            cistopic_obj.binary_matrix = cistopic_obj.binary_matrix[keep_cell_indices, :]
        print(f"✓ Filtered binary_matrix to shape: {cistopic_obj.binary_matrix.shape}")

    # --- FIXED: KEEP ORIGINAL BARCODES - DO NOT REPLACE WITH cell_000000 ---
    print("\n" + "="*80)
    print("PRESERVING ORIGINAL BARCODES")
    print("="*80)
    
    # Synchronize cell_names with the (potentially deduplicated) cell_data index
    cistopic_obj.cell_names = cistopic_obj.cell_data.index.tolist()
    
    print(f"✅ Preserved original barcode format")
    print(f"Sample barcodes: {cistopic_obj.cell_names[:5]}")

    # --- FINAL VERIFICATION ---
    print("\n" + "="*80)
    print("FINAL VERIFICATION")
    print("="*80)

    print(f"cell_data shape: {cistopic_obj.cell_data.shape}")
    print(f"cell_data index unique: {cistopic_obj.cell_data.index.is_unique}")
    print(f"cell_names length: {len(cistopic_obj.cell_names)}")
    print(f"cell_data index matches cell_names: {cistopic_obj.cell_data.index.tolist() == cistopic_obj.cell_names}")

    if hasattr(cistopic_obj, 'fragment_matrix'):
        print(f"fragment_matrix shape: {cistopic_obj.fragment_matrix.shape}")
        print(f"Matrix rows match cell_data: {cistopic_obj.fragment_matrix.shape[0] == cistopic_obj.cell_data.shape[0]}")

    if hasattr(cistopic_obj, 'binary_matrix'):
        print(f"binary_matrix shape: {cistopic_obj.binary_matrix.shape}")
        print(f"Binary matrix rows match cell_data: {cistopic_obj.binary_matrix.shape[0] == cistopic_obj.cell_data.shape[0]}")

    # Check celltype distribution
    celltype_counts = cistopic_obj.cell_data['celltype'].value_counts()
    print(f"\n📊 FINAL CELLTYPE DISTRIBUTION:")
    for celltype, count in celltype_counts.items():
        print(f"  {celltype}: {count} cells")

    # --- Save cleaned object ---
    with open(output_pickle, "wb") as f:
        pickle.dump(cistopic_obj, f)

    print("\n💾 UPDATED AND CLEANED CistopicObject saved to:", output_pickle)
    print(f"📊 Final cell_data shape: {cistopic_obj.cell_data.shape}")
    if hasattr(cistopic_obj, 'fragment_matrix'):
        print(f"📊 Final fragment_matrix shape: {cistopic_obj.fragment_matrix.shape}")
    if hasattr(cistopic_obj, 'binary_matrix'):
        print(f"📊 Final binary_matrix shape: {cistopic_obj.binary_matrix.shape}")
    print(f"📊 Barcodes preserved: {cistopic_obj.cell_names[:3]}")

    print("🎉 SUCCESS: Object is now properly formatted for topic modeling!")
    print("✅ Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FIXED: Add scRNA metadata to merged CistopicObject")
    parser.add_argument("--cistopic_pickle", required=True, help="Merged cisTopic pickle file")
    parser.add_argument("--scrna_csv", required=True, help="scRNA CSV with 'barcode' and 'celltype'")
    parser.add_argument("--output_pickle", required=True, help="Output path for updated pickle")
    args = parser.parse_args()

    add_scrna_metadata_fixed(args.cistopic_pickle, args.scrna_csv, args.output_pickle)
