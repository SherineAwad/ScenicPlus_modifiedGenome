#!/usr/bin/env python3
import os
import argparse
import pickle
import numpy as np
from pycisTopic.lda_models import run_cgs_models_mallet

def run_mallet_on_object(cistopic_obj, mallet_path, n_topics, n_cpu, n_iter,
                         tmp_path, save_path, mallet_memory, random_state,
                         alpha, alpha_by_topic, eta, eta_by_topic):
    """
    Run MALLET on a single merged CistopicObject and select the best model.
    """
    sample_name = "merged_sample"
    print(f"Running MALLET for {sample_name}...")

    # Ensure tmp and save directories exist
    sample_tmp_path = os.path.join(tmp_path, sample_name)
    sample_save_path = os.path.join(save_path, sample_name)
    os.makedirs(sample_tmp_path, exist_ok=True)
    os.makedirs(sample_save_path, exist_ok=True)

    # Set Java memory for Mallet
    os.environ['MALLET_MEMORY'] = mallet_memory

    # Run MALLET CGS models
    models = run_cgs_models_mallet(
        cistopic_obj,
        n_topics=[n_topics] if isinstance(n_topics, int) else n_topics,
        n_cpu=n_cpu,
        n_iter=n_iter,
        random_state=random_state,
        alpha=alpha,
        alpha_by_topic=alpha_by_topic,
        eta=eta,
        eta_by_topic=eta_by_topic,
        tmp_path=sample_tmp_path,
        save_path=sample_save_path,
        mallet_path=mallet_path,
    )

    if not models or len(models) == 0:
        raise ValueError("No MALLET models were generated.")

    # Select best model by log likelihood
    best_model = max(models, key=lambda m: getattr(m, 'log_likelihood', -1e10))
    
    # FIX: Align model vocabulary with original DTM
    best_model = align_model_vocabulary(best_model, cistopic_obj)

    # Save best model immediately
    best_model_file = os.path.join(sample_save_path, "best_model.pkl")
    with open(best_model_file, "wb") as f:
        pickle.dump(best_model, f)
    print(f"[INFO] Saved best LDA model -> {best_model_file}")

    # Attach best model to object
    cistopic_obj.add_LDA_model(best_model)
    cistopic_obj.selected_model = best_model

    print(f"Finished MALLET. Model added to merged CistopicObject.\n")
    return cistopic_obj

def align_model_vocabulary(model, cistopic_obj):
    """
    Align MALLET model vocabulary with original DTM vocabulary.
    MALLET may filter out some regions/features, causing shape mismatch.
    """
    # Get original vocabulary from cistopic object
    original_vocab = cistopic_obj.region_names
    original_vocab_set = set(original_vocab)
    
    # Get model's current vocabulary (from MALLET)
    # Check what attributes are available in the model
    model_vocab = None
    
    # Try different possible attribute names for vocabulary
    if hasattr(model, 'vocabulary'):
        model_vocab = model.vocabulary
    elif hasattr(model, 'feature_names'):
        model_vocab = model.feature_names
    elif hasattr(model, 'region_names'):
        model_vocab = model.region_names
    elif hasattr(model, 'model') and hasattr(model.model, 'vocabulary'):
        model_vocab = list(model.model.vocabulary)
    
    if model_vocab is None:
        print("[WARNING] Could not find vocabulary in model, skipping alignment")
        return model
    
    model_vocab_set = set(model_vocab)
    
    # Check if alignment is needed
    if len(original_vocab) == len(model_vocab) and set(original_vocab) == set(model_vocab):
        print(f"[INFO] Vocabulary already aligned: {len(original_vocab)} regions")
        return model
    
    print(f"[INFO] Aligning vocabulary: original={len(original_vocab)}, model={len(model_vocab)}")
    print(f"[INFO] Missing in model: {len(original_vocab_set - model_vocab_set)} regions")
    print(f"[INFO] Extra in model: {len(model_vocab_set - original_vocab_set)} regions")
    
    # Create mapping from model vocabulary to indices
    model_vocab_to_idx = {region: idx for idx, region in enumerate(model_vocab)}
    
    # Get topic-word distribution (beta matrix)
    if hasattr(model, 'beta'):
        beta_matrix = model.beta  # Shape: (n_topics, n_features_in_model)
    elif hasattr(model, 'topic_word_distrib'):
        beta_matrix = model.topic_word_distrib
    else:
        print("[WARNING] Could not find topic-word distribution, skipping alignment")
        return model
    
    n_topics = beta_matrix.shape[0]
    
    # Create aligned beta matrix with zeros for missing regions
    aligned_beta = np.zeros((n_topics, len(original_vocab)))
    
    for orig_idx, region in enumerate(original_vocab):
        if region in model_vocab_to_idx:
            model_idx = model_vocab_to_idx[region]
            aligned_beta[:, orig_idx] = beta_matrix[:, model_idx]
    
    # Update the model's beta matrix
    if hasattr(model, 'beta'):
        model.beta = aligned_beta
    if hasattr(model, 'topic_word_distrib'):
        model.topic_word_distrib = aligned_beta
    
    # Update vocabulary/feature names
    if hasattr(model, 'vocabulary'):
        model.vocabulary = original_vocab
    if hasattr(model, 'feature_names'):
        model.feature_names = original_vocab
    if hasattr(model, 'region_names'):
        model.region_names = original_vocab
    
    print(f"[INFO] Vocabulary aligned successfully")
    return model

def main(cistopic_obj_pickle, mallet_path, n_topics, n_cpu, n_iter,
         tmp_path, save_path, mallet_memory="300G", random_state=555,
         alpha=5.0, alpha_by_topic=True, eta=0.1, eta_by_topic=True):

    # Load merged CistopicObject
    with open(cistopic_obj_pickle, "rb") as f:
        cistopic_obj = pickle.load(f)

    # Run MALLET
    updated_obj = run_mallet_on_object(
        cistopic_obj=cistopic_obj,
        mallet_path=mallet_path,
        n_topics=n_topics,
        n_cpu=n_cpu,
        n_iter=n_iter,
        tmp_path=tmp_path,
        save_path=save_path,
        mallet_memory=mallet_memory,
        random_state=random_state,
        alpha=alpha,
        alpha_by_topic=alpha_by_topic,
        eta=eta,
        eta_by_topic=eta_by_topic
    )

    # Save updated object
    output_pickle = os.path.join(save_path, "merged_cistopic_with_models.pkl")
    with open(output_pickle, "wb") as f:
        pickle.dump(updated_obj, f)
    print(f"Updated merged CistopicObject saved to: {output_pickle}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MALLET topic modeling on a merged CistopicObject")
    parser.add_argument("--cistopic_obj_pickle", required=True)
    parser.add_argument("--mallet_path", required=True)
    parser.add_argument("--n_topics", type=int, nargs="+", required=True)
    parser.add_argument("--n_cpu", type=int, default=12)
    parser.add_argument("--n_iter", type=int, default=500)
    parser.add_argument("--tmp_path", required=True)
    parser.add_argument("--save_path", required=True)
    parser.add_argument("--mallet_memory", default="300G")
    parser.add_argument("--random_state", type=int, default=555)
    parser.add_argument("--alpha", type=float, default=5.0)
    parser.add_argument("--alpha_by_topic", action="store_true", default=True)
    parser.add_argument("--eta", type=float, default=0.1)
    parser.add_argument("--eta_by_topic", action="store_true", default=True)
    args = parser.parse_args()

    main(
        cistopic_obj_pickle=args.cistopic_obj_pickle,
        mallet_path=args.mallet_path,
        n_topics=args.n_topics,
        n_cpu=args.n_cpu,
        n_iter=args.n_iter,
        tmp_path=args.tmp_path,
        save_path=args.save_path,
        mallet_memory=args.mallet_memory,
        random_state=args.random_state,
        alpha=args.alpha,
        alpha_by_topic=args.alpha_by_topic,
        eta=args.eta,
        eta_by_topic=args.eta_by_topic
    )
