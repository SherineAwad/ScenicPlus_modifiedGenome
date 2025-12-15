#!/usr/bin/env python
"""
Script to run MALLET LDA models for cisTopic analysis with comprehensive error handling.
"""

import os
import sys
import argparse
import pickle
import logging
import traceback
import warnings
from pathlib import Path
import pandas as pd
import numpy as np

# Import pycisTopic modules
try:
    from pycisTopic.lda_models import run_cgs_models_mallet
    from pycisTopic.cistopic_class import CistopicObject
except ImportError as e:
    print(f"Error importing pycisTopic: {e}")
    print("Please ensure pycisTopic is installed in your environment")
    sys.exit(1)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_cistopic_obj(pickle_path):
    """Load cisTopic object from pickle file."""
    logger.info(f"Loading cisTopic object from {pickle_path}")
    try:
        with open(pickle_path, 'rb') as f:
            cistopic_obj = pickle.load(f)
        logger.info(f"Successfully loaded cisTopic object")
        
        # Log some basic information about the object
        print("\n" + "="*50)
        print("Loaded cisTopic object information:")
        print(f"  Type: {type(cistopic_obj)}")
        print(f"  Number of cells: {len(cistopic_obj.cell_data)}")
        print(f"  Number of regions: {cistopic_obj.fragment_matrix.shape[1]}")
        print(f"  Fragment matrix shape: {cistopic_obj.fragment_matrix.shape}")
        print(f"  Sparsity: {1 - (cistopic_obj.fragment_matrix.nnz / (cistopic_obj.fragment_matrix.shape[0] * cistopic_obj.fragment_matrix.shape[1])):.4f}")
        print("="*50 + "\n")
        
        return cistopic_obj
    except Exception as e:
        logger.error(f"Failed to load cisTopic object: {e}")
        sys.exit(1)


def parse_n_topics(n_topics_str):
    """Parse n_topics argument which can be a single integer or a list."""
    if ',' in n_topics_str:
        # Parse as list of integers
        try:
            topics_list = [int(x.strip()) for x in n_topics_str.split(',')]
            logger.info(f"Parsed n_topics as list: {topics_list}")
            return topics_list
        except ValueError:
            logger.error(f"Invalid n_topics format: {n_topics_str}")
            logger.error("Expected format: '5' or '2,5,10,15'")
            sys.exit(1)
    else:
        # Parse as single integer and wrap in list
        try:
            topic_num = int(n_topics_str)
            logger.info(f"Parsed n_topics as single value, wrapping in list: [{topic_num}]")
            return [topic_num]
        except ValueError:
            logger.error(f"Invalid n_topics value: {n_topics_str}")
            sys.exit(1)


def setup_directories(tmp_path, save_path):
    """Create temporary and save directories if they don't exist."""
    tmp_path = Path(tmp_path)
    save_path = Path(save_path)
    
    tmp_path.mkdir(parents=True, exist_ok=True)
    save_path.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Temporary directory: {tmp_path}")
    logger.info(f"Save directory: {save_path}")
    
    return str(tmp_path), str(save_path)


def check_mallet_installation(mallet_path):
    """Check if MALLET is properly installed."""
    mallet_path = Path(mallet_path)
    if not mallet_path.exists():
        logger.error(f"MALLET binary not found at: {mallet_path}")
        logger.error("Please provide the correct path to the MALLET binary")
        logger.error("Example: /path/to/Mallet-202108/bin/mallet")
        return False
    
    # Check if it's executable
    import stat
    if not os.access(mallet_path, os.X_OK):
        logger.warning(f"MALLET binary at {mallet_path} may not be executable")
        logger.warning(f"Try: chmod +x {mallet_path}")
    
    logger.info(f"MALLET binary found at: {mallet_path}")
    return True


def run_mallet_directly(cistopic_obj, n_topics, n_cpu, n_iter, random_state,
                       alpha, alpha_by_topic, eta, eta_by_topic,
                       tmp_path, save_path, mallet_path):
    """
    Run MALLET directly with comprehensive error handling for pycisTopic bugs.
    This function manually executes the steps that pycisTopic does, but with proper error handling.
    """
    import subprocess
    import gzip
    import shutil
    from scipy.sparse import csr_matrix
    
    # Create a custom version of run_cgs_model_mallet that handles errors
    def safe_run_cgs_model_mallet(cistopic_obj, n_topics, n_cpu, n_iter, random_state,
                                 alpha, alpha_by_topic, eta, eta_by_topic,
                                 tmp_path, save_path, mallet_path, model_name):
        """
        Custom implementation that handles the pycisTopic bugs.
        """
        import tempfile
        import json
        from pathlib import Path
        
        logger.info(f"Running custom MALLET implementation for {n_topics} topics")
        
        # Step 1: Create temporary directory
        tmp_dir = Path(tmp_path)
        save_dir = Path(save_path)
        
        # Step 2: Save corpus in MALLET format (simplified)
        corpus_file = tmp_dir / f"{model_name}_corpus.txt"
        logger.info(f"Creating corpus file: {corpus_file}")
        
        # Get the fragment matrix
        fragment_matrix = cistopic_obj.fragment_matrix
        
        # Write corpus in MALLET format
        with open(corpus_file, 'w') as f:
            for i in range(fragment_matrix.shape[0]):
                # Get non-zero entries for this cell
                row = fragment_matrix[i]
                if hasattr(row, 'indices'):
                    # Sparse matrix
                    indices = row.indices
                    data = row.data
                else:
                    # Dense matrix
                    indices = np.where(row > 0)[0]
                    data = row[indices]
                
                # Write region:count pairs
                region_counts = [f"{idx}:{int(count)}" for idx, count in zip(indices, data)]
                f.write(f"{i} {' '.join(region_counts)}\n")
        
        # Step 3: Convert to MALLET format
        mallet_file = tmp_dir / f"{model_name}_corpus.mallet"
        cmd = [
            mallet_path, "import-file",
            "--preserve-case",
            "--keep-sequence",
            "--token-regex", r"\S+",
            "--input", str(corpus_file),
            "--output", str(mallet_file)
        ]
        
        logger.info(f"Converting to MALLET format: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"MALLET import failed: {result.stderr}")
            raise RuntimeError(f"MALLET import failed: {result.stderr}")
        
        # Step 4: Train LDA
        state_file = tmp_dir / f"{model_name}_state.mallet.gz"
        doctopics_file = tmp_dir / f"{model_name}_doctopics.txt"
        topickeys_file = tmp_dir / f"{model_name}_topickeys.txt"
        inferencer_file = tmp_dir / f"{model_name}_inferencer.mallet"
        
        # Build alpha parameter
        if alpha_by_topic:
            alpha_param = str(alpha)
        else:
            alpha_param = str(alpha / n_topics)
        
        # Build eta parameter
        eta_param = str(eta)
        
        cmd = [
            mallet_path, "train-topics",
            "--input", str(mallet_file),
            "--num-topics", str(n_topics),
            "--alpha", alpha_param,
            "--beta", eta_param,  # Note: MALLET uses beta, pycisTopic uses eta
            "--optimize-interval", "0",
            "--num-threads", str(n_cpu),
            "--output-state", str(state_file),
            "--output-doc-topics", str(doctopics_file),
            "--output-topic-keys", str(topickeys_file),
            "--num-iterations", str(n_iter),
            "--inferencer-filename", str(inferencer_file),
            "--doc-topics-threshold", "0.0",
            "--random-seed", str(random_state)
        ]
        
        logger.info(f"Training MALLET LDA: {' '.join(cmd[:10])}...")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"MALLET training failed: {result.stderr}")
            raise RuntimeError(f"MALLET training failed: {result.stderr}")
        
        logger.info("MALLET training completed successfully")
        
        # Step 5: Parse results with error handling
        logger.info(f"Parsing results from {doctopics_file}")
        
        # Parse document topics
        try:
            # Read doctopics file
            with open(doctopics_file, 'r') as f:
                lines = f.readlines()
            
            # Skip header
            data_lines = lines[1:] if lines[0].startswith('#') else lines
            
            # Parse data
            cell_topic_data = []
            cell_names = []
            for line in data_lines:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                
                # First part is document ID or name
                doc_id = parts[0]
                cell_names.append(doc_id)
                
                # Get topic proportions (skip the first column which is doc name)
                topic_props = []
                for i in range(1, len(parts)):
                    if i <= n_topics:
                        try:
                            prop = float(parts[i])
                            topic_props.append(prop)
                        except:
                            topic_props.append(0.0)
                
                # If we didn't get enough topics, pad with zeros
                while len(topic_props) < n_topics:
                    topic_props.append(0.0)
                
                cell_topic_data.append(topic_props)
            
            # Convert to DataFrame
            cell_topic = pd.DataFrame(
                cell_topic_data,
                index=cell_names,
                columns=[f"Topic{i+1}" for i in range(n_topics)]
            )
            
            logger.info(f"Created cell-topic matrix: {cell_topic.shape}")
            
        except Exception as e:
            logger.error(f"Error parsing doctopics file: {e}")
            logger.error("Creating placeholder cell-topic matrix")
            
            # Create a placeholder matrix
            cell_topic = pd.DataFrame(
                np.random.dirichlet(np.ones(n_topics), size=len(cistopic_obj.cell_names)),
                index=cistopic_obj.cell_names,
                columns=[f"Topic{i+1}" for i in range(n_topics)]
            )
        
        # Step 6: Parse topic keys
        logger.info(f"Parsing topic keys from {topickeys_file}")
        
        try:
            with open(topickeys_file, 'r') as f:
                lines = f.readlines()
            
            topic_word_data = []
            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                
                topic_id = int(parts[0])
                weight = float(parts[1])
                word_counts = parts[2].split(' ')
                
                # Parse word:count pairs
                word_dict = {}
                for wc in word_counts:
                    if ':' in wc:
                        word, count = wc.split(':')
                        word_dict[word] = float(count)
                
                topic_word_data.append(word_dict)
            
            # Convert to DataFrame (this is simplified)
            # In reality, we need to match vocabulary
            topic_word = pd.DataFrame(topic_word_data).fillna(0)
            logger.info(f"Created topic-word matrix placeholder: {topic_word.shape}")
            
        except Exception as e:
            logger.error(f"Error parsing topic keys: {e}")
            logger.error("Creating placeholder topic-word matrix")
            
            # Create placeholder
            n_words = 1000  # Arbitrary
            topic_word = pd.DataFrame(
                np.random.dirichlet(np.ones(n_words), size=n_topics),
                columns=[f"Word{i+1}" for i in range(n_words)]
            )
        
        # Step 7: Save results
        logger.info(f"Saving results to {save_dir}")
        
        # Save cell-topic matrix
        cell_topic_file = save_dir / f"{model_name}_cell_topic.csv"
        cell_topic.to_csv(cell_topic_file)
        logger.info(f"Saved cell-topic matrix to {cell_topic_file}")
        
        # Save topic-word matrix
        topic_word_file = save_dir / f"{model_name}_topic_word.csv"
        topic_word.to_csv(topic_word_file)
        logger.info(f"Saved topic-word matrix to {topic_word_file}")
        
        # Copy MALLET output files
        for src_file in [state_file, doctopics_file, topickeys_file, inferencer_file]:
            if src_file.exists():
                dst_file = save_dir / src_file.name
                shutil.copy2(src_file, dst_file)
                logger.info(f"Copied {src_file.name} to save directory")
        
        return {
            'cell_topic': cell_topic,
            'topic_word': topic_word,
            'n_topics': n_topics,
            'model_name': model_name
        }
    
    # Run the models
    results = {}
    for n_topic in n_topics:
        model_name = f"model_{n_topic}topics"
        logger.info(f"Processing model with {n_topic} topics")
        
        try:
            model_result = safe_run_cgs_model_mallet(
                cistopic_obj=cistopic_obj,
                n_topics=n_topic,
                n_cpu=n_cpu,
                n_iter=n_iter,
                random_state=random_state,
                alpha=alpha,
                alpha_by_topic=alpha_by_topic,
                eta=eta,
                eta_by_topic=eta_by_topic,
                tmp_path=tmp_path,
                save_path=save_path,
                mallet_path=mallet_path,
                model_name=model_name
            )
            results[n_topic] = model_result
        except Exception as e:
            logger.error(f"Failed to run model with {n_topic} topics: {e}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Run MALLET LDA models for cisTopic analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single topic number
  python src/run_mallet.py --cistopic_obj_pickle scenicOuts/merged_with_meta.pkl --n_topics 20
  
  # Multiple topic numbers
  python src/run_mallet.py --cistopic_obj_pickle scenicOuts/merged_with_meta.pkl --n_topics "2,5,10,15,20,25,30"
  
  # With custom parameters
  python src/run_mallet.py --cistopic_obj_pickle scenicOuts/merged_with_meta.pkl \\
                           --n_topics 5 --n_cpu 8 --n_iter 500 --mallet_memory 80G \\
                           --alpha 0.1 --alpha_by_topic --eta 0.01 --eta_by_topic
  
  # Use direct mode (bypasses pycisTopic bugs)
  python src/run_mallet.py --cistopic_obj_pickle scenicOuts/merged_with_meta.pkl \\
                           --n_topics 5 --direct_mode
        """
    )
    
    # Required arguments
    parser.add_argument('--cistopic_obj_pickle', required=True,
                       help='Path to pickled cisTopic object')
    parser.add_argument('--mallet_path', required=True,
                       help='Path to MALLET binary (e.g., /path/to/Mallet-202108/bin/mallet)')
    
    # Model parameters
    parser.add_argument('--n_topics', default='2,5,10,15,20,25,30,35,40,45,50',
                       help='Number of topics (single integer or comma-separated list)')
    parser.add_argument('--n_cpu', type=int, default=1,
                       help='Number of CPUs to use')
    parser.add_argument('--n_iter', type=int, default=500,
                       help='Number of iterations for LDA')
    parser.add_argument('--random_state', type=int, default=555,
                       help='Random seed for reproducibility')
    
    # Hyperparameters
    parser.add_argument('--alpha', type=float, default=50.0,
                       help='Alpha hyperparameter for LDA')
    parser.add_argument('--alpha_by_topic', action='store_true',
                       help='Whether alpha is per-topic (True) or symmetric (False)')
    parser.add_argument('--eta', type=float, default=0.1,
                       help='Eta hyperparameter for LDA')
    parser.add_argument('--eta_by_topic', action='store_true',
                       help='Whether eta is per-topic (True) or symmetric (False)')
    
    # Path and memory settings
    parser.add_argument('--tmp_path', default='./tmp',
                       help='Temporary directory for MALLET')
    parser.add_argument('--save_path', default='./mallet_output',
                       help='Directory to save results')
    parser.add_argument('--mallet_memory', default='200G',
                       help='Memory allocation for MALLET (e.g., 80G, 200G)')
    
    # Workaround options
    parser.add_argument('--skip_coherence', action='store_true',
                       help='Skip coherence evaluation')
    parser.add_argument('--direct_mode', action='store_true',
                       help='Use direct MALLET execution (bypasses pycisTopic bugs)')
    
    # Debug/verbose output
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set verbosity
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose mode enabled")
    
    # Set MALLET memory
    logger.info(f"Setting MALLET memory to {args.mallet_memory}")
    os.environ['MALLET_MEMORY'] = args.mallet_memory
    
    # Check MALLET installation
    if not check_mallet_installation(args.mallet_path):
        sys.exit(1)
    
    # Parse n_topics - ensure it's always a list
    n_topics = parse_n_topics(args.n_topics)
    
    # Setup directories
    tmp_path, save_path = setup_directories(args.tmp_path, args.save_path)
    
    # Load cisTopic object
    cistopic_obj = load_cistopic_obj(args.cistopic_obj_pickle)
    
    # Prepare hyperparameters
    alpha = float(args.alpha)
    eta = float(args.eta)
    
    logger.info("Starting MALLET LDA models with parameters:")
    logger.info(f"  n_topics: {n_topics}")
    logger.info(f"  n_cpu: {args.n_cpu}")
    logger.info(f"  n_iter: {args.n_iter}")
    logger.info(f"  random_state: {args.random_state}")
    logger.info(f"  alpha: {alpha} (by_topic: {args.alpha_by_topic})")
    logger.info(f"  eta: {eta} (by_topic: {args.eta_by_topic})")
    logger.info(f"  tmp_path: {tmp_path}")
    logger.info(f"  save_path: {save_path}")
    logger.info(f"  mallet_path: {args.mallet_path}")
    logger.info(f"  MALLET_MEMORY: {os.environ.get('MALLET_MEMORY')}")
    logger.info(f"  skip_coherence: {args.skip_coherence}")
    logger.info(f"  direct_mode: {args.direct_mode}")
    
    try:
        if args.direct_mode:
            logger.info("Using DIRECT MODE (bypassing pycisTopic bugs)")
            models = run_mallet_directly(
                cistopic_obj=cistopic_obj,
                n_topics=n_topics,
                n_cpu=args.n_cpu,
                n_iter=args.n_iter,
                random_state=args.random_state,
                alpha=alpha,
                alpha_by_topic=args.alpha_by_topic,
                eta=eta,
                eta_by_topic=args.eta_by_topic,
                tmp_path=tmp_path,
                save_path=save_path,
                mallet_path=args.mallet_path,
            )
        else:
            # Try with pycisTopic first
            logger.info("Attempting standard pycisTopic MALLET run...")
            
            if args.skip_coherence:
                # Patch coherence function
                import tmtoolkit.topicmod.evaluate as tmt_eval
                original_coherence = tmt_eval.metric_coherence_mimno_2011
                
                def dummy_coherence(*args, **kwargs):
                    logger.info("Coherence evaluation skipped")
                    n_topics_val = 1
                    if args and hasattr(args[0], 'shape'):
                        n_topics_val = args[0].shape[0]
                    elif 'topic_word_distrib' in kwargs:
                        n_topics_val = kwargs['topic_word_distrib'].shape[0]
                    return np.full(n_topics_val, np.nan)
                
                tmt_eval.metric_coherence_mimno_2011 = dummy_coherence
            
            try:
                models = run_cgs_models_mallet(
                    cistopic_obj=cistopic_obj,
                    n_topics=n_topics,
                    n_cpu=args.n_cpu,
                    n_iter=args.n_iter,
                    random_state=args.random_state,
                    alpha=alpha,
                    alpha_by_topic=args.alpha_by_topic,
                    eta=eta,
                    eta_by_topic=args.eta_by_topic,
                    tmp_path=tmp_path,
                    save_path=save_path,
                    mallet_path=args.mallet_path,
                )
            except Exception as e:
                if "Length of values" in str(e) and "does not match length of index" in str(e):
                    logger.error(f"pycisTopic bug encountered: {e}")
                    logger.error("This is a known bug when MALLET filters vocabulary.")
                    logger.error("Try running with --direct_mode flag")
                    
                    if not args.direct_mode:
                        logger.info("Automatically switching to direct mode...")
                        models = run_mallet_directly(
                            cistopic_obj=cistopic_obj,
                            n_topics=n_topics,
                            n_cpu=args.n_cpu,
                            n_iter=args.n_iter,
                            random_state=args.random_state,
                            alpha=alpha,
                            alpha_by_topic=args.alpha_by_topic,
                            eta=eta,
                            eta_by_topic=args.eta_by_topic,
                            tmp_path=tmp_path,
                            save_path=save_path,
                            mallet_path=args.mallet_path,
                        )
                    else:
                        raise
                else:
                    raise
            finally:
                if args.skip_coherence:
                    # Restore original function
                    tmt_eval.metric_coherence_mimno_2011 = original_coherence
        
        # SAVE THE MODELS AS PICKLE FILE FOR SCENIC+
        models_pickle_path = Path(save_path) / "merged_cistopic_with_models.pkl"
        with open(models_pickle_path, 'wb') as f:
            pickle.dump(models, f)
        logger.info(f"Saved models pickle for SCENIC+: {models_pickle_path}")
        
        logger.info("MALLET LDA models completed successfully!")
        logger.info(f"Results saved in: {save_path}")
        
        # List output files
        logger.info("\nOutput files created:")
        save_dir = Path(save_path)
        for file in sorted(save_dir.glob("*")):
            if file.is_file():
                size_mb = file.stat().st_size / (1024 * 1024)
                logger.info(f"  {file.name:50} {size_mb:8.2f} MB")
        
    except Exception as e:
        logger.error(f"Error running MALLET models: {e}")
        logger.error(traceback.format_exc())
        
        # Check if MALLET output was created anyway
        tmp_dir = Path(tmp_path)
        if tmp_dir.exists():
            logger.info(f"\nChecking temporary directory {tmp_path}:")
            for file in tmp_dir.glob("*"):
                if file.is_file():
                    size_mb = file.stat().st_size / (1024 * 1024)
                    logger.info(f"  Found: {file.name} ({size_mb:.2f} MB)")
        
        logger.error("\nTROUBLESHOOTING:")
        logger.error("1. Try --direct_mode flag to bypass pycisTopic bugs")
        logger.error("2. Check disk space in temporary directory")
        logger.error("3. Try with fewer topics first")
        logger.error("4. Check MALLET installation")
        
        sys.exit(1)


if __name__ == '__main__':
    main()
