#!/usr/bin/env python3
"""
Script to identify BioSamples with both Long Reads and Short Reads data
from a list of studies provided in a JSON file.
Uses multiprocessing to speed up metadata fetching.
"""

import gzip
import json
import logging
import time
import argparse
import pandas as pd
from typing import List, Dict, Any
from multiprocessing import Pool, cpu_count
from pysradb.sraweb import SRAweb

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("find_hybrid_samples.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def load_studies(filepath: str) -> List[str]:
    """Load unique study accessions from the gzipped JSON file."""
    logger.info(f"Loading studies from {filepath}...")
    try:
        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            data = json.load(f)

        studies = set()
        for entry in data:
            if 'study_accession' in entry:
                studies.add(entry['study_accession'])

        logger.info(f"Found {len(studies)} unique studies.")
        return list(studies)
    except Exception as e:
        logger.error(f"Error loading studies: {e}")
        return []

def classify_platform(instrument_model: str) -> str:
    """Classify instrument model into 'LONG', 'SHORT', or 'OTHER'."""
    if not isinstance(instrument_model, str):
        return 'OTHER'

    model = instrument_model.upper()

    # Long Reads
    if any(x in model for x in ['NANOPORE', 'MINION', 'GRIDION', 'PROMETHION', 'PACBIO', 'SEQUEL']):
        return 'LONG'

    # Short Reads
    if any(x in model for x in ['ILLUMINA', 'HISEQ', 'MISEQ', 'NEXTSEQ', 'NOVASEQ', 'ION TORRENT', 'BGISEQ', 'DNBSEQ', 'SOLID', '454', 'AB 5500', 'HELIOS']):
        return 'SHORT'

    return 'OTHER'

def process_batch_worker(studies: List[str]) -> List[Dict]:
    """Worker function to process a batch of studies."""
    # Create a new SRAweb instance for each worker
    db = SRAweb()
    hybrid_samples = []

    # Retry logic
    max_retries = 3
    for attempt in range(max_retries):
        try:
            # Fetch metadata for all studies in the batch
            df = db.sra_metadata(studies, detailed=True)
            break # Success
        except Exception as e:
            if attempt < max_retries - 1:
                sleep_time = 2 * (attempt + 1)
                time.sleep(sleep_time)
            else:
                logger.error(f"Failed to process batch {studies[:3]}... after {max_retries} attempts: {e}")
                return []

    try:
        if df is None or df.empty:
            return []

        # Ensure required columns exist
        required_cols = ['sample_accession', 'run_accession', 'instrument_model', 'study_accession']
        for col in required_cols:
            if col not in df.columns:
                if col == 'instrument_model' and 'instrument' in df.columns:
                    df['instrument_model'] = df['instrument']
                else:
                    return []

        # Group by sample_accession
        for sample_acc, group in df.groupby('sample_accession'):
            if pd.isna(sample_acc) or sample_acc == 'N/A':
                continue

            long_reads = []
            short_reads = []

            for _, row in group.iterrows():
                platform = classify_platform(row['instrument_model'])
                run_info = {
                    'run_accession': row['run_accession'],
                    'instrument_model': row['instrument_model'],
                    'study_accession': row['study_accession']
                }

                if platform == 'LONG':
                    long_reads.append(run_info)
                elif platform == 'SHORT':
                    short_reads.append(run_info)

            if long_reads and short_reads:
                hybrid_samples.append({
                    'biosample': sample_acc,
                    'short_reads': short_reads,
                    'long_reads': long_reads,
                    'study_accession': list(set(r['study_accession'] for r in long_reads + short_reads))
                })

    except Exception as e:
        logger.error(f"Error processing batch {studies[:3]}...: {e}")

    return hybrid_samples

def main():
    parser = argparse.ArgumentParser(description="Find hybrid BioSamples in SRA studies.")
    parser.add_argument("--limit", type=int, help="Limit the number of studies to process (for testing).")
    parser.add_argument("--workers", type=int, default=4, help="Number of worker processes.")
    args = parser.parse_args()

    input_file = "data_metagenome.json.gz"
    output_file = "hybrid_biosamples.json"
    batch_size = 50

    studies = load_studies(input_file)
    if not studies:
        logger.error("No studies found. Exiting.")
        return

    if args.limit:
        logger.info(f"Limiting to first {args.limit} studies.")
        studies = studies[:args.limit]

    all_hybrid_samples = []

    # Prepare batches
    batches = [studies[i:i + batch_size] for i in range(0, len(studies), batch_size)]
    total_batches = len(batches)

    logger.info(f"Processing {len(studies)} studies in {total_batches} batches using {args.workers} workers...")

    start_time = time.time()

    try:
        with Pool(processes=args.workers) as pool:
            for i, result in enumerate(pool.imap_unordered(process_batch_worker, batches)):
                if result:
                    all_hybrid_samples.extend(result)

                # Log progress every batch
                logger.info(f"Processed {i + 1}/{total_batches} batches. Found {len(all_hybrid_samples)} hybrid samples so far.")

                # Incremental save every 5 batches
                if (i + 1) % 5 == 0:
                     try:
                        with open(output_file, 'w') as f:
                            json.dump(all_hybrid_samples, f, indent=2)
                        logger.info(f"Incremental save to {output_file}")
                     except Exception as e:
                        logger.error(f"Error saving incremental results: {e}")

    except KeyboardInterrupt:
        logger.warning("Interrupted by user. Saving partial results...")
    except Exception as e:
        logger.error(f"An error occurred: {e}")

    end_time = time.time()
    duration = end_time - start_time
    logger.info(f"Processing complete in {duration:.2f} seconds.")
    logger.info(f"Total hybrid samples found: {len(all_hybrid_samples)}")

    # Save results
    try:
        with open(output_file, 'w') as f:
            json.dump(all_hybrid_samples, f, indent=2)
        logger.info(f"Results saved to {output_file}")
    except Exception as e:
        logger.error(f"Error saving results: {e}")

if __name__ == "__main__":
    main()
