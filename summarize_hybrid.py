import json
import pandas as pd
from pysradb.sraweb import SRAweb
import time
import sys
import argparse

def summarize_hybrid():
    parser = argparse.ArgumentParser(description="Summarize hybrid BioSamples.")
    parser.add_argument("input_file", nargs="?", default="hybrid_biosamples.json", help="Input JSON file path.")
    parser.add_argument("--output", default="hybrid_data_summary.tsv", help="Output TSV file path.")
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output

    try:
        with open(input_file, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {input_file} not found.")
        sys.exit(1)

    # Extract unique biosamples
    biosamples = list(set([entry.get('biosample') for entry in data if 'biosample' in entry]))
    print(f"Found {len(biosamples)} unique BioSamples.")

    db = SRAweb()
    results = []

    batch_size = 50
    max_retries = 3

    for i in range(0, len(biosamples), batch_size):
        batch = biosamples[i:i+batch_size]
        print(f"Processing batch {i//batch_size + 1}/{(len(biosamples) + batch_size - 1)//batch_size} (Samples {i} to {i+len(batch)})...")

        success = False
        df = None
        for attempt in range(max_retries):
            try:
                df = db.sra_metadata(batch, detailed=True)
                success = True
                break
            except Exception as e:
                print(f"  Attempt {attempt+1} failed: {e}")
                time.sleep(2 * (attempt + 1))

        if not success:
            print(f"  Failed to process batch after {max_retries} attempts. Skipping.")
            continue

        if df is not None and not df.empty:
            # Group by biosample to aggregate instruments and pick representative metadata
            for biosample, group in df.groupby('biosample'):
                # Sample Type
                sample_type = "N/A"
                if 'organism_name' in group.columns:
                    sample_type = group['organism_name'].dropna().iloc[0] if not group['organism_name'].dropna().empty else "N/A"

                # Environment
                env = "N/A"
                # Priority list for environment
                env_cols = ['env_local_scale', 'env_broad_scale', 'isolation_source', 'env_medium', 'sample_name', 'study_title']
                for col in env_cols:
                    if col in group.columns:
                        vals = group[col].dropna().astype(str).tolist()
                        # Filter out empty strings, "nan", "not applicable", "missing"
                        valid_vals = [v for v in vals if v.lower() not in ['nan', '', 'not applicable', 'missing', 'none']]
                        if valid_vals:
                            env = valid_vals[0]
                            break

                # Instruments
                instruments = "N/A"
                if 'instrument_model' in group.columns:
                    inst_list = sorted(list(set(group['instrument_model'].dropna().astype(str).tolist())))
                    instruments = ", ".join(inst_list)

                results.append({
                    "BioSample ID": biosample,
                    "Sample Type": sample_type,
                    "Environment": env,
                    "Instruments": instruments
                })

        time.sleep(1) # Rate limit courtesy

    if not results:
        print("No results collected.")
        return

    # Create DataFrame and save
    summary_df = pd.DataFrame(results)

    # Ensure columns order
    cols = ["BioSample ID", "Sample Type", "Environment", "Instruments"]
    # Add any missing cols with N/A
    for c in cols:
        if c not in summary_df.columns:
            summary_df[c] = "N/A"

    summary_df = summary_df[cols]

    # Remove duplicates just in case
    summary_df = summary_df.drop_duplicates(subset=["BioSample ID"])

    summary_df.to_csv(output_file, sep='\t', index=False)
    print(f"Summary saved to {output_file}")

if __name__ == "__main__":
    summarize_hybrid()
