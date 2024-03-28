# Adapted from https://github.com/nf-core/rnaseq/blob/3.14.0/bin/filter_gtf.py

import os
import sys
import re
import statistics
from typing import Set

def extract_fasta_seq_names(fasta_name: str) -> Set[str]:
    """Extracts the sequence names from a FASTA file."""
    with open(fasta_name) as fasta:
        return {line[1:].split(None, 1)[0] for line in fasta if line.startswith(">")}

def tab_delimited(file: str) -> float:
    """Check if file is tab-delimited and return median number of tabs."""
    with open(file, "r") as f:
        data = f.read(102400)
        return statistics.median(line.count("\t") for line in data.split("\n"))

def filter_gtf(fasta: str, gtf_in: str, filtered_gtf_out: str, skip_transcript_id_check: bool) -> None:
    """Filter GTF file based on FASTA sequence names."""
    if tab_delimited(gtf_in) != 8:
        raise ValueError("Invalid GTF file: Expected 9 tab-separated columns.")
    seq_names_in_genome = extract_fasta_seq_names(fasta)
    print(f"Extracted chromosome sequence names from {fasta}")
    print("All sequence IDs from FASTA: " + ", ".join(sorted(seq_names_in_genome)))
    seq_names_in_gtf = set()
    try:
        with open(gtf_in) as gtf, open(filtered_gtf_out, "w") as out:
            line_count = 0
            for line in gtf:
                seq_name = line.split("\t")[0]
                seq_names_in_gtf.add(seq_name)  # Add sequence name to the set
                if seq_name in seq_names_in_genome:
                    if skip_transcript_id_check or re.search(r'transcript_id "([^"]+)"', line):
                        out.write(line)
                        line_count += 1
            if line_count == 0:
                raise ValueError("All GTF lines removed by filters")
    except IOError as e:
        print(f"File operation failed: {e}")
        return

    print("All sequence IDs from GTF: " + ", ".join(sorted(seq_names_in_gtf)))
    print(f"Extracted {line_count} matching sequences from {gtf_in} into {filtered_gtf_out}")

filter_gtf(par["fasta"], par["gtf"], par["filtered_gtf"], par["skip_transcript_id_check"])

# Get versions
text = f"{meta['functionality_name']}:\n  python: {sys.version.split()[0]}"
if par['versions'] and os.path.exists(par['versions']):
    with open(par['versions'], 'a') as f:
        f.write(text + '\n')
    os.rename(par['versions'], par['updated_versions'])
else:
    with open(par['updated_versions'], 'w') as f:
        f.write(text + '\n')