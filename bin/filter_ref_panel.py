#!/usr/bin/env python3
"""
Filter a full reference panel VCF to keep:
  1. Variants at (CHROM, POS) positions present in the array positions file
  2. Structural variants whose ID starts with 'Sniffles2'

Reads the positions file into a set, then streams the ref VCF sequentially
(no random access / indexing required) and writes matching lines to stdout.

Usage:
    filter_ref_panel.py --positions all_array_positions.tsv --ref panel.vcf.gz \
        | bgzip -c > reduced_ref.vcf.gz
"""

import sys
import gzip
import argparse


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--positions', required=True,
                   help='TSV file with CHROM<tab>POS lines (all batches merged)')
    p.add_argument('--ref', required=True,
                   help='Full reference panel VCF or VCF.gz')
    args = p.parse_args()

    # Load all array positions into a set of (chrom, pos) tuples
    keep_positions = set()
    with open(args.positions) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t', 1)
            if len(parts) == 2:
                keep_positions.add((parts[0], parts[1]))

    sys.stderr.write(f"Loaded {len(keep_positions):,} unique array positions\n")

    # Stream ref VCF and write matching records to stdout
    n_written = 0
    opener = gzip.open if args.ref.endswith('.gz') else open
    with opener(args.ref, 'rt') as fin:
        for line in fin:
            if line.startswith('#'):
                sys.stdout.write(line)
                continue
            fields = line.split('\t', 4)
            chrom, pos, vid = fields[0], fields[1], fields[2]
            if (chrom, pos) in keep_positions or vid.startswith('Sniffles2'):
                sys.stdout.write(line)
                n_written += 1

    sys.stderr.write(f"Wrote {n_written:,} variants to reduced ref panel\n")


if __name__ == '__main__':
    main()
