#!/usr/bin/env python3

import argparse
from collections import defaultdict

import pysam
import pandas as pd
import  matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


BASES = "ACGT"

def dict_to_dataframe(mismatches, len_limit):
    """Convert dictionary of substitution counts at different read positions
    into a tidy DataFrame.
    """
    counts = defaultdict(list)
    for mismatch in [f"{b1}{b2}" for b1 in BASES for b2 in BASES if b1 != b2]:
        for pos in range(len_limit):
            counts[mismatch].append(mismatches[(mismatch, pos)])
    df = pd.DataFrame(counts)
    df.index = df.index + 1
    return df


def process_read(mismatches, ref_bases, read_bases, len_limit):
    """Calculate the number of mismatches along a read and update
    the mismatch table. Do not analyze sites within a read which are
    further from the beginning than len_limit.
    """
    for pos, (b_ref, b_read) in enumerate(zip(ref_bases, read_bases)):
        # skip the rest of the read if already beyond the limit
        if pos >= len_limit: break

        # skip the position if it's not A, C, G or T
        if b_ref not in BASES or b_read not in BASES: continue

        # if there' a mismatch on this position, increment the counter
        if b_ref != b_read:
            mismatches[(b_ref + b_read, pos)] += 1
    return mismatches


def count_mismatches(bam, len_limit=20):
    """Count substitution patterns in a BAM file and return the results
    as a pair of pandas DataFrames (one for the forward direction, another
    one for the reverse direction).
    """   
    # initialize counters of reads and mismatch-counting data frames
    n_fwd_reads, n_rev_reads = 0, 0
    fwd_mismatches = defaultdict(int)
    rev_mismatches = defaultdict(int)
    
    for i, read in enumerate(bam, 1):
        ref_bases = read.get_reference_sequence().upper()
        read_bases = read.query_alignment_sequence.upper()
        mismatches = fwd_mismatches

        # if the read is in a reversed orientation, reverse the fetched
        # sequences and change reference to the table to be updated
        if read.is_reverse:
            read_bases = reversed(read_bases)
            ref_bases = reversed(ref_bases)
            mismatches = rev_mismatches

        process_read(mismatches, ref_bases, read_bases, len_limit)
        
        if read.is_reverse:
            n_rev_reads += 1
        else:
            n_fwd_reads += 1
        
        if i % 100000 == 0: print(f"{i} reads processed", end="\r")
    
    print(f"\nTotal number of forward reads: {n_fwd_reads}")
    print(f"Total number of reverse reads: {n_rev_reads}")

    fwd_df = dict_to_dataframe(fwd_mismatches, len_limit) / n_fwd_reads
    rev_df = dict_to_dataframe(rev_mismatches, len_limit) / n_rev_reads

    fwd_df.insert(0, "strand", "forward")
    rev_df.insert(0, "strand", "reverse")

    return fwd_df, rev_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", help="Path to a BAM file", required=True)
    parser.add_argument("--len_limit", help="How far into each read to look for damage?", type=int, default=20)
    parser.add_argument("--all-mismatches", help="Report all mismatches?", action="store_true")
    parser.add_argument("--figure", help="Path to an output figure (suffix determines format - png, svg, etc.)")
    parser.add_argument("--dpi", help="DPI of an output figure", type=int, default=300)
    parser.add_argument("--table", help="Path to the output figure")
    args = parser.parse_args()

    if not (args.figure or args.table):
        parser.error("At least one output option must be specified (--figure or --table)")

    bam = pysam.AlignmentFile(args.bam)
    fwd_mismatches, rev_mismatches = count_mismatches(bam, args.len_limit)

    # keep only the major aDNA mismatches
    if not args.all_mismatches:
        fwd_mismatches = fwd_mismatches[["strand", "CT", "GA"]] 
        rev_mismatches = rev_mismatches[["strand", "CT", "GA"]] 

    # write results as a TSV file
    if args.table:
        pd.concat([
            fwd_mismatches.reset_index().rename(columns={"index": "pos"}),
            rev_mismatches.reset_index().rename(columns={"index": "pos"})
        ]).to_csv(args.table, sep="\t", index=False, float_format="%.5f")

    # plot the results
    if args.figure:
        rev_mismatches.index = -rev_mismatches.index

        fig, axes = plt.subplots(nrows=1, ncols=2)

        axes[0].xaxis.set_major_locator(MaxNLocator(integer=True))
        axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
        axes[0].set_xlim(xmin=1, xmax=args.len_limit)
        axes[1].set_xlim(xmin=-args.len_limit, xmax=-1)

        fwd_mismatches.plot(ax=axes[0], figsize=(15, 7), title="Forward reads", )
        rev_mismatches.plot(ax=axes[1], title="Reverse reads")

        axes[0].set(xlabel="Position in the read", ylabel="Proportion of reads with a mismatch")
        axes[1].set(xlabel="Position in the read", ylabel="Proportion of reads with a mismatch")

        plt.savefig(args.figure, dpi=args.dpi)