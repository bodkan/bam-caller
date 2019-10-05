#!/usr/bin/env python3

import argparse
import functools
import math
import random
import signal
import sys
import os.path
from collections import Counter

import pysam
import pandas as pd

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def tolerance(i, tol):
    """Calculate the number of alleles required to agree at a site
    for a consensus to pass."""
    return i - math.trunc((i - 1) * tol)


def consensus(bases, tol):
    """Call consensus on piled-up bases."""
    counts = Counter(bases).most_common()
    base = counts[0][0]
    count = counts[0][1]
    if count >= tolerance(len(bases), tol):
        return base
    else:
        return "N"

def majority(bases):
    """Select the most common allele. Handling cases of multiple
    alleles of the same proportion relies on the fact that the
    method `most_common()` returns such elements in an arbitrary
    order.
    """
    return Counter(bases).most_common()[0][0]

def flush(i, calls, out_fun):
    """Save genotype calls to a file (either a VCF or a pileup TSV)."""
    print(f"\r{i + 1} positions processed", end="")
    calls = pd.DataFrame(
        calls, columns=["chrom", "pos", "ref", "coverage", "call"]
    ).query('ref != "N" & call != "N"')
    out_fun(calls)


def get_ref_base(col):
    tuples = col.pileups[0].alignment.get_aligned_pairs(with_seq=True)
    for read_pos, ref_pos, ref_base in tuples:
        if ref_pos == col.reference_pos:
            return ref_base.upper()


def call_bases(call_fun, out_fun, bam, mincov, minbq, minmq, chrom):
    """Sample bases in a given region of the genome based on the pileupx
    of reads. If no coordinates were specified, sample from the whole BAM file.
    """
    calls = []
    for i, col in enumerate(bam.pileup(contig=chrom, compute_baq=False,
                                       min_base_quality=minbq,
                                       min_mapping_quality=minmq)):
        bases = col.get_query_sequences(add_indels=True)

        # filter out sites with no reads and sites with indels
        if bases and "*" not in bases and all(len(i) == 1 for i in bases):
            bases = [b.upper() for b in col.get_query_sequences()
                     if b and b.upper() in "ACGT"]

            if len(bases) >= mincov:
                calls.append((
                    col.reference_name,
                    col.reference_pos + 1,
                    get_ref_base(col),
                    len(bases),
                    call_fun(bases)
                ))

        if i % 1000000 == 0:
            flush(i, calls, out_fun)
            calls = []

    flush(i, calls, out_fun)


def write_vcf(calls, output, sample_name):
    filename = output + ".vcf"
    new_file = not os.path.isfile(filename)
    with open(filename, "w" if new_file else "a") as vcf:
        if new_file:
            print(
                "##fileformat=VCFv4.1\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}".
                format(sample=sample_name), file=vcf
            )
        for i in calls.itertuples():
            alt, gt = (".", 0) if i.ref == i.call else (i.call, 1)
            print(f"{i.chrom}\t{i.pos}\t.\t{i.ref}\t{alt}\t.\t.\t.\tGT:DP\t{gt}:{i.coverage}", file=vcf)


def write_pileup(pileups, output):
    filename = output + ".txt"
    new_file = not os.path.isfile(filename)
    with open(filename, "w" if new_file else "a") as tsv:
        if new_file: print("chrom\tpos\tref\tpileup\tA\tC\tG\tT", file=tsv)
        for i in pileups.itertuples():
            counts = Counter(i.call)
            counts_str = '\t'.join(str(counts[i]) for i in 'ACGT')
            print(f"{i.chrom}\t{i.pos}\t{i.ref}\t{''.join(i.call)}\t{counts_str}", file=tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sample random alleles from a given BAM file")
    parser.add_argument("--bam", help="BAM file to sample from", required=True)
    parser.add_argument("--chrom", help="Chromosome to sample from")
    parser.add_argument("--strategy", help="How to 'genotype'?", choices=["random", "consensus", "majority", "pileup"], required=True)
    parser.add_argument("--tolerance", help="What proportion of discordant alleles to allow for consensus?", type=float, default=0.0)
    parser.add_argument("--mincov", help="Minimum coverage", type=int, default=1)
    parser.add_argument("--minbq", help="Minimum base quality", type=int, default=13)
    parser.add_argument("--minmq", help="Minimum read mapping quality", type=int, default=0)
    parser.add_argument("--sample-name", help="Sample name to put in a VCF header")
    parser.add_argument("--output", help="Output file prefix")

    args = parser.parse_args()

    if args.strategy != "pileup" and not args.sample_name:
        parser.error(f"Sample name has to be specified when writing a VCF file")

    bam = pysam.AlignmentFile(args.bam)

    if not bam.has_index():
        print("BAM file index is missing", file=sys.stderr)
        sys.exit(1)

    if args.strategy == "pileup":
        call_fun = lambda x: "".join(x)
        out_fun = functools.partial(write_pileup, output=args.output)
    elif args.strategy == "random":
        call_fun = lambda x: random.choice(x)
        out_fun = functools.partial(write_vcf, output=args.output,
                                    sample_name=args.sample_name)
    elif args.strategy == "consensus":
        call_fun = functools.partial(consensus, tol=args.tolerance)
        out_fun = functools.partial(write_vcf, output=args.output,
                                    sample_name=args.sample_name)
    elif args.strategy == "majority":
        call_fun = majority
        out_fun = functools.partial(write_vcf, output=args.output,
                                    sample_name=args.sample_name)

    call_bases(call_fun, out_fun, bam, args.mincov,
               args.minbq, args.minmq, args.chrom)
