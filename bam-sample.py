#!/usr/bin/env python3

import argparse
import math
import random
import signal
import sys
from collections import Counter

import pysam
import pandas as pd

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def tolerance(i, tol):
    """Calculate the number of alleles required to agree at a site
    for a consensus to pass."""
    return i - math.trunc((i - 1) * tol)

def pileup(bam, ref, minbq, minmq):
    """Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    """
    pileups = []
    for i, col in enumerate(bam.pileup(compute_baq=False, min_base_quality=minbq, min_mapping_quality=minmq)):
        bases = col.get_query_sequences(add_indels=True)
        # filter out sites with no reads and sites with indels or close to indels
        if bases and "*" not in bases and all(len(i) == 1 for i in bases):
            bases = [b.upper() for b in col.get_query_sequences() if b and b.upper() in "ACGT"]
            pileups.append((
                col.reference_name,
                col.reference_pos + 1,
                ref.fetch(col.reference_name, col.reference_pos, col.reference_pos + 1),
                bases
            ))
        if i % 50000 == 0: print(f"\r{i + 1} positions processed", end="")
    print()
    pileups = pd.DataFrame(pileups, columns=["chrom", "pos", "ref", "pileup"]).query('ref != "N"')
    pileups["coverage"] = pileups.pileup.apply(lambda x: len(x))

    return pileups


def write_vcf(output, sites, sample_name):
    with open(output + ".vcf", "w") as vcf:
        print(
            "##fileformat=VCFv4.1\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}".
            format(sample=sample_name), file=vcf
        )
        for i in sites.itertuples():
            alt, gt = (".", 0) if i.ref == i.base else (i.base, 1)
            print(f"{i.chrom}\t{i.pos}\t.\t{i.ref}\t{alt}\t.\t.\t.\tGT:DP\t{gt}:{i.coverage}", file=vcf)


def write_eigenstrat(output, sites, sample_name):
    with open(output + ".ind", "w") as ind, \
         open(output + ".snp", "w") as snp, \
         open(output + ".geno", "w") as geno:
        print(f"{sample_name}\tU\t{sample_name}", file=ind)
        for i in sites.itertuples():
            print(f".\t{i.chrom}\t0.0\t{i.pos}\t{i.ref}\t{i.base}", file=snp)
            print(f"{2 * int(i.ref == i.base)}", file=geno)


def write_pileup(output, pileups):
    with open(output + ".txt", "w") as tsv:
        print("chrom\tpos\tref\tpileup\tA\tC\tG\tT", file=tsv)
        for i in pileups.itertuples():
            counts = Counter(i.pileup)
            counts_str = '\t'.join(str(counts[i]) for i in 'ACGT')
            print(f"{i.chrom}\t{i.pos}\t{i.ref}\t{''.join(i.pileup)}\t{counts_str}", file=tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sample random alleles from a given BAM file")
    parser.add_argument("--bam", help="BAM file to sample from", required=True)
    parser.add_argument("--ref", help="FASTA reference sequence", required=True)
    parser.add_argument("--strategy", help="How to 'genotype'?", choices=["random", "consensus"])
    parser.add_argument("--tolerance", help="What proportion of discordant alleles to allow for consensus?", type=float, default=0.0)
    parser.add_argument("--mincov", help="Minimum coverage", type=int, default=0)
    parser.add_argument("--minbq", help="Minimum base quality", type=int, default=13)
    parser.add_argument("--minmq", help="Minimum read mapping quality", type=int, default=0)
    parser.add_argument("--sample-name", help="Sample name to put in VCF/EIGENSTRAT")
    parser.add_argument("--format", help="Output formats", nargs="+", choices=["vcf", "eigenstrat", "pileup"], required=True)
    parser.add_argument("--output", help="Output file prefix", required=True)

    args = parser.parse_args()

    if args.format in ["vcf", "eigenstrat"] and not args.sample_name:
            parser.error(f"Sample name has to be specified when writing {args.format}")

    if args.format in ["vcf", "eigenstrat"] and not args.strategy:
            parser.error(f"Sampling strategy has to be specified when writing {args.format}")

    bam = pysam.AlignmentFile(args.bam)
    ref = pysam.FastaFile(args.ref)

    if not bam.has_index():
        print("An indexed BAM file is required, please run 'samtools index' first", file=sys.stderr)
        sys.exit(1)

    pileups = pileup(bam, ref, args.minbq, args.minmq).query(f"coverage >= {args.mincov}")

    if "pileup" in args.format:
        write_pileup(args.output, pileups)

    if args.strategy == "random":
        pileups["base"] = pileups["pileup"].apply(lambda x: random.choice(x))
    elif args.strategy == "consensus":
        pileups["base"] = pileups.pileup.apply(lambda x: Counter(x).most_common()[0][0])
        pileups["count"] = pileups.pileup.apply(lambda x: Counter(x).most_common()[0][1])
        pileups["tolerance"] = pileups.coverage.apply(tolerance, args = (args.tolerance, ))
        pileups = pileups[pileups["count"] >= pileups["tolerance"]].drop(["count", "tolerance"], axis=1)

    if "eigenstrat" in args.format:
        write_eigenstrat(args.output, pileups, args.sample_name)
    if "vcf" in args.format:
        write_vcf(args.output, pileups, args.sample_name)
