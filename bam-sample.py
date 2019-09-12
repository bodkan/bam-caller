#!/usr/bin/env python3

import argparse
import functools
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


def consensus(bases, tol):
    """Call consensus on piled-up bases."""
    counts = Counter(bases).most_common()
    base = counts[0][0]
    count = counts[0][1]
    if count >= tolerance(len(bases), tol):
        return base
    else:
        return "N"


def call_bases(fun, bam, ref, minbq, minmq, chrom):
    """Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    """
    calls = []
    for i, col in enumerate(bam.pileup(contig=chrom, compute_baq=False, min_base_quality=minbq, min_mapping_quality=minmq)):
        bases = col.get_query_sequences(add_indels=True)
        # filter out sites with no reads and sites with indels or close to indels
        if bases and "*" not in bases and all(len(i) == 1 for i in bases):
            bases = [b.upper() for b in col.get_query_sequences() if b and b.upper() in "ACGT"]
            calls.append((
                col.reference_name,
                col.reference_pos + 1,
                ref.fetch(col.reference_name, col.reference_pos, col.reference_pos + 1),
                len(bases),
                fun(bases)
            ))
        if i % 1000000 == 0: print(f"\r{i + 1} positions processed", end="")
    print()
    calls = pd.DataFrame(calls, columns=["chrom", "pos", "ref", "coverage", "call"]).query('ref != "N" & call != "N"')

    return calls


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
            alt, gt = (".", 0) if i.ref == i.call else (i.call, 1)
            print(f"{i.chrom}\t{i.pos}\t.\t{i.ref}\t{alt}\t.\t.\t.\tGT:DP\t{gt}:{i.coverage}", file=vcf)


def write_eigenstrat(output, sites, sample_name):
    with open(output + ".ind", "w") as ind, \
         open(output + ".snp", "w") as snp, \
         open(output + ".geno", "w") as geno:
        print(f"{sample_name}\tU\t{sample_name}", file=ind)
        for i in sites.itertuples():
            print(f".\t{i.chrom}\t0.0\t{i.pos}\t{i.ref}\t{i.call}", file=snp)
            print(f"{2 * int(i.ref == i.call)}", file=geno)


def write_pileup(output, pileups):
    with open(output + ".txt", "w") as tsv:
        print("chrom\tpos\tref\tpileup\tA\tC\tG\tT", file=tsv)
        for i in pileups.itertuples():
            counts = Counter(i.call)
            counts_str = '\t'.join(str(counts[i]) for i in 'ACGT')
            print(f"{i.chrom}\t{i.pos}\t{i.ref}\t{''.join(i.call)}\t{counts_str}", file=tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sample random alleles from a given BAM file")
    parser.add_argument("--bam", help="BAM file to sample from", required=True)
    parser.add_argument("--ref", help="FASTA reference sequence", required=True)
    parser.add_argument("--chrom", help="Chromosome to sample from")
    parser.add_argument("--strategy", help="How to 'genotype'?", choices=["random", "consensus", "pileup"])
    parser.add_argument("--tolerance", help="What proportion of discordant alleles to allow for consensus?", type=float, default=0.0)
    parser.add_argument("--mincov", help="Minimum coverage", type=int, default=0)
    parser.add_argument("--minbq", help="Minimum base quality", type=int, default=13)
    parser.add_argument("--minmq", help="Minimum read mapping quality", type=int, default=0)
    parser.add_argument("--sample-name", help="Sample name to put in VCF/EIGENSTRAT")
    parser.add_argument("--format", help="Output formats", choices=["vcf", "eigenstrat"])
    parser.add_argument("--output", help="Output file prefix")

    args = parser.parse_args()
    #args = parser.parse_args("--bam ../ychr/data/bam/exome_den4.bam --ref /mnt/solexa/Genomes/hg19_evan/whole_genome.fa --chrom Y --strategy consensus --minbq 20 --minbq 25 --sample-name den4 --format vcf --output test_output".split())

    if args.format in ["vcf", "eigenstrat"] and not args.sample_name:
            parser.error(f"Sample name has to be specified when writing {args.format}")

    if args.format in ["vcf", "eigenstrat"] and not args.strategy:
            parser.error(f"Sampling strategy has to be specified when writing {args.format}")

    if args.format in ["vcf", "eigenstrat"] and args.strategy == "pileup":
        parser.error("VCF or EIGENSTRAT output does not make sense for pileup")

    bam = pysam.AlignmentFile(args.bam)
    ref = pysam.FastaFile(args.ref)

    if not bam.has_index():
        print("An indexed BAM file is required, please run 'samtools index' first", file=sys.stderr)
        sys.exit(1)

    if args.strategy == "pileup":
        fun = lambda x: "".join(x)
    elif args.strategy == "random":
        fun = lambda x: random.choice(x)
    elif args.strategy == "consensus":
        fun = functools.partial(consensus, tol=args.tolerance)

    calls = call_bases(fun, bam, ref, args.minbq, args.minmq, args.chrom).query(f"coverage >= {args.mincov}")

    if args.strategy == "pileup":
        write_pileup(args.output, calls)
    elif args.format == "eigenstrat":
        write_eigenstrat(args.output, calls, args.sample_name)
    elif args.format == "vcf":
        write_vcf(args.output, calls, args.sample_name)
