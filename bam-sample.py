#!/usr/bin/env python3

import argparse
import random
import re
import signal
import subprocess
import sys
from collections import Counter

import pysam
import pandas as pd

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def pileup(bam, ref):
    """Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    """
    pileups = []
    for col in bam.pileup():
        bases = [b.upper() for b in col.get_query_sequences() if b and b.upper() in "ACGT"]
        if bases:
            pileups.append((
                col.reference_name,
                col.reference_pos + 1,
                ref.fetch(col.reference_name, col.reference_pos, col.reference_pos + 1),
                bases
            ))
    pileups = pd.DataFrame(pileups, columns=["chrom", "pos", "ref", "pileup"])
    pileups["coverage"] = pileups.pileup.apply(lambda x: len(x))

    return pileups


def write_vcf(output, sites, sample_name):
    filename = re.sub(".gz", "", output) if output.endswith(".gz") else output
    with open(filename, "w") as vcf:
        print(
            "##fileformat=VCFv4.1\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}".
            format(sample=sample_name), file=vcf
        )
        for i in sites.itertuples():
            alt, gt = (".", 0) if i.ref == i.base else (i.base, 1)
            print(f"{i.chrom}\t{i.pos}\t.\t{i.ref}\t{alt}\t.\t.\t.\tGT\t{gt}", file=vcf)

    if filename != output:
        subprocess.run(["bgzip", "-f", filename])
        subprocess.run(["tabix", "-f", output])


def write_eigenstrat(output, sites, sample_name):
    with open(output + ".ind", "w") as ind, \
         open(output + ".snp", "w") as snp, \
         open(output + ".geno", "w") as geno:
        print(f"{sample_name}\tU\t{sample_name}", file=ind)
        for i in sites.itertuples():
            print(f".\t{i.chrom}\t0.0\t{i.pos}\t{i.ref}\t{i.base}", file=snp)
            print(f"{2 * int(i.ref == i.base)}", file=geno)


def write_pileup(output, pileups):
    with open(output, "w") as tsv:
        for i in pileups.itertuples():
            print(f"{i.chrom}\t{i.pos}\t{i.ref}\t{''.join(i.pileup)}", file=tsv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sample random alleles from a given BAM file")
    parser.add_argument("--bam", help="BAM file to sample from", required=True)
    parser.add_argument("--ref", help="FASTA reference sequence", required=True)
    parser.add_argument("--strategy", help="How to 'genotype'?", choices=["random", "consensus"], required=True)
    parser.add_argument("--sample-name", help="Sample name to put in VCF/EIGENSTRAT")
    parser.add_argument("--format", help="Output format", choices=["vcf", "eigenstrat", "pileup"], required=True)
    parser.add_argument("--output", help="Output filename or EIGENSTRAT prefix", required=True)

    args = parser.parse_args()

    if args.format in ["vcf", "eigenstrat"]:
        if not args.sample_name:
            parser.error(f"Sample name has to be specified when writing {args.format}")
        if args.strategy == "consensus":
            parser.error(f"Only random-calling strategy allowed for {args.format} output format")


    bam = pysam.AlignmentFile(args.bam)
    ref = pysam.FastaFile(args.ref)

    if not bam.has_index():
        print("An indexed BAM file is required, please run 'samtools index' first", file=sys.stderr)
        sys.exit(1)

    pileups = pileup(bam, ref)

    if args.strategy == "random":
        pileups["base"] = pileups["pileup"].apply(lambda x: random.choice(x))
    elif args.strategy == "consensus":
        base_counts = pileups.pileup.apply(lambda x: len(Counter(x)))
        pileups = pileups[base_counts == 1]

    if args.format == "eigenstrat":
        write_eigenstrat(args.output, pileups, args.sample_name)
    elif args.format == "vcf":
        write_vcf(args.output, pileups, args.sample_name)
    else:
        write_pileup(args.output, pileups)
