#!/usr/bin/env python3

import argparse
import math
import random
import re
import signal
import subprocess
import sys

import pysam
import pandas as pd

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def regions(bed_path):
    """Iterator over a BED file."""
    if not bed_path: # will iterate over the whole BAM if BED not specified
        yield None, None, None
    else:
        with open(bed_path, "r") as f:
            for line in f:
                chrom, start, end, *_ = line.split()
                yield chrom, int(start), int(end)


def pileup(bam, ref, bed):
    """Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    """
    pileups = []
    for chrom, start, end in bed:
        for col in bam.pileup(chrom, start, end):
            # check first if the current column lies within this region (pysam
            # pileup returns whole reads overlapping a queried region, not
            # just bases in this region)
            if chrom and start and end and not (chrom == col.reference_name and start <= col.pos < end):
                continue
            else:
                bases = [b.upper() for b in col.get_query_sequences() if b.upper() in "ACGT"]
                if bases:
                    pileups.append((
                        col.reference_name,
                        col.reference_pos + 1,
                        ref.fetch(col.reference_name, col.reference_pos, col.reference_pos + 1),
                        bases
                    ))
    return pd.DataFrame(pileups, columns=["chrom", "pos", "ref", "pileup"])


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
            print(f"{i.chrom}\t{i.pos}\t.\t{i.ref}\t{i.base}\t.\t.\t.\tGT\t{int(i.ref != i.base)}", file=vcf)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sample random alleles from a given BAM file")
    parser.add_argument("--bam", help="BAM file to sample from", required=True)
    parser.add_argument("--ref", help="FASTA reference sequence", required=True)
    parser.add_argument("--bed", help="BED file with coordinates of regions/sites to sample")
    parser.add_argument("--output", help="Output filename or EIGENSTRAT prefix", required=True)
    parser.add_argument("--format", help="Output format", choices=["VCF", "EIGENSTRAT", "pileup"], required=True)
    parser.add_argument("--sample-name", help="Sample name to put in VCF/EIGENSTRAT", required=True)

    args = parser.parse_args()

    if args.format in ["VCF", "EIGENSTRAT"] and not args.sample_name:
        parser.error(f"Sample name has to be specified when writing {args.format}")

    bam = pysam.AlignmentFile(args.bam)
    ref = pysam.FastaFile(args.ref)
    bed = regions(args.bed)

    if not bam.has_index():
        print("An indexed BAM file is required, please run 'samtools index' first", file=sys.stderr)
        sys.exit(1)

    pileups = pileup(bam, ref, bed)
    pileups["base"] = pileups["pileup"].apply(lambda x: random.choice(x))

    if args.format == "VCF":
        write_vcf(args.output, pileups, args.sample_name)
    elif args.format == "EIGENSTRAT":
        write_eigenstrat(args.output, pileups, args.sample_name)
    else:
        print(pileups)