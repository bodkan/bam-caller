#!/usr/bin/env python3

import sys
import argparse
import random
from collections import Counter
from itertools import chain
import signal
import functools

import math
import pysam
from pybedtools import BedTool

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def check_position(pos, read_len, is_reverse, strand_check):
    '''Test if a position in a read is potentially informative of ancient
    DNA damage (C->T or G->A substitutions) given a library preparation
    method.
    '''
    if strand_check == 'USER':
        if not is_reverse and (pos == 0 or read_len - pos <= 2): return True
        if     is_reverse and (pos  < 2 or read_len - pos == 1): return True

    elif strand_check == 'USER_term5' and (pos < 5 or (read_len - pos <= 5)):
        return True

    elif strand_check == 'non-USER_term3' and (pos < 3 or (read_len - pos <= 3)):
        return True

    elif strand_check == 'non-USER_all':
        return True

    else:
        return False


def damage_at_site(pileup_info, strand_check=None):
    '''Ignore the base in the pileup of reads in case that:
       A) reference is C
          and:
           - USER treated: read has T on the first position or on the last two
           - non-USER treated:
               a) read has T on first three or last three positions
               b) read has T anywhere on the forward read
       B) reference is G
          and:
           - USER treated: read has A on the first two positions or on the last
           - non-USER treated:
               a) read has A at first three or last three positions
               b) read has A anywhere on the reverse read
    '''
    if not strand_check: return False

    ref_base, read_base, pos_in_read, read_len, reverse_strand, _ = pileup_info

    # if there is a C->T (on forward strand) or G->A (on reverse strand)
    # substitution at this site...
    if ((read_base == 'T' and not reverse_strand) or \
        (read_base == 'A' and     reverse_strand)):
        # ... check if it occured on a position likely to carry aDNA damage
        return check_position(pos_in_read, read_len, reverse_strand,
                              strand_check)
    return False

 
def filter_damage(pileup_column, strand_check):
    '''Filter out bases in a given pileup column that are likely result
    of DNA damage (C->T on forward strand, G->A on reverse strand)    '''
    return [pileup_info for pileup_info in pileup_column
                        if not damage_at_site(pileup_info, strand_check)]


def filter_bqual(pileup_column, minbq):
    '''Filter out bases in a given pileup column that are bellow a given
    base quality cut-off.
    '''
    return [pileup_info for pileup_info in pileup_column
                        if minbq <= pileup_info[5]]


def call_base(pileup_info, sampling_method):
    """Return the most frequently occuring element of a list."""
    bases = [base for _, base, _, _, _, _ in pileup_info]

    if sampling_method == 'majority':
        counts = Counter(bases).most_common()
        # take all bases with the highest count
        max_freq = max(c[1] for c in counts)
        bases = [c[0] for c in counts if c[1] == max_freq]
        # if there is more than one "best" choice, return nothing
        if len(bases) > 1:
            return None
        else: # return the majority allele
            return bases[0]
    elif sampling_method == 'consensus':
        return bases[0] if len(set(bases)) == 1 else None
    elif sampling_method == 'random':
        return random.choice(bases)


def bases_in_column(column, ref_base):
    '''Return a list of bases in a given pileup column.
    '''
    pileup = []

    # walk through all reads overlapping the current column
    # and accumulate bases at that position
    for pileup_read in column.pileups:
        # skip deletions
        if pileup_read.is_del: continue

        pos_in_read = pileup_read.query_position
        read_len = pileup_read.alignment.query_length
        read_base = pileup_read.alignment.query_sequence[pos_in_read]
        is_reverse = pileup_read.alignment.is_reverse
        baseq = ord(pileup_read.alignment.qual[pos_in_read]) - 33

        if read_base in "ACGT":
            pileup.append((ref_base,
                           read_base,
                           pos_in_read,
                           read_len,
                           is_reverse,
                           baseq))

    return pileup


def sample_bases(bam, ref, sampling_method, print_fn, minbq, mincov,
                 maxcov, strand_check=None, chrom=None, start=None, end=None):
    '''Sample bases in a given region of the genome based on the pileup
    of reads. If no coordinates were specified, sample from the whole BAM file.
    '''
    for col in bam.pileup(chrom, start, end):
        # if coordinates were specified, check first if a current column
        # lies within this region (pysam pileui return whole reads overlapping
        # a requested region, not just bases in this region)
        if chrom and start and end and not (start <= col.pos < end):
            continue
        else:
            ref_base = ref.fetch(col.reference_name, col.pos, col.pos + 1)

            pileup_bases = bases_in_column(col, ref_base)

            # filter out this position if it does not pass the coverage filter
            if not (mincov <= len(pileup_bases) <= maxcov): continue
            # filter out low qual bases or those that are most likely damage
            if strand_check: pileup_bases = filter_damage(pileup_bases, strand_check)
            # filter out low quality bases
            if minbq > 0: pileup_bases = filter_bqual(pileup_bases, minbq)

            # if there is any base in the pileup left, call one allele
            if len(pileup_bases) > 0:
                called_base = call_base(pileup_bases, sampling_method)
                if called_base:
                    print_fn(col.reference_name, col.pos + 1, ref_base,
                             called_base, len(pileup_bases))


def sample_in_regions(bam, bed, ref, sampling_method, print_fn,
                      minbq, mincov, maxcov, strand_check=None):
    '''Sample alleles from the BAM file at each position specified in a BED
    file. Return the result as a list of tuples in the form of
    (chromosome, position, ref_base, called_base).
    '''
    for region in bed:
        sample_bases(bam, ref, sampling_method, print_fn, minbq, mincov,
                     maxcov, strand_check, region.chrom, region.start, region.end)


def print_vcf_header(sample_name, handle):
    '''Print VCF header.'''
    print('##fileformat=VCFv4.1\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}'.
          format(sample=sample_name), file=handle)


def print_record(chrom, pos, ref, called, length, rec_fmt, handle):
    '''Print information about sampled site in a given string format.'''
    alt, gt = ('.', 0) if ref == called else (called, 1)
    print(rec_fmt.format(chrom=chrom, start=pos - 1, end=pos, pos=pos,
          ref=ref, allele=called, alt=alt, gt=gt, dp=length))


def main(argv=None):
    parser = argparse.ArgumentParser(description='Sample alleles from a given'
        ' BAM file based on pileup of reads, either by drawing random bases or'
        ' by performing a majority call at each position (from the whole BAM'
        ' or limited to regions specified by a BED file).')
    parser.add_argument('--bam', help='BAM file to sample from', required=True)
    parser.add_argument('--chrom', help='Chromosome to sample')
    parser.add_argument('--bed', help='BED file with coordinates of regions'
                        '/sites to sample')
    parser.add_argument('--ref', help='FASTA reference', required=True)
    parser.add_argument('--output', help='Name of the output file '
                        '(direct output to stdout if missing)', default=None)
    parser.add_argument('--format', help='Output as VCF or BED?',
                        choices=['VCF', 'BED'], required=True)
    parser.add_argument('--sample-name', help='Sample name to put in VCF')
    parser.add_argument('--method', help='How to sample alleles?',
                        choices=['majority', 'consensus', 'random'],
                        required=True)
    parser.add_argument('--strand-check', help='How and where to check for '
                        'damage? If not specified, no checks are performed.',
                        choices=['USER', 'USER_term5', 'non-USER_term3',
                                 'non-USER_all'], default=None)
    parser.add_argument('--minbq', help='Minimal quality of a base to be '
                        'considered for sampling (inclusive)', type=int,
                        default=0)
    parser.add_argument('--mincov', help='Required minimal coverage at '
                        'a position (inclusive)', type=int, default=1)
    parser.add_argument('--maxcov', help='Required maximal coverage at '
                        'a position (inclusive)', type=int, default=math.inf)

    # if there were no arguments supplied to the main function, use sys.argv
    # (skipping the first element, i.e. the name of this script)
    args = parser.parse_args(argv if argv else sys.argv[1:])

    if args.format == 'VCF' and not args.sample_name:
        parser.error('Sample has to be specified when outputting to VCF')

    bam = pysam.AlignmentFile(args.bam)
    if not bam.has_index():
        print("BAM file does not seem to index. An index is required for sampling.",
              file=sys.stderr)
        sys.exit(1)
    ref = pysam.FastaFile(args.ref)

    # output the results as specified by user
    handle = open(args.output, 'w') if args.output else sys.stdout

    if args.format == 'VCF':
        print_vcf_header(args.sample_name, handle)
        rec_fmt = '{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:DP\t{gt}:{dp}'
    else:
        rec_fmt = '{chrom}\t{start}\t{end}\t{ref}\t{allele}'

    print_fn = functools.partial(print_record, rec_fmt=rec_fmt, handle=handle)

    # if user specified a BED file, perform pileup on each region in that file
    if args.bed:
        bed = BedTool(args.bed)
        sample_in_regions(bam, bed, ref, args.method, print_fn,
                          args.minbq, args.mincov, args.maxcov,
                          args.strand_check)
    else: # otherwise scan the whole BAM file directly
        sample_bases(bam, ref, args.method, print_fn, args.minbq,
                     args.mincov, args.maxcov, args.strand_check,
                     args.chrom)

    if args.output:
        handle.close()


if __name__ == "__main__":
    main()
