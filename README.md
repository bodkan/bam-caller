# bam-sample
`bam-sample` is a simple Python CLI tool for drawing random alleles from a given BAM file.
It generates output in a VCF format, as a simple tab-separated table or simply writes its
output to `stdout`.

It provides a set of options for processing ancient DNA data, such as filtering out
potential false SNPs caused by DNA damage, etc.

The program has several aDNA damage filtering switches that are confusing and don't make
a lot of sense, because it evolved while I was working on a very tricky ancient DNA dataset
which required a constant re-evaluation of different filtering strategies. It will be cleaned
up soon.
