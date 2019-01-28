# Introduction

ClinCNV does not perform low-level analysis and requires files to be prepared in specific formats. This guide will walk you through the data preparation.

We use [`ngs-bits`](https://github.com/imgag/ngs-bits) library for our preliminary preparation procedures which is fairly easy to install and fast. It can be replaced by standard tools such as `samtools`, however we don't provide examples with tools other than `ngs-bits`.


# Files you can have before CNV analysis

We assume that you have `.bam` files. For targeted sequencing you may also have `.bed` file provided by the manufacturer that describes design of your enrichment system.

## ClinCNV's input files you can get if you use targeted sequencing

**Important:** don't mix samples seuqenced with different enrichment kits - coverage patterns are different between them and ClinCNV is not supposed to remove such biases.

Depending on your goal, you may want to extract:

* on-target read coverage (has to be supported by on-target `.bed` file)

* off-target read coverage (has to be supported by off-target `.bed` file, generation of such file from on-target `.bed` file will be explained below)

* B-allele frequency files (for somatic framework)


## ClinCNV's input files you can get if you use whole genome sequencing

* on-target read coverage (has to be supported by on-target `.bed` file, generation of such file from coordinates of reference genome's chromosomes will be explained below)

* B-allele frequency files (for somatic framework)


## Calcuiation of read coverage