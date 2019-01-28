# Introduction

ClinCNV does not perform low-level analysis and requires files to be prepared in specific formats. This guide will walk you through the data preparation.

We use [`ngs-bits`](https://github.com/imgag/ngs-bits) library for our preliminary preparation procedures which is fairly easy to install and fast. It can be replaced by standard tools such as `samtools`, however we don't provide examples with tools other than `ngs-bits`.


# Files you should have before CNV analysis

We assume that you have `.bam` files. For targeted sequencing you may also have `.bed` file provided by the manufacturer that describes design of your enrichment system. We also assume that you know which reference was used for `.bam` files generation. Files for __hg19__ are already provided, files for other genomes such as __hg38__ or __mm10__ may have to be prepared by yourself. 

## ClinCNV's input files you can get if you use targeted sequencing

**Important:** don't mix samples seuqenced with different enrichment kits - coverage patterns are different between them and ClinCNV is not supposed to remove such biases.

Depending on your goal, you may want to extract:

* on-target read coverage (has to be supported by on-target `.bed` file)

* off-target read coverage (has to be supported by off-target `.bed` file, generation of such file from on-target `.bed` file will be explained below)

* B-allele frequency files (for **somatic** framework)

* text files with comma separated sample IDs where each line contains IDs of mother, father, kid __(order is important)__ for **trio** framework

* text files with comma separated sample IDs where each line contains IDs of tumor and normal samples __(order is important)__for **somatic** framework


## ClinCNV's input files you can get if you use whole genome sequencing (below "WGS")

* on-target read coverage (has to be supported by on-target `.bed` file, generation of such file from coordinates of reference genome's chromosomes will be explained below)

* B-allele frequency files (for somatic framework)


## Generation of `.bed` files (on- and off-target)

### Targeted sequencing

You should already have on-target `.bed` file. To generate off-target `.bed` file from your on-target file `panel.bed` you need to run following commands:

### WGS

For WGS you need to prepare file with chromosomes' starts and ends based on reference genome you've used. 

## Calcuiation of read coverage (both on- and off- target)


## B-allele frequency for Somatic framework

## Text files with sample IDs

Due to variability of data storage formats we can not suggest you a single recipy on how to get files with sample IDs. 



