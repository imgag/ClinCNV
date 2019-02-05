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

* text files with comma separated sample IDs where each line contains IDs of mother, father, kid _(order is important)_ for **trio** framework

* text files with comma separated sample IDs where each line contains IDs of tumor and normal samples _(order is important)_ for **somatic** framework


## ClinCNV's input files you can get if you use whole genome sequencing (below "WGS")

* on-target read coverage (has to be supported by on-target `.bed` file, generation of such file from coordinates of reference genome's chromosomes will be explained below)

* B-allele frequency files (for somatic framework)


## Generation of `.bed` files (on- and off-target)

### Targeted sequencing

You should already have on-target `.bed` file. To generate off-target `.bed` file from your on-target file `panel.bed` you need to run following commands (making offset of 400 bp from the sides of targeted regions, window sike of 50KB and removing all regions that are less than 25KB in the end):

```
# Determine offtarget with offset of 400 to the left and to the right of the targeted region
BedExtend -in $bedFile -n 400 | BedMerge -out "extended_"$bedFile
# hg19.bed here is from the previous paragraph
BedSubtract -in hg19.bed -in2 "extended_"$bedFile -out offtarget.bed
# Chunk offtarget into pieces of 50kbps
BedChunk -in offtarget.bed -n 50000 -out offtarget_chunks.bed
# Remove regions <25k
BedShrink -in offtarget_chunks.bed -n 12500 | BedExtend -n 12500 -out "offtarget_chunks_"$bedFile
```
### .bed file annotation

`ngsbits` has to be installed. 

```
BedAnnotateGC -in $bedFile -out "gcAnnotated."$bedFile -ref reference.fa
BedAnnotateGenes -in "gcAnnotated."$bedFile -out "annotated."$bedFile
rm "gcAnnotated.""$bedFile
```


### WGS

For WGS you need to prepare file with chromosomes' starts and ends based on reference genome you've used. Then you can use the same strategy as for off-target `.bed` file generation:

```
BedChunk -in startAndEndOfChromosomes.bed -n 1000 -out chunks.bed
```

You may not care about centromeric regions - they will be excluded by `ClinCNV`. 1000bp is OK for 30x and more coverage (for 40x you can go for 500bp - if your biologists did a good job and the level of noise is low) which means CNV length of 3KB-1.5KB. But for shallow WGS you should choose 5000 or 10000 bp (0.1x would require 10KB, 5x may go up to 5KB).

A rule of thumb for CNV detection resolution with `ClinCNV` - the smallest CNV you may try to detect should contain 3 windows, mainly because the windows on the left and on the right can be partially affected by CNVs so you need to have at least one "central" window to be able to genotype variant. (Imagine: you have a copy number 4 variant that affect only 2 neighboring windows and only partially - you may see a copy number 3 in the end, depicted below)

```
Variant:
__________^^^^^^^^^^___________
2222222222444444444422222222222

What coverage we see in our windows:
_______|_______|_______|_______
   2       3       3       2   
```

## Calculation of read coverage (both on- and off- target)

Let's say you want to calculate read coverage, using `.bed` file with path specified with `bed_file`, `.bam` file specified with `BAM` variable and sample name `output`. This command has to be executed:

```BedCoverage -bam $BAM -in $bed_file -min_mapq 5 -decimals 4 > $output".cov"```

Please keep files that you obtain from normals/tumors or on-target/off-target or of course different sequencing kits in separate folders - you will need to merge this files later!

Then you need to merge your ".cov" files into one table. To do this, you can use script `mergeFilesFromFolder.R` script provided with `ClinCNV` using `input_folder` and `output_folder` as variables to keep your absolute paths:

```Rscript mergeFilesFromFolder.R -i $input_folder -o $output_folder```



## B-allele frequency for Somatic framework

B-allele frequency extraction works in 2 steps - first, you extract high quality positions from your `.vcf` file from the normal sample, then you calculate coverage of this position in your tumor sample using BAF file obtained at previous step.

You should perform standard germline variant calling on you normal sample and get `.vcf` as output. `.vcf` files may differ - so you may need to slightly change the script.


## Text files with sample IDs

Due to variability of data storage formats we can not suggest you a single recipy on how to get files with sample IDs. But in the end the files you should get need to look like:

**Trios:**

```
Mother_ID1,Fater_ID1,Kid_ID1
Mother_ID2,Fater_ID2,Kid_ID2
etc
```

**Normal-Tumor pairs:**
```
Tumor_ID1,Normal_ID1
Tumor_ID1,Normal_ID2
etc
```

"Mother_ID1" is just an example, it can be anything (better without special characters such as `#`). If you have sequenced multiple tumor probes and only one normal per one sample, specify it as it is. Main condition is: *those IDs has to be presented in the file with coverages with exactly same names.*



