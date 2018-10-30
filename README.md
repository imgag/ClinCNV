# ClinCNV

A tool for large-scale CNV and CNA detection.

Authors: G. Demidov, S. Ossowski.

## About this software

ClinCNV is supposed to detect CNVs in germline and somatic context (no mosaicism for now, but it can easily be implemented on request) in NGS data (targeted and whole-genome). We work in cohorts, so it makes sense to try ClinCNV if you have more than 10 samples (recommended amount - 40 since we estimate variances from the data). Note: by "cohort" we mean samples sequenced with the same enrichment kit with approximately the same depth (ie 1x WGS and 30x WGS better be analysed in a separate runs of ClinCNV). Of course it is better if your samples were sequenced within the same sequencing facility. Currently we work with hg19 only. For hg38 or mouse genome or any other diploid organism you have to replace *cytobands.txt* with the corresponding file. ClinCNV do not work with small panels (hundreds of regions) since GC-correction can not be performed accurately for samples sequenced with such panels.

## Pre-requisites

We expect you to install ClinCNV on Linux or MacOS platforms. We expect you to install R (as new version as possible, we used ClinCNV with R 3.2.3, but you may experience problems installing libraries using the old version) and the following libraries: 
```
install.packages("robustbase")
install.packages("MASS")
install.packages("data.table")
install.packages("foreach")
install.packages("doParallel")
```
For now we do not provide our own tool for pre-processing of the data. We recommend you to use *ngs-bits* (https://github.com/imgag/ngs-bits), however, as soon as your data match the format expected by ClinCNV you may proceed with any tool of choice (eg, *samtools*).

You should also have .bed file with the coordinates of targeted regions.

## Quick launch

You can try to start ClinCNV as follows:

```
	Rscript firstStep.R --normal normal.cov --out outputFolder --bed annotatedBedFile --folderWithScript $PWD
```

for *germline* samples and for *somatic* as

```
	Rscript firstStep.R --normal normal.cov --tumor tumor.cov  --out outputFolder --pair fileWithPairs --bed annotatedBedFile --folderWithScript $folderWithScript 
```

If it does not work, check if your files (.cov, .bed, file with pairs) are concordant with the descriptions below.

## File formats

Current version of ClinCNV works with 3 possible types of data: on-target reads, off-target reads, B-allele frequencies. For WGS obviously we work with only 2 of them (on-target and B-allele)

### .bed format
We expect .bed file annotated with GC-content and (optionally) intersecting genes. Header should be removed or commented with # symbol.
```
chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t gcContent[real, from 0 to 1] \t genesName[character comma delimited] \n
```

Example of .bed (here and below we provide only one line, assuming that there are as many as needed):

```
chr1    12171   12245   0.4595  DDX11L1
```

### .cov format
We expect AVERAGE coverage depths of samples to be written as (starting from header): 
```
chr \t start \t end \t sampleName1 \t sampleName2 \n
chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t averageCoverageDepth1[real] \t averageCoverageDepth1[real] \n
```
Example: 
```
chr   start   end     Sam1     Sam2
chr1    11166636        11166864        2374.32 1224.54
```

*Note1:* you may create such files for your samples separately and use the mergeFilesFromFolder.R script to merge them together.

*Note2:* if you suffer a lot with calculating average coverage, but you have the raw coverage depths, you can change the function
```
gc_and_sample_size_normalise <- function(info, coverages, averageCoverage=T, allowedChroms=NULL)
```
to 
```
gc_and_sample_size_normalise <- function(info, coverages, averageCoverage=F, allowedChroms=NULL)
```
in the file generalHelpers.R.

*Note3:* on-target and off-target reads should be pre-processed in .cov formats. If you do not have off-target reads for some samples - don't worry, ClinCNV will work with available data only.

*Note4:* Please be sure that you do not round your coverage of shallow-sequenced samples too much (e.g., the average coverage of the region is 0.0005, and you round it to 0.00).

### B-allele frequency format (expected file extension is .tsv)

Without header:
```
chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t uniqueID[char, can be used as chrom + coord merged] \t frequencyOfAltAllele[real from 0 to 1] \t coverageDepthOfThisPosition[int] \n
```
Example: 
```
chr1    2488153 2488153 chr1_2488153    0.4913  289
```
*Note1:* despite the fact we have start and end coordinates, B-allele frequency are expected to be calculated only from SNVs, not from indels.

*Note2:* if you do not have B-allele frequencies for some samples - don't worry, ClinCNV will work with available data only.


### File with information about pairs (normal vs tumor from the same sample)

We require presence of both normal and tumor input files to work in somatic context. To explain ClinCNV the connection between normals and tumors, you need to prepare file with the following format:

```
TumorSampleFromPatient1,NormalSampleFromPatient1
TumorSampleFromPatient2,NormalSampleFromPatient2
```

Please take care - sample names such as "TumorSampleFromPatient1" should match column name in .cov files and file name in .baf (if you want to use B-allele frequencies for this sample). The file can have any extension, we use "pairs.txt" to name such files.









## Hints and advices

### How to create .bed file for WGS

You will need a .bed file with start and end of each chromsome. For hg19 lengths of chromosomes can be found at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes , just add 0s as a second column.

Segmentation of whole genome with *ngs-bits*:

```
BedChunk -in hg19.bed -n $sizeOfBin -out "preparedBedHg19.bin"$sizeOfBin".bed"
```

where $sizeOfBin means pre-specified size of the segment (see below).

### How to create .bed file for off-target regions

Assume you have a .bed file $bedFile. This is how you create offtarget .bed:

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


### Off-target or WGS region size - how to choose?

To use ClinCNV in WGS and off-target contexts you need to choose a size of the window you want to segment your genome with.

From our experience, if you have shallow coverage (0.5x on average) - the window size should be 25kb at least, for 30x 1kb windows are totally OK and you can probably go to smaller window sizes. But if you have something intermediate, the criteria to choose the window size of off-target reads is: 1) not a lot of zero coverage regions (then the distributions will become zero-inflated and ClinCNV's results will be inaccurate), 2) approximate normality of coverage. To check the 2nd assumption, we recommend you to choose your desired window size, calculate coverage for ~30 samples, choose like 10-20 regions on random from the autosomes and built a density plot (plot(density(coverages)) in R). If you will see something that does not even remind you a bell shape, but has a large tail - you should increase the window size.

### How to annotate your .bed file with *ngs-bits*

```
BedAnnotateGC -in $bedFile -out "gcAnnotated."$bedFile
BedAnnotateGenes -in "gcAnnotated."$bedFile -out "annotated."$bedFile
rm "gcAnnotated.""$bedFile
```

### How to calculate coverage and form .cov file

With *ngs-bits*:
```
BedCoverage -bam $bamPath -in $bedPath -min_mapq 3 -out $sampleName".cov"
```


With *samtools*:

```
samtools bedcov $bedFilePath -Q 3 $bamPath > $sampleName".cov"
```

### How to calculate BAF-files

```
VariantAnnotateFrequency -in $nameOfNormalSample".vcf" -bam $nameOfSample".bam" -out $nameOfSample".tsv" -depth
 grep "^[^#]" $nameOfSample".tsv" | nawk -F'\t' '(length($4) == 1) && (length($5) == 1) {print $1 "\t" $2 "\t" $3 "\t" $1 "_" $2 "\t" $(NF-1) "\t" $NF}' > $nameOfSample".tsv"
```

*Not1:* ClinCNV uses BAF files only in somatic context now so you have to perform this procedure for both somatic and normal samples.

## Citation

ClinCNV is not published for now so it is not possible to properly cite the paper. However, you can: 

1) cite it as an unpublished tool in the text of your paper;

2) ask us to help you with the analysis for the co-authorship.

Paper is coming soon, stay tuned!