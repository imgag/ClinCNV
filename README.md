# ClinCNV

A tool for large-scale CNV and CNA detection.

Authors: G. Demidov, S. Ossowski.

## About this software

ClinCNV is supposed to detect CNVs in germline and somatic context (no mosaicism for now, but it can easily be implemented on request) in NGS data (targeted and whole-genome). We work in cohorts, so it makes sense to try ClinCNV if you have more than 10 samples (recommended amount - 40 since we estimate variances from the data). Note: by "cohort" we mean samples sequenced with the same enrichment kit with approximately the same depth (ie 1x WGS and 30x WGS better be analysed in a separate runs of ClinCNV). Of course it is better if your samples were sequenced within the same sequencing facility. Currently we work with hg19 only. For hg38 or mouse genome or any other diploid organism you have to replace *cytobands.txt* with the corresponding file.

## Pre-requisites

We expect you to install ClinCNV on Linux or MacOS platforms. We expect you to install R (as new version as possible, we used ClinCNV with R 3.2.3, but you may experience problems installing libraries using the old version) and the following libraries: 

install.packages("robustbase")
install.packages("MASS")
install.packages("data.table")
install.packages("foreach")
install.packages("doParallel")

For now we do not provide our own tool for pre-processing of the data. We recommend you to use *ngs-bits* (https://github.com/imgag/ngs-bits), however, as soon as your data match the format expected by ClinCNV you may proceed with any tool of choice (eg, samtools).

You should also have .bed file with the coordinates of targeted regions.

## File formats

Current version of ClinCNV works with 3 possible types of data: on-target reads, off-target reads, B-allele frequencies. For WGS obviously we work with only 2 of them (on-target and B-allele)

### .bed format
We expect .bed file annotated with GC-content and (optionally) intersecting genes. Header should be removed or commented with # symbol.

chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t gcContent[real, from 0 to 1] \t genesName[character comma delimited] \n

### .cov format
We expect AVERAGE coverage depths of samples to be written as (starting from header): 

chr \t start \t end \t sampleName1 \t sampleName2 \n

chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t averageCoverageDepth1[real] \t averageCoverageDepth1[real] \n

Example: 

chr   start   end     Sam1     Sam2

chr1    11166636        11166864        2374.32 1224.54

Please be sure that you do not round your coverage of shallow-sequenced samples too much (e.g., the average coverage of the region is 0.0005, and you round it to 0.00).

*Note1:* you may create such files for your samples separately and use the mergeFilesFromFolder.R script to merge them together.

*Note2:* if you suffer a lot with calculating average coverage, but you have the raw coverage depths, you can change the function

gc_and_sample_size_normalise <- function(info, coverages, averageCoverage=T, allowedChroms=NULL)

to 

gc_and_sample_size_normalise <- function(info, coverages, averageCoverage=F, allowedChroms=NULL)

in the file generalHelpers.R.

*Note3:* on-target and off-target reads should be pre-processed in .cov formats. If you do not have off-target reads for some samples - don't worry, ClinCNV will work with available data only.

### B-allele frequency format (expected file extension is .tsv)

Without header:

chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t uniqueID[char, can be used as chrom + coord merged] \t frequencyOfAltAllele[real from 0 to 1] \t coverageDepthOfThisPosition[int] \n

Example: 

chr1    2488153 2488153 chr1_2488153    0.4913  289

*Note1:* despite the fact we have start and end coordinates, B-allele frequency are expected to be calculated only from SNVs, not from indels.

*Note2:* if you do not have B-allele frequencies for some samples - don't worry, ClinCNV will work with available data only.



## Hints and advices

### Off-target or WGS region size - how to choose?

To use ClinCNV in WGS and off-target contexts you need to choose a size of the window you want to segment your genome with.

From our experience, if you have shallow coverage (0.5x on average) - the window size should be 25kb at least, for 30x 1kb windows are totally OK and you can probably go to smaller window sizes. But if you have something intermediate, the criteria to choose the window size of off-target reads is: 1) not a lot of zero coverage regions (then the distributions will become zero-inflated and ClinCNV's results will be inaccurate), 2) approximate normality of coverage. To check the 2nd assumption, we recommend you to choose your desired window size, calculate coverage for ~30 samples, choose like 10-20 regions on random from the autosomes and built a density plot (plot(density(coverages)) in R). If you will see something that does not even remind you a bell shape, but has a large tail - you should increase the window size.

### How to annotate your .bed file with *ngs-bits*

BedAnnotateGC -in $bedFile -out "gcAnnotated."$bedFile
BedAnnotateGenes -in "gcAnnotated."$bedFile -out "annotated."$bedFile
rm "gcAnnotated.""$bedFile




## Citation

ClinCNV is not published for now so it is not possible to properly cite the paper. However, you can: 

1) cite it as an unpublished tool in the text of your paper;

2) ask us to help you with the analysis for the co-authorship.

Paper is coming soon, stay tuned!