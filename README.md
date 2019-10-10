# ClinCNV
[![Build Status](https://travis-ci.org/imgag/ClinCNV.svg?branch=master)](https://travis-ci.org/imgag/ClinCNV)

A tool for large-scale CNV and CNA detection.

Authors: G. Demidov, S. Ossowski.

Any issues should be reported to: *german dot demidov at medizin dot uni-tuebingen dot de *.

This software is distributed under [MIT licence](./LICENSE).

## About this software

ClinCNV is supposed to detect CNVs in germline and somatic context (no mosaicism for now, but it can easily be implemented on request) in NGS data (targeted and whole-genome). We work in cohorts, so it makes sense to try ClinCNV if you have more than 10 samples (recommended amount - 40 since we estimate variances from the data). By "cohort" we mean samples sequenced with the same enrichment kit with approximately the same depth (ie 1x WGS and 30x WGS better be analysed in separate runs of ClinCNV). Of course it is better if your samples were sequenced within the same sequencing facility. Currently we work with hg19 only. For hg38 or mouse genome or any other diploid organism you have to replace *cytobands.txt* with the corresponding file. ClinCNV do not work with small panels (hundreds of regions) since GC-correction can not be performed accurately for samples sequenced with such panels.

NOTE: Folder `PCAWG` was used for CNVs detection in PanCancer Analysis of Whole Genomes cohort and is *research* only version. It is located here for historical reasons. Feel free to remove it.

## Pre-requisites

We expect you to install ClinCNV on Linux or MacOS platforms. Copy all the files using `git clone https://github.com/imgag/ClinCNV.git`. We expect you to install `R` (as new version as possible, we used ClinCNV with `R 3.2.3`, but you may experience problems installing libraries using the old version) and the following libraries: 

```
install.packages("optparse")
install.packages("robustbase")
install.packages("MASS")
install.packages("data.table")
install.packages("foreach")
install.packages("doParallel")
install.packages("mclust")
install.packages("R.utils")
install.packages("RColorBrewer")
install.packages("party")
```

**Test run:**
```
fold=/folder/with/the/cloned/ClinCNV
mkdir $fold"/results"
Rscript $fold"/clinCNV.R" --bed $fold"/samples/bed_file.bed" --normal $fold"/samples/coverages_normal.cov" --out $fold"/results"
```

If `ClinCNV` fails with the test run, set `chmod 755` to the output folder.

You will find result in `$fold/ClinCNV/result/` folder.

For now we do not provide our own tool for pre-processing of the data. We recommend you to use *ngs-bits* (https://github.com/imgag/ngs-bits), however, as soon as your data match the format expected by ClinCNV you may proceed with any tool of choice (eg, *samtools*).

You should also have `.bed` file with the coordinates of targeted regions and reference genome in `.fasta` format for annotation of `.bed` file with *ngs-bits*.

## Quick launch

More informative manuals are located in the `doc` folder.

You can try to start ClinCNV as follows:

```
Rscript firstStep.R --normal normal.cov --out outputFolder --bed annotatedBedFile --folderWithScript $PWD
```

for *germline* samples and for *somatic* as

```
Rscript firstStep.R --normal normal.cov --tumor tumor.cov  --out outputFolder --pair fileWithPairs --bed annotatedBedFile --folderWithScript $PWD 
```

If it does not work, check if your files (.cov, .bed, file with pairs) are concordant with the descriptions below.

## Use cases

ClinCNV now can work in 3 different contexts: germline calling, somatic calling (normal / tumor pairs) and trios. To run ClinCNV in germline context, you should minimally specify `--bed` and `--normal`. For trios, parameter `--triosFile` has to be specified (sample names for a child, his/hers mother and father, divided by comma). For somatic `--tumor` and `--pairs` has to be specified (pairs should contain tumor and normal sample names, divided by comma). It is highly recommended to use B-allele frequencies in somatic context. A folder with `.tsv` files has to be specified to switch to B-allele frequencies mode. Files need to have same sample names as column names in `--normal` and `--tumor` files. If ClinCNV does not find B-allele frequencies for a particular sample, it tries to detect CNVs using read depths only.

## File formats

Current version of ClinCNV works with 3 possible types of data: on-target reads, off-target reads, B-allele frequencies. For WGS obviously we work with only 2 of them (on-target and B-allele). For shallow WGS only "on-target" coverage is informative.

To perform operations with these types of data, we specify several formats of files: annotated `.bed` file, `.cov` file, `.tsv` file with information about BAF.

### Annotated .bed format
We expect `.bed` file annotated with GC-content and (optionally) intersecting genes. Header should be removed or commented with # symbol.
```
chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t gcContent[real, from 0 to 1] \t genesName[character comma delimited] \n
```

Example of `.bed` (here and below we provide only one line, assuming that there are as many as needed):

```
chr1    12171   12245   0.4595  DDX11L1
```

### .cov format
We expect AVERAGE coverage depths of samples to be written as (starting from header): 
```
chr \t start \t end \t sampleName1 \t sampleName2 \n
chrI[char, "chr" is a prefix] \t startCoord[int] \t endCoord[int] \t averageCoverageDepth1[real] \t averageCoverageDepth2[real] \n
```
Example: 
```
chr   start   end     Sam1     Sam2
chr1    11166636        11166864        2374.32 1224.54
```

*Note1:* you may create such files for your samples separately and use the `mergeFilesFromFolder.R` script to merge them together.

*Note2:* if you suffer a lot with calculating average coverage, but you have the raw coverage depths, you can change the function
```
gc_and_sample_size_normalise <- function(info, coverages, averageCoverage=T, allowedChroms=NULL)
```
to 
```
gc_and_sample_size_normalise <- function(info, coverages, averageCoverage=F, allowedChroms=NULL)
```
in the file `generalHelpers.R`.

*Note3:* on-target and off-target reads should be pre-processed in `.cov` formats. If you do not have off-target reads for some samples - don't worry, ClinCNV will work with available data only.

*Note4:* Please be sure that you do not round your coverage of shallow-sequenced samples too much (e.g., the average coverage of the region is 0.0005, and you round it to 0.00).

*Note5:* Names of columns (sample names) are meaningful and should match between `normal.cov`, `tumor.cov`, `pairs.txt` files.

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

Please take care - sample names such as "TumorSampleFromPatient1" should match column name in `.cov` files and file name in BAF `TumorSampleFromPatient1.tsv` files (if you want to use B-allele frequencies for this sample). The file can have any extension, we use "pairs.txt" to name such files.









## Hints and advices

### How to create .bed file for WGS

You will need a `.bed` file with start and end of each chromsome. For hg19 lengths of chromosomes can be found at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes , just add 0s as a second column.

Segmentation of whole genome with *ngs-bits*:

```
BedChunk -in hg19.bed -n $sizeOfBin -out "preparedBedHg19.bin"$sizeOfBin".bed"
```

where $sizeOfBin means pre-specified size of the segment (see below).

### How to create .bed file for off-target regions

Assume you have a `.bed` file $bedFile. This is how you create offtarget .bed:

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

From our experience, if you have shallow coverage (0.5x on average) - the window size should be 25kb at least, for 30x 1kb windows are totally OK and you can probably go to smaller window sizes. But if you have something intermediate, the criteria to choose the window size of off-target reads is: 1) not a lot of zero coverage regions (then the distributions will become zero-inflated and ClinCNV's results will be inaccurate), 2) approximate normality of coverage. To check the 2nd assumption, we recommend you to choose your desired window size, calculate coverage for ~30 samples, choose like 10-20 regions on random from the autosomes and built a density plot (`plot(density(coverages)`) in `R`). If you will see something that does not even remind you a bell shape, but has a large tail - you should increase the window size.

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

Using *ngs-bits*:

```
VariantAnnotateFrequency -in $nameOfNormalSample".vcf" -bam $nameOfSample".bam" -out $nameOfSample".tsv" -depth
 grep "^[^#]" $nameOfSample".tsv" | nawk -F'\t' '(length($4) == 1) && (length($5) == 1) {print $1 "\t" $2 "\t" $3 "\t" $1 "_" $2 "\t" $(NF-1) "\t" $NF}' > $nameOfSample".tsv"
```

*Note1:* ClinCNV uses BAF files only in somatic context now so you have to perform this procedure for both somatic and normal samples.

## Full list of parameters of ClinCNV

`--normal` - file with normal samples' coverages depths in `.cov` format. On-target or WGS coverage.

`--tumor` - file with tumor samples' coverages depths in `.cov` format. On-target or WGS coverage.

`--normalOfftarget` - file with normal samples' coverages depths in `.cov` format. Off-target coverage.

`--tumorOfftarget` - file with tumor samples' coverages depths in `.cov` format. Off-target coverage.

`--out` - output folder, default: `./result/`

`--pair` - comma separated table. Each row = "tumorSampleName,normalSampleName"

`--bed` - `.bed` file with coordinates of on-target regions, GC-annotated. Gene anotations in the last column can be used in output.

`--bedOfftarget` - `.bed` file with coordinates of off-target regions, GC-annotated. Gene anotations in the last column can be used in output.

`--colNum` - number of column in `.cov` file where coverage depth columns start. Normally equal to 4.

`--folderWithScript` - obsolete. Basically it was $PWD

`--reanalyseCohort` - if equal to "T", all the samples are analysed. If equal to "F", only samples which subfolders do not exist in output folder will be analysed

`--scoreG` - germline CNV score threshold (loglikelihood difference). For additional information: https://en.wikipedia.org/wiki/Bayes_factor#Interpretation . As a guideline: 30 is usually a sensitive, but not specific option, 50 is balanced, 100 is quite strict and leads to high specificity, but (possibly) low sensitivity. 

`--lengthG` - minimum length of germline CNV (number of markers = on- and off-target regions). Depends on your desired purposes. Recomended to keep up not smaller than 2.

`--scoreS` - somatic score threshold. Since CNAs (copy-number changes in caner) are usually long and quite significant, but FFPE introduce huge noise into the tumor sequencing data, we do not recommend to keep it lower than 50. As a rule of thumb, we use threshold of 100 for panel of ~500 cancer genes and 200 for WES samples.

`--lengthS` - minimum length of somatic CNA (# of markers, for WES it can be ~5-10 since CNAs are usually long).

`--maxNumGermCNVs` - maximum number of allowed germline CNVs. If sample has too many variants, thresholds will be automatically increased and the sample will be reanalysed several times.

`--maxNumSomCNAs` - maximum number of allowed somatic CNAs.

`--maxNumIter` - maximum number of iterations "increase threshold - detect CNVs - check if the number of detected variants less than specified above numbers"

`--bafFolder` - folder with `.tsv` files containing information about B-allele frequencies in tumor and normal samples. Note: both tumor and normal sample have to be presented if you want to use it as a predictor.

`--normalSample` - name of the sample (if only one germline sample is expected to be calculated). All the parameters are estimated for the whole cohort. Has to be presented in a header of the file with normal samples. Otherwise, the tool will fail with assert message.

`--tumorSample` - name of the tumor sample (if only one sample is expected to be calculated). Has to be presented in 1) file pairs.txt in pair with the sample specified by --normalSample option, 2) in a header of the file with tumor samples. Otherwise, the tool will fail with assert message.

`--fdrGermline` - number of iteration of FDR control procedure performed (more is better for FDR control, but slower). If 0, no permutations and FDR control performed (default).

`--numberOfThreads` - maximum number of threads allowed to use by the tool.


## Citation

ClinCNV is not published for now so it is not possible to properly cite the paper. However, you can: 

1) cite it as an unpublished tool in the text of your paper;

2) ask us to help you with the analysis for the co-authorship.

Paper is coming soon, stay tuned!
