# ClinCNV
[![Build Status](https://travis-ci.org/imgag/ClinCNV.svg?branch=master)](https://travis-ci.org/imgag/ClinCNV)

A tool for large-scale CNV and CNA detection.

Authors: G. Demidov, S. Ossowski.

Any issues should be reported to: *german dot demidov at medizin dot uni-tuebingen dot de *.

This software is distributed [under MIT licence](./LICENSE).

## About this software

ClinCNV is supposed to detect CNVs in germline and somatic context (no mosaicism for now, but it can easily be implemented on request) in NGS data (targeted and whole-genome). We work in cohorts, so it makes sense to try ClinCNV if you have more than 10 samples (recommended amount - 40 since we estimate variances from the data). By "cohort" we mean samples sequenced with the same enrichment kit with approximately the same depth (ie 1x WGS and 30x WGS better be analysed in separate runs of ClinCNV). Of course it is better if your samples were sequenced within the same sequencing facility. Currently we work with hg19 only. For hg38 or mouse genome or any other diploid organism you have to replace *cytobands.txt* with the corresponding file. ClinCNV do not work with small panels (hundreds of regions) since GC-correction can not be performed accurately for samples sequenced with such panels.

NOTE: version of ClinCNV we are talking about is located in the folder `somatic`. Folder `PCAWG` was used for CNVs detection in PanCancer Analysis of Whole Genomes cohort and is *research* only version. It can be useful if you want to analyse several thousands of whole genomes, but better contact us in this case.

## Install and run ClinCNV
To run ClinCNV please refer to it's [user manual](./somatic/README.md).
