# How to run `ClinCNV` for germline analysis

Germline CNVs detection pipeline of `ClinCNV` may utilize on-target reads coverage (reads that were aligned to the pre-specified and enriched locations in the genome or, in case of whole genome sequencing (WGS), reads that were mapped somewhere in the genome) and off-target reads coverage (even for the perfectly prepared enrichment there is typically a large number of reads that come from regions outside of your pre-specified intervals) as evidence of genomic ploidy. Off-targed reads coverage is usually summarised in windows of large size, such as 50KBps or even bigger, since the coverage outside of enriched areas is extremely low and sufficient number of reads is required for usage of number of mapped reads as a marker of underying copy-number.

We suggest the following step by step guide into CNVs detection in different use cases (data after enrichment and WGS, idealy with the smallest amount of PCR cycles as possible):

## Tagret sequencing

1. Simple analysis, having only on-target coverage from the cohort of samples (generation of such files described in the documentation):
`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed`

2. Specifying output folder **(here and until the end of the numbered list instruction you should add the line to the line you got on the previous step)**:
`--out /your/folder`. Folder `normal` will be created in your output folders, results will be put into subfolders, one sample per subfolder.

3. Adding offtarget coverages (does not increase _Sensitivity_ a lot since germline CNVs are rarely as long as off-target window, but improves _Specificity_ efficiently removing CNVs formed by probes standing far away from each other, off-target coverage has to be extracted for *sufficiently large number of samples*, but not necessarily for all the samples from `normal.cov` file):
`--normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed`

4. Playing with _Sensitivity_ and _Specificity_ balance (increase of the threshold `--scoreG` leads to higher _Specificity_ and lower _Sensitivity_ and vice versa, default value is 20):
`--scoreG 60`

5. Increasing or decreasing minimum length of detected variant (default = 2 data points or bigger, the real number of data points is actually +1 to you specify here, 0 means at least 1 data point):
`--lengthG 0`

6. Increasing the number of CNVs you'd expect to get from that sample (default = 10000, but for WGS samples better make it 1000 or even bigger, when such threshold is exceeded, quality threshold is increased and sample is re-analysed as many times as specified with `--maxNumIter` flag):
`--maxNumGermCNVs 2000 --maxNumIter 5`


7. Dividing your large cohort into several clusters of similar samples for more accurate parameter estimation and speeding the tool up. You specify the minimum size of the cluster (we recommend to keep it in between 20 and 100):
` --minimumNumOfElemsInCluster 50`

8. Speeding up the tool by using several cores (only some parts of the pipeline are parallelised):
`--numberOfThreads 4`

In the end you will have a command line like:

`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 0 --maxNumGermCNVs 2000 --maxNumIter 5 --minimumNumOfElemsInCluster 50 --numberOfThreads 4`

## WGS

For high-coverage (>10X) WGS you should pay special attention to step 6 (may be put the maximum number of expected CNVs equal to 100.000 if you run your dataset for the first time and don't know what to expect). 

For low-coverage WGS (<10X) we would recommend to use window size at least 5KB and minimal size of CNV as 3.

## Polymorphic CNVs 

CNVs whose frequency is larger than 2.5% are usually more difficult to detect using conventional methods. Thus, we detect them using Gaussian Mixture Models. This is a separate routine within the same script. To call such variants, you need to specify

`--polymorphicCalling YES`

If instead of YES you provide the path to the `.bed` file with the coordinates of polymorphic regions, `ClinCNV` will ignore them during the calling.

Otherwise you may get such pictures:

![Polymorphic region][polymorph_call]

Here horizontal lines show different copy-numbers and each dot denotes one sample from the cohort.




## Mosaic CNVs calling

Add `--mosaicism` option to your command line. Mosaic CNVs can take values from 1.1 to 2.9 with the step of 0.05.

Several times this mode was used for CNAs calling in tumor samples under the absence of matching normal. It can be done only if CNAs affect less than a half of your sample and still will be less accurate and much less informative than the paired calling. To do the same trick, you need to manually chage `germlineSolver.R` file - there is a line starting with  `cn_states_mosaicism <- seq(from=`. Modify it accordingly (from should be the minimal copy-number you want to detect, to is the maximum, step is the step of discrete grid). Same can be done for looking for aneuploidies in prenatal testing settings - we've never tested that. You also may change the value `fineForMosaicism = 0.05` to something smaller. This fine downweights all the mosaic CNVs and thus the caller prefers integer copy-numbers, but if you want something different - go ahead. 

## Analysis of germline CNVs in Whole Exome Sequencing using Agilent SureSelect v6 and v7 panels and read depth 100-140

We have validated a script for the annotation of CNVs, namely - how probable is that this particular CNV is a False Positive. It is located in `https://github.com/imgag/ClinCNV/tree/master/helper_scripts/FDR_annotation`. An example run is:

`Rscript germlineCallsAnnotatorFDR.R --pathToClassifier /helper_scripts/FDR_annotation/randomForests.RObj --input your_clinCNV_file_with_CNVs.tsv --output annotated_file_with_cnvs.tsv`

It will add a column with FDR. FDR of 0.01 means that this CNV has only 1% chance to be a False Positive. 0.9 - it has 90% chance to be False Positive, however, if you see that this CNV is affecting a gene that may be responsible for the disease, that means this CNV deserves a second look, even if it has high chances to be a False Discovery.

All CNVs are evaluated as they are singletons (happen once per studies cohort) - it is possible to achieve higher power for recurrent CNVs, however, usually pathogenic CNVs are rare.
 


# How to interpret results

Your output files will be saved in `/your/folder/normal` directory - `ClinCNV` will make a separate folder per sample and put resulting files there. 

There (if you did not specify flag `--visulizationIGV` as `F`) you can find 2 [IGV tracks](http://software.broadinstitute.org/software/igv/SEG) - one for real coverage values, second for CNVs track. Both tracks heights denote the underlying copy number (height = 0 denote homozygous deletion, height = 1 denote heterozygous deletion, height = 2 = diploid, height = 3 heterozyous duplication, etc.). Height is cut at level 6 since othervise differences would be much less illustrative (so, even if you have values with copy-number >6, they will be depicted just by points at height 6). You can just drag and drop these tracks into your [IGV browser](http://software.broadinstitute.org/software/igv/home).

**Coverage values** are GC-normalised and median-normalised (that means that "average" coverage for cohort is centered around copy-number 2 in autosomes). But coverage is not square root normalised as it is inside `ClinCNV`'s algorithm - so sometimes even if it seems that the dot is closer to alternative copy-number it may not be true since variances are different for different copy numbers.

**CNV-segments** here are mainly used for algorithmic diagnostic purposes - if the segments are too fragmented or if you clearly see some CNVs in non-polymorphic (variance is usually overestimated for regions that are affected by CNVs frequently) regions that were not detected by ClinCNV, then you've chosen quality threshold as too small/big, respectively, and you need to increase it and re-do the analysis. CNVs are described much more informatively in files "sampleName_cnvs.tsv"

![IGV tracks for one sample (exome seq)][IGV_track]



You can scroll into your data at chromosome level:


![IGV tracks for one sample (exome seq)][IGV_track_chr]

Or deeper at a single event level:

![IGV tracks for one sample (exome seq)][IGV_track_cnv]

## "sampleName_cnvs.tsv" files

You will get tab delimited files with 2 header lines (starting with hashtag), results' table column names and the results table itself.

![Table with results][table_of_results]

The first 3 lines are technical. 4th line shows the inferred gender, 5th line - how many interations of quality score increasing passed since the number of CNVs in the sample became acceptable, quality used is the quality score used, was outlier in clustering? - shows if this sample was an outlier after batch effect clustering.

The last line -- fraction of outliers -- shows how many coverage data points had p-value less than 0.05. It can not really be used for the QC control since having large aneuploidy will lead to high fraction of outliers.

### Columns in the table

`#chr	start	end` - are self explainable - they contain coordinates of detected variants.

`CN_change` denotes a copy number of detected variant. 

`loglikelihood` actually says "how more probable is the detected copy number of the particular region comparing to the baseline (diploid or 1 copy for sex chromosomes in males)?". You can find more detailed description on [wikipedia](https://en.wikipedia.org/wiki/Bayes_factor#Interpretation) (use the second table). According to the table, loglikelihood > 10 means that the strength of evidence in favor of alternative copy number is "Very strong", however, we would recommend you to keep the values even higher (40 or at least 20) - genomic coverage is affected by numerous events such as short indels causing alignment problems, batch effects, sequencing depth, technical artifacts, etc., and the genome is quite long so seeing large loglikelihood changes is not so rare there.

`no_of_regions` shows how many datapoints are included into the variant. Longer variants usually are more credible (but not always, e.g., long variants with small `loglikelihood` value are likely to be false positives).

`length_KB` is just difference between `end` and `start` of the variant.

`potential_AF` shows, how often coverage in particular regions is unusually lower/higher (in case when the variant is deletion/duplication respectively).

`genes` is filled with information only if you used annotated `.bed` file. We recommend to annotate only on-target `.bed` files from the panel since you are mainly interested in ploidy of genes you've included into your panel. If you don't annotate your `.bed` file - you will see just zeros in this column.

`qvalue` is basically p-value from Z-test, corrected for multiple test correction. High p-values are expected to be in polymorphic short regions.


[polymorph_call]: ./images/polymorphic.png "Polymorphic CNV region"
[IGV_track]: https://github.com/imgag/ClinCNV/raw/master/doc/images/germline_tracks.png "IGV tracks for germline sample"
[IGV_track_chr]: https://github.com/imgag/ClinCNV/raw/master/doc/images/germline_tracks_chrom_level.png "IGV tracks for germline sample (chromosome level)"
[IGV_track_cnv]: https://github.com/imgag/ClinCNV/raw/master/doc/images/germline_tracks_cnv_level.png "IGV tracks for germline sample (one CNV level)"
[table_of_results]: ./images/normal_calling.png "Table with results (simulated data)"