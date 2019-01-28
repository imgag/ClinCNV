# How to run `ClinCNV` for germline analysis

## Tagret sequencing

1. Simple analysis, having only on-target coverage:
`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed`

2. Specifying output folder **(here and until the end of the numbered list instruction you should add the line to the line you got on the previous step)**:
`--out /your/folder`

3. Adding offtarget coverages (does not increase _Sensitivity_ a lot since germline CNVs are rarely as long as off-target window, but improves _Specificity_ efficiently removing CNVs formed by probes standing far away from each other):
`--normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed`

4. Playing with _Sensitivity_ and _Specificity_ balance (increase of the threshold `--scoreG` leads to higher _Specificity_ and lower _Sensitivity_ and vice versa, default value is 40):
`--scoreG 60`

5. Increasing or decreasing minimum length of detected variant (default = 3 data points or bigger):
`--lengthG 1`

6. Increasing the number of CNVs you'd expect to get from that sample (default = 100, but for WGS samples better make it 1000 or even bigger, when such threshold is exceeded, quality threshold is increased and sample is re-analysed as many times as specified with `--maxNumIter` flag):
`--maxNumGermCNVs 2000 --maxNumIter 5`

7. Running the FDR correction (you specify number of iterations, higher the number is - more accurate FDR correction will be, but we do not recommend to make it bigger than 20 since it will slow down the calling):
`--fdrGermline 20`

8. Dividing your large cohort into several clusters of similar samples for more accurate parameter estimation and speeding the tool up. You specify the minimum size of the cluster (we recommend to keep it in between 20 and 100):
` --minimumNumOfElemsInCluster 50`

9. Speeding up the tool by using several cores (only some parts of the pipeline are parallelised):
`--numberOfThreads 4`

In the end you will have a command line like:

`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 1 --maxNumGermCNVs 2000 --maxNumIter 5 --fdrGermline 20 --minimumNumOfElemsInCluster 50 --numberOfThreads 4`

## WGS

For WGS you should skip step 3 and pay special attention to step 6 (may be put the maximum number of expected CNVs equal to 100.000 if you run your dataset for the first time and don't know what to expect). Probably it is also better to keep step 7 (number of iterations for FDR control at low level) since permutations of WGS data may take a while and even several permutations may be a good estimator of your data level of noise.


# How to interpret results

Your output files will be saved in `/your/folder/normal` directory - `ClinCNV` will make a separate folder per sample. 

There (if you did not specify flag `--visulizationIGV` as `F`) you can find 2 IGV tracks - one for real coverage values, second for CNVs track. Both tracks heights denote the underlying copy number (height = 0 denote homozygous deletion, height = 1 denote heterozygous deletion, height = 2 = diploid, height = 3 heterozyous duplication, etc.). Height is cut at level 6 since othervise differences would be much less illustrative (so, even if you have values with copy-number >6, they will be depicted just by points at height 6). You can just drag and drop these tracks into your IGV browser.

**Coverage values** are GC-normalised and median-normalised (that means that "average" coverage for cohort is around copy-number 2). But coverage is not square root normalised - so sometimes even if it seems that the dot is closer to alternative copy-number it may not be true since variances are different for different copy numbers at this plot.

**CNV-segments** here are mainly used for algorithmic diagnostic purposes - if the segments are too fragmented or if you clearly see some CNVs in non-polymorphic regions that were not detected by ClinCNV, then you've chosen quality threshold as too small/big, respectively, and you need to increase it and re-do the analysis. CNVs are described much more informatively in files "sampleName_cnvs.tsv"

![IGV tracks for one sample (exome seq)][IGV_track]



You can scroll into your data at chromosome level:


![IGV tracks for one sample (exome seq)][IGV_track_chr]

Or deeper at a single event level:

![IGV tracks for one sample (exome seq)][IGV_track_cnv]

## "sampleName_cnvs.tsv" files

Usually you will get tab delimited files with 2 header lines (starting with hashtag), column names description and the results itself.

First line - number of iterations. If the number of detected CNVs exceeded values specified with the flag `--maxNumGermCNVs`, your sample is re-analysed with stricter thresholds. The number here denotes number of re-runs.

Second line - number of outliers in autosomes. It is designed in a way that 5% of dots have to be "outliers" (their Z-score exceed the corresponding quantiles of standard normal) in a diploid sample. If you see a bigger number (such as 0.1 or bigger), that means that either your sample is largely affected by CNVs (or eg there is one aneuploidy) or the tool determined variances (or locations) wrongly. Unfortunately, you can do almost nothing in this case except dropping me an email. Value much below 5% may indicate same problems with parameters' estimation, but now variances were overestimated. It may lead to low sensitivity of the tool.

### Columns in the table




[IGV_track]: https://github.com/imgag/ClinCNV/raw/master/doc/images/germline_tracks.png "IGV tracks for germline sample"
[IGV_track_chr]: https://github.com/imgag/ClinCNV/raw/master/doc/images/germline_tracks_chrom_level.png "IGV tracks for germline sample (chromosome level)"
[IGV_track_cnv]: https://github.com/imgag/ClinCNV/raw/master/doc/images/germline_tracks_chrom_level.png "IGV tracks for germline sample (one CNV level)"