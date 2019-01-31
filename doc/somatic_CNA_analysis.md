# How to run `ClinCNV` for somatic analysis

Saying _somatic_ we assume that for each tumor sample sequenced you also sequence a normal sample (it can be that several tumor samples were taken from one patient and only one normal was sequenced). We don't detect such events as _mosaicism_ even if it is theoretically possible - because we don't have a request for implementing such mode in `ClinCNV`, but if you need to detect mosaic events in "germline" samples - drop us an email.

Same as for _germline_ samples `ClinCNV` utilize on-target and (optionally) off-target coverage, but also `ClinCNV` may use _B-allele frequencies_ (BAFs) in _somatic_ mode. In brief, a CNV event may lead to changes in allele balance of point mutations (SNVs). The typical strategy on how to use signals from SNVs is 1) calculate read depths of SNVs for reference allele and alternative allele in both normal and tumor samples, 2) sub-select SNVs that are likely to be heterozgyous in normal tissue (which allele balance is close to 0.5), 3) measure the deviation of allele balance from normal to tumor. For example, if a sample has normally diploid region and it has a SNV with _C_ as reference and _T_ as alternative, we may expect 50% of reads, aligned to that position, to have _C_ in their sequence and 50% should have _T_ - thus we have something close to 0.5 as a B-allele frequency at this position (more accurately, allele balance _can be modelled_ with binomial distribution with probability of success close to 0.5). But if the region was duplicated (e.g. region containing _T_ was copied and we end up with 2 genomic regions containing _T_ and 1 genomic regions containing _C_ in each cancer cell), we expect ~33% of reads containing _C_ and ~66% containing _T_. Signal from B-allele frequencies helps a lot in somatic diagnositcs, it can distinguish between complex copy-number alterations and normal ones, only B-allele frequency change may indicate Loss of Heterozygosity events and it significantly improves the accuracy. So we highly recommend you to obtain BAF data from your samples, even if some additional efforts are required to prepare such data.

## Targeted sequencing

1. Simple analysis, having only on-target coverage from the cohort of tumor and paired normal samples with the information about pairing described in `pairs.txt` file (generation of such files described in the documentation):
`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed --tumor tumor.cov --pairs pairs.txt`

2. Specifying output folder **(here and until the end of the numbered list instruction you should add the line to the line you got on the previous step)**:
`--out /your/folder`

The folder `somatic` will be created and results for each samples pair will be put in the corresponding subfolder.

3. Adding off-target coverage (especially important for panel sequencing and overall playing bigger role in the analysis since sizes of CNAs are typically large)

`--normalOfftarget normalOff.cov --tumorOfftarget tumorOff.cov --bedOfftarget bedFileOff.bed`

4. Playing with _Sensitivity_ and _Specificity_ balance (increase of the threshold `--scoreS` leads to higher _Specificity_ and lower _Sensitivity_ and vice versa, default value is 100, however due to size of CNAs and high level of noise due to usage of FFPE samples it is recommended to keep this value high):
`--scoreS 150`

5. Increasing or decreasing minimum length of detected variant (measured in data points, default = 5 data points or bigger):
`--lengthS 10`

6. Including B-allele frequencies. Files with BAFs need to be located into same folder `/example/BAFs`:
`--bafFolder /example/BAFs`

Due to described above problems not all the samples may have BAF file, but you need to be sure that **at least one pair of tumor/normal samples has its BAF track** in the corresponding folder.

7. peeding up the tool by using several cores (only some parts of the pipeline are parallelised):
`--numberOfThreads 4`

## WGS

The same, but step 3 has to be skipped. 





# How to interpret results

You will have two output folders 1) IGV tracks and summary of BAF deviations on chromosome's arms level (located in your BAF folder in subfolder `/results`), 2) CNA calling results with IGV tracks and clonality plot (located in your output folder in subfolder `/somatic/`).

## BAF visualisation

### Barplots of deviated positions per chromosome arm
You can find there barplots showing number of chromosome arms with significant deviations in BAFs between tumor and normal samples. They show proportion of SNVs with significantly (p-value<0.05) different BAFs within the tumor-normal pair, thus, 1.0 indicates that 100% of SNVs have deviations in BAFs (thus can indicate presence of aneuploidy with high tumor content in the particular sample and chromosome arm), 0.0 indicated that 0% of SNVs have deviations. By random, we expect approximately 5% of SNVs having deviations, even if no CNAs happen in particular region.

![Barplot of deviated BAFs on chromosome arm level][BAF_barplot]

In this particular sample we can see some chromosome arms marked with red bars, which means they are not used for normalisation, but we can see that 1) those chromosome arms have small amount of SNVs (less than 10), 2) they are from chromosomes which short p arms are "empty" (below we will call "p" arm as left arm and "q" arm as right arm since they are oriented like this in human reference genome). We can conclude that it is highly likely that this sample does not have long CNAs that affect SNVs or such CNVs are presented only in clones with low clonality so the deviation is undistinguishable from normal. The advantage of such samples is that it is most likely efficiently normalised since all the chromosome arms are used for normalisation.

![Barplot of deviated BAFs on chromosome arm level][BAF_barplot_more]

In this particular example we can see that part of chr1 left arm is affected by CNVs (or there could be a aneuploidy with small tumor content) and even larger part of chr1 right arm affected by CNVs (or it could be aneuploidy with higher tumor content). Other chromosome arms may be analysed in a similar manner (chr13 right, chr14 right, chr16 right, chr17 left and right, chr18 right, chr19 left, etc.).

![Barplot of deviated BAFs on chromosome arm level][BAF_barplot_even_more]

In this sample almost all the chromosomes were likely to be damaged by CNAs. We take the least damaged chromosomes for the internal normalization (some of chromosome arms are indicated with orange color which means that they have >5% of BAFs significantly deviated, but we still had to take them into normalization procedure since the amount of "normal" material was too low), however you may expect higher level of noise in such samples due to small amount of "seemingly normal" material for normalization.

### IGV tracks

You can find 3 IGV BAF tracks per pair tumor/normal in the BAF folder.

First two tracks (named according to tumor and normal samples) contain just germline heterozygous positions and how did they change in tumor:
![IGV track of BAFs][BAF_track]

At the chromosome level it may look like depicted - you can see that dark red dots are signficantly deviated from germline ones in the right arm and somehow different in the left arm of the chromosome which may be caused by 1) clonality of the event happened with the left arm is lower than clonality of the event on the right, 2) event on the left is duplication with 3 alleles except 2 in normal while the event on the left is a heterozygous deletion or LOH event, 3) the event on the left may be a duplication of high copy number with the low clonality while the event on the right can be a duplication of extra high copy number with almost 100% purity of the tumor sample. Only looking at coverage plots (and running `ClinCNV` of course) may help you to understand what actually happened.

![IGV track of BAFs][BAF_track_chr]


The third track shows p-values less than 0.05. Red segment means "BAF for particular SNVs in this segment are significantly different between normal and tumor samples", blue segment shows "the difference is not significant at 0.05 level". As you can see, only bunch of SNVs may show us a CN change - there are a lot of CNAs not significant at 0.05 level while in general we can clearly see a shift. This plot is especially useful for sanity checks of `ClinCNV` results together with coverage plots. If red segments (or "saw" of small red pikes) match with the detected CNAs, then you can trust the results.

![IGV track of BAFs][BAF_track_plus_pvals]

For the "normal" genomic regions you may expect this plot to look like below. There is some evidence of a CNA on the right arm of the chromosome, however it is inconclusive without checking the coverage plot.

![IGV track of BAFs][BAF_track_plus_pvals_norm]

An investigator has to be especially careful with balanced copy number changes (thus, both alleles were duplicated or even triplicated) - they show no deviation in BAF patterns and only coverage may help you to solve this case. However **`ClinCNV` can not detect whole genome du- and triplications if no normal regions left in the genome**.



## IGV tracks and results folder

In folder `/your_output_folder/somatic/tumor_sample-paired_normal_sample` you can find IGV tracks, clonality plot and `.tsv` file with the results table.

### IGV tracks

There are 2 tracks - first is `.seg` file with CNAs detected and second is file with tumor copy number relative to normal copy number. In other words, both plots are "centered" around copy-number 2 since we assume that germline sample is diploid (except males and sex chromosomes where we assume copy number 1) and both plots depict tumor copy number. Again, internally we apply log-transformation, that's why vertical distances between segment may be misleading since depicted values are just ratios.

For tumor with high purity you could see something like this:

![IGV track of copy number][IGV_tracks_high_clonality]

For low purity tumor segments are much less distant from the horizontal line `y=2`.

![IGV track of copy number][IGV_tracks_low_clonality]

Sometimes you can see copy number segments which are located exactly at `y=2` line. These are lines indicating copy-neutral loss-of-heterozygosity. Usage of BAF tracks and coverage tracks together is necessary for visual inspection of LOH events.



### Clonality plots

Each tumor sample may contain sub-clones - parts of tumor that share some unique variants between them which are not presented in samples that do not belong to this particular sub-clone. Further sub-clones may be divided into more sub-clones, etc. `ClinCNV` does not try to build a tree of somatic evolution for sub-clones, however it is crucially important to identify all these sub-clones and classify detected CNAs according to the abundance of such sub-clones.

`ClinCNV` divides all the possible purities into intervals (by default from 0.05 to 1.0 with the step 0.025) and identifies from 1 up to 6 distinct sub-clone purities that help to explain big parts of variation in detected CNAs. `ClinCNV` does not try to estimate clonal structure if no CNAs were detected.

The first plot - heatmap of the likelihood landscape of variants detected in the sample, red color denotes "pretty unrealistic", blue - the opposite, darker is more likely. On y-axis you can see possible purities of the biggest clone, on x-axis - possible purity of the second biggest clone. You may see 2 main types of heatmaps: "stripe" and "ellipsoid":

![Heatmap of likelihood landscape][stripe_low_clone]  ![Heatmap of likelihood landscape][stripe_high_clone]


When more than 2 clones are presented 2D heatmap representation




[BAF_barplot]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_barplot1.png "Barplot of deviated BAFs on chromosome arm level"
[BAF_barplot_more]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_barplot2.png "Barplot of deviated BAFs on chromosome arm level"
[BAF_barplot_even_more]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_barplot3.png "Barplot of deviated BAFs on chromosome arm level"
[BAF_barplot_even_more]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_barplot3.png "Barplot of deviated BAFs on chromosome arm level"
[BAF_track]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_BAFs_tumor_normal.png "IGV track of BAFs"
[BAF_track_chr]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_BAF_chr_level.png "IGV track of BAFs"
[BAF_track_plus_pvals]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_BAF_pvalue.png "IGV track of BAFs"
[BAF_track_plus_pvals_norm]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_BAF_pvalue_normal.png "IGV track of BAFs"
[IGV_tracks_high_clonality]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_IGV_tracks_high_clonality.png "IGV track of copy number"
[IGV_tracks_low_clonality]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_IGV_tracks_low_clonality.png "IGV track of copy number"
[stripe_low_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/stripe_low_clone.png "Heatmap of likelihood landscape"
[stripe_high_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/stripe_high_clone.png "Heatmap of likelihood landscape"

