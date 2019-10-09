# How to run `ClinCNV` for somatic analysis

Saying _somatic_ we assume that for each tumor sample sequenced you also sequence a normal sample (it can be that several tumor samples were taken from one patient and only one normal was sequenced). We don't detect such events as _mosaicism_ even if it is theoretically possible - because we don't have a request for implementing such mode in `ClinCNV`, but if you need to detect mosaic somatic events - drop us an email.

Same as for _germline_ samples `ClinCNV` utilize on-target and (optionally) off-target coverage (explanation can be found [there](https://www.cell.com/trends/genetics/fulltext/S0168-9525(13)00127-3 ), but also `ClinCNV` may use _B-allele frequencies_ (BAFs) in _somatic_ mode. In brief, a CNV event may lead to changes in allele balance of point mutations (SNVs). The typical strategy on how to use signals from SNVs is 1) calculate read depths of SNVs for reference allele and alternative allele in both normal and tumor samples, 2) sub-select SNVs that are likely to be heterozgyous in normal tissue (which allele balance is close to 0.5 - half of reads support one allele and half support the alternative), 3) measure the deviation of allele balance from normal to tumor. For example, if a sample has normally diploid region and it has a SNV with _C_ as reference and _T_ as alternative, we may expect 50% of reads, aligned to that position, to have _C_ in their sequence and 50% should have _T_ - thus we have something close to 0.5 as a B-allele frequency at this position (more accurately, allele balance _can be modelled_ with binomial distribution with probability of success close to 0.5). But if the region was duplicated (e.g. region containing _T_ was copied and we end up with 2 genomic regions containing _T_ and 1 genomic regions containing _C_ in each cell), we expect ~33% of reads containing _C_ and ~66% containing _T_. Signal from B-allele frequencies helps a lot in somatic diagnositcs, it can help to detect allele-specific copy number changes, only B-allele frequency change may indicate Loss of Heterozygosity events and it significantly improves the accuracy of somatic calling. So we highly recommend you to obtain BAF data from your samples, even if some additional efforts are required to prepare this data.

## Targeted sequencing

1. Simple analysis, having only on-target coverage from the cohort of tumor and paired normal samples with the information about pairing described in `pairs.txt` file (generation of such files described in the documentation):
`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed --tumor tumor.cov --pairs pairs.txt`

2. Specifying output folder **(here and until the end of the numbered list instruction you should add the line to the line you got on the previous step)**:
`--out /your_output_folder`

The folder `somatic` will be created and results for each samples pair will be put in the corresponding subfolder.

3. Adding off-target coverage (off-target reads are beneficial in somatic context comparing to germline)

`--normalOfftarget normalOff.cov --tumorOfftarget tumorOff.cov --bedOfftarget bedFileOff.bed`

4. Playing with _Sensitivity_ and _Specificity_ balance (increase of the threshold `--scoreS` leads to higher _Specificity_ and lower _Sensitivity_ and vice versa, default value is 100, however due to size of CNAs and high level of noise due to usage of FFPE samples it is recommended to keep this value high):
`--scoreS 150`

5. Increasing or decreasing minimum length of detected variant (measured in data points, default = 10 data points or bigger):
`--lengthS 4`

6. Including B-allele frequencies. Files with BAFs need to be located into same folder `/example/BAFs` (tumor and normal mixed, no need to keep a special directory structure):
`--bafFolder /example/BAFs`. Due to different issues not all the samples may have BAF file, but you need to be sure that **at least one pair of tumor/normal samples has its BAF track** in the corresponding folder.

7. Speeding up the tool by using several cores (only some parts of the pipeline are parallelised):
`--numberOfThreads 4`

8. In case if you have prior knowledge on number of clones in your tumor and the clonal inference of `ClinCNV` was not accurate you can specify the penalty for additional clones (less penalty = more clones but up to 5, default penalty is 300):
`--clonePenalty 500`

9. Variants' QC filtering. This parameter may have 2 values: 0, 1 or 2. 0 means that no QC filtering of variants is performed. 1 means that QC filtering is performed only at the first step of the algorithm (clonal structure inference). 2 means that even final calls will be QC filtered. It is recommended to set it equal to 2 for WES and to 1 for gene panel sequencing.
`--filterStep 0`

10. Usually, starting from some cancer cell fraction and below this value the results become non trustable (most of them are false positives). This value depends on the quality of your data: for high coverage WES variants of 5% CCF can be called while for low-coverage targeted panel sequencing better keep it above 20-30%. For example, if you do not trust variants with cancer cell fraction below 20%, you need to specify it as:
`--clonalityForChecking 0.2`
`ClinCNV` will apply additional filterings for such variants. No allele balanced variants will be called below 20% since we need BAF signal to perform QC control of variants. If you absolutely do not trust variants with cancer cell fraction below 20%, you may specify:
`--minimumPurity 20`

11. If you want to call only one particular pair of samples, you need to specify their IDs:
`--normalSample IDofNormal --tumorSample IDofTumor`

12. In case some particular sample was called incorrectly and you see a large stretch of homozygous deletion (>10MBs), you may point `ClinCNV` to this mistake by providing a parameter guiding a diploid baseline:
`--guideBaseline chrN:X-Y`
where chrN:X-Y is the coordinate of this long stretch of homozygous deletion. `ClinCNV` will consider this part as diploid and re-call the sample. (Don't forget to specify samples' IDs from step 11).

13. If you do not want to bother with the interpretation of 2 rounds of possible CNVs (or your tumor samples are untreated and early stage), you may turn the support of the 2nd round of CNVs off:
`--notComplexTumor`

## WGS

The same, but step 3 has to be skipped. **WE HAVE NOT TESTED CLINCNV FOR WGS YET**. It may be time consuming. Try to increase the clonality step:

`--purityStep 5` instead of default `--purityStep 2.5`.





# How to interpret results

If you **have time** (e.g., you are analysing a single sample in clinical context), you go through the pictures and table with variants. You may recall the sample if you do not like the results.

If you **do not have time** but need to analyse thousands of samples, you identify suspicious samples and then re-analyse them or remove from the analysis. Suspicious samples are the ones with large FDR (QC value provided in the header of `ClinCNV` output files, more than 10% is a lot), the ones that have a lot of CNAs large part of them is LOH variants (it may happen when samples' labels are mixed up), samples with many short homozygous deletions. If a sample has large FDR, but less than 20%, it is usually possible just to filter out false positive variants (the ones with no deviation in B-allele frequencies). If it is bigger, may be it is easier to remove the sample. If there is a long (>10MB) homozygous deletion in the results, it is a sign that this sample is actually tetraploid (or more) and you need to apply step 12 from the list above to re-call this sample. And if you need to analyse hundreds of samples - scroll this doc until the header "I do not have time".

## I have time!

You will have two output folders 1) IGV tracks and summary of BAF deviations on chromosome's arms level (located in your BAF folder in subfolder `/results`), 2) CNA calling results with IGV tracks and clonality plot (located in your output folder in subfolder `/somatic/`).


### BAF visualisation

_This part is explanation of the weird plots you will get in the folder with `baf` files you've specified. If you run `ClinCNV` for the first time, you may skip this part._

#### Barplots of deviated positions per chromosome arm
You can find there barplots showing number of chromosome arms with significant deviations in BAFs between tumor and normal samples. They show proportion of SNVs with significantly (p-value<0.01) different BAFs within the tumor-normal pair, thus, 1.0 indicates that 100% of SNVs have deviations in BAFs (thus can indicate presence of aneuploidy with high tumor content in the particular sample and chromosome arm), 0.0 indicated that 0% of SNVs have deviations. By random, we expect less than 5% of SNVs having deviations, even if no CNAs happen in particular region - since quite a lot of SNVs are false positives due to alignment problems, this number is not equal to 5%.

![Barplot of deviated BAFs on chromosome arm level][BAF_barplot]

In this particular sample we can see some chromosome arms marked with red bars, which means they are not used for normalisation, but we can see that 1) those chromosome arms have small amount of SNVs (less than 10 and we mark such chromosome arms as "not suitable for normalization" since we know nothing about their CNA status - as well as sex chromosomes but for simplicity reasons), 2) they are from chromosomes which short p arms are "empty" (below we will call "p" arm as left arm and "q" arm as right arm since they are oriented like this in human reference genome). We can conclude that it is highly likely that this sample does not have long CNAs that affect SNVs or such CNVs are presented only in clones with low clonality so the deviation is undistinguishable from normal. The advantage of such samples is that they can be (most likely) efficiently normalised since all the chromosome arms (except "empty" ones and sex chromosomes) are used for normalisation, and more data leads to more accurate estimation of statistical parameters.

![Barplot of deviated BAFs on chromosome arm level][BAF_barplot_more]

In this particular example we can see that part of chr1 left arm is affected by CNVs (or there could be a aneuploidy with small tumor content) and even larger part of chr1 right arm affected by CNAs (or it could be aneuploidy with higher tumor content). Other chromosome arms may be analysed in a similar manner (chr13 right, chr14 right, chr16 right, chr17 left and right, chr18 right, chr19 left, etc.).

![Barplot of deviated BAFs on chromosome arm level][BAF_barplot_even_more]

In this sample almost all the chromosomes were likely to be damaged by CNAs. We take the least damaged chromosomes for the internal normalization (some of chromosome arms are indicated with orange color which means that they have >5% of BAFs significantly deviated, but we still had to take them into normalization procedure since the amount of "normal" material was too low), however you may expect higher level of noise in such samples due to small amount of "seemingly normal" material for normalization.

#### IGV tracks of BAFs

_This section describes the BAF tracks generated for IGV visualization. May be skipped in the beginning._

You can find 3 IGV BAF tracks per pair tumor/normal in the BAF folder.

First two tracks (named according to tumor and normal samples) contain just germline heterozygous positions (green) and how did they change in tumor (dark red):
![IGV track of BAFs][BAF_track]

At the chromosome level it may look like depicted below - you can see that dark red dots are signficantly deviated from germline ones in the right arm and somehow different in the left arm of the chromosome which may be caused by 1) clonality of the event happened with the left arm is lower than clonality of the event on the right, 2) event on the left is duplication with 3 alleles except 2 in normal while the event on the right is a heterozygous deletion or LOH event, 3) the event on the left may be a duplication of high copy number with the low clonality while the event on the right can be a duplication of extra high copy number with almost 100% purity of the tumor sample. Only looking at coverage plots (and running `ClinCNV` of course) may help you to understand what actually happened.

![IGV track of BAFs][BAF_track_chr]


The third track shows p-values less than 0.05. Red segment means "BAFs for particular SNVs in this segment are significantly different between normal and tumor samples", blue segment shows "the difference is not significant at 0.05 level". As you can see on the left part of the track, only bunch of SNVs may show us a CN change - there are a lot of CNAs not significant at 0.05 level while in general we can clearly see a shift. This plot is especially useful for sanity checks of `ClinCNV` results together with coverage plots. If red segments (or "saw" of small red pikes) match with the detected CNAs, then you can trust the results.

![IGV track of BAFs][BAF_track_plus_pvals]

For the "normal" genomic regions you may expect this plot to look like below - almost straight line (5% of SNVs will show p-value < 0.05 at random). There is some evidence of a focal CNA on the right arm of the chromosome (red pike), however it is inconclusive without checking the coverage plot.

![IGV track of BAFs][BAF_track_plus_pvals_norm]

An investigator has to be especially careful with balanced copy number changes (thus, both alleles were duplicated or even triplicated) - they show no deviation in BAF patterns and only coverage may help you to solve this case. However **`ClinCNV` can not detect whole genome du- and triplications if no normal regions left in the genome**.



### IGV tracks of coverage and BAFs and results folder

_This step is useful if you go into in-depth analysis of each CNA in IGV. Useless otherwise (if you want to assess the general quality of the calls)._

In folder `/your_output_folder/somatic/tumor_sample-paired_normal_sample` you can find IGV tracks, clonality plot and `.tsv` file with the results table.

#### IGV tracks

There are 2 tracks - first is `.seg` file with CNAs detected and second is file with tumor copy number relative to normal copy number. In other words, both plots are "centered" around copy-number 2 since we assume that germline sample is diploid (except males and sex chromosomes where we assume copy number 1) and both plots depict tumor copy number. Again, internally we apply log-transformation, that's why vertical distances between segment may be misleading since depicted values are just ratios.

For tumor with high purity you could see something like this:

![IGV track of copy number][IGV_tracks_high_clonality]

For low purity tumor segments are much less distant from the horizontal line `y=2`.

![IGV track of copy number][IGV_tracks_low_clonality]

Sometimes you can see copy number segments which are located exactly at `y=2` line. These are lines indicating copy-neutral loss-of-heterozygosity. Usage of BAF tracks and coverage tracks together is necessary for visual inspection of LOH events.

If you visualise both BAF and coverage tracks you should see something like:

![IGV track of copy number][IGV_tracks_all]

It is important to check visually if the predictions match both coverage and BAF patterns. At this plot not the most complicated tumor is shown - sometimes literally whole sample is detected as copy number alterated (with different events there), and, unfortunately, `ClinCNV` does not contain all the answers to the tumor structure. You still have to investigate complex cases yourself.

At the chromosome level a complex case may look like:

![IGV track of copy number][IGV_tracks_chrom]

Note - at the p-arm (left) B-allele frequency changed from 0.5 to almost 0/1, but the coverage increased. If a "classic" duplication happens, this region has to be duplicated around 10 times (BAF balance will move from 0.5 to ~0.1-0.9), and we don't see such a huge increase in coverage. Thus, this coverage/BAF pattern may indicate loss-of-heterozygosity event followed by a low-copy duplication. `ClinCNV` determines this event as 

#chr	|start|	end|	tumor_CN_change|	tumor_clonality|	CN_change|	loglikelihood|	number_of_regions|	state
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
chr7	|688268|	55208953|	3|	0.8|	2.8|	2646|	642|	LOHDup

More on this table and potentially detectable events below.


### Sub-clonal structure of samples

#### Clonality plots

_This part describes clonal structure of sample and the plots describing the sub-clonal structure generated by ClinCNV_

Each tumor sample may contain sub-clones - parts of tumor that share some unique variants between them which are not presented in samples that do not belong to this particular sub-clone. Further sub-clones may be divided into more sub-clones, etc. `ClinCNV` does not try to build a tree of somatic evolution for sub-clones, however it is crucially important to identify all these sub-clones and classify detected CNAs according to the abundance of such sub-clones.

`ClinCNV` divides all the possible purities into intervals (by default from 0.05 to 1.0 with the step 0.025) and identifies from 1 up to 6 distinct sub-clone purities that help to explain big parts of variation in detected CNAs. `ClinCNV` does not try to estimate clonal structure if no CNAs were detected.

The first plot - heatmap of the likelihood landscape of variants detected in the sample, red color denotes "pretty unrealistic", blue - the opposite, darker is more likely. On y-axis you can see possible purities of the biggest clone, on x-axis - possible purity of the second biggest clone. You may see 2 main types of heatmaps: "stripe" and "ellipsoid":

Low purity "stripe" pattern            |  High purity "stripe" pattern
:-------------------------:|:-------------------------:
![Heatmap of likelihood landscape][stripe_low_clone]  |  ![Heatmap of likelihood landscape][stripe_high_clone]
 
Stripe pattern basically means that, according to `ClinCNV`, this tumor contains only 1 clone or 1 major clone + 1 subclone of that major clone, and subclone inside the major clone is large (at least 80% of the major clone) - if the sublone is smaller, we will see "ellipsoid" pattern:

Low purity "ellipsoid" pattern            |  High purity "ellipsoid" pattern
:-------------------------:|:-------------------------:
![Heatmap of likelihood landscape][ellipsoid_low_clone]  |  ![Heatmap of likelihood landscape][ellipsoid_high_clone]

Ellipsoid pattern means that there are at least 2 clones and 1 clone is not enough to explain the tumor's variants. All these rules described are not always true. You may see multiple clones and elliplsoid pattern, however it will most probably mean that the minor clones have small length and low likelihood.

Rarely you may see other complex patterns, most often they happen when 2 clones is not enough. In this case you should look at the barplot in your results directory.

Example of complex pattern         |  Example of complex pattern
:-------------------------:|:-------------------------:
![Heatmap of likelihood landscape][strange_pattern]  |  ![Heatmap of likelihood landscape][strange_pattern2]

When more than 2 clones are presented 2D heatmap representation becomes non relevant. You can check the full clonal structure at barplots.

At X axis of this barplot you can see all the potential purities we've investigated. At Y axis - length in MB. Each purity has 4 bars: brown for copy-number neutral changes (loss of heterozygosity), blue - for duplication of copy number 3 and 4, dark blue - duplications of high copy number, red - for deletions. White marks at barplot shows lengths of variants. Number on top of the bar shows the number of variants.

In this particular tumor we can see one major clone with purity of 0.325. It contains 3 duplication of copy number 3 or 4 (blue bar) and 3 deletions (red bar).
![Barplot of clonal landscape][barplot_one_clone] 

In this tumor we can see 2 clones: major with purity 0.775 and minor with purity 0.175. Small clone contains only duplications of small copy number (7 in total) and big clone - only deletions (four cases).
![Barplot of clonal landscape][barplot_two_clone]

And in this example we see 3 clones at 0.3, 0.4 and 0.725, all of them have different variants, only the smallest clone has duplications of high copy numbers (>4) and all clones have copy-number neutral LOH events.
![Barplot of clonal landscape][barplot_multiple_clone]

`ClinCNV` was never intended to be a tool for clonal evolution reconstruction, however, it makes corrections according to number of clones. But `ClinCNV` can still make mistakes. If `ClinCNV` tells you that there are 2 clones and you are absolutely sure that there should be only 1 clone, you can increase the penalty for a new clone with the flag, e.g.: `--clonePenalty 500` instead of default 200. Decrease this parameter if you want the tool to be more sensitive to different clones. This parameter does not affect the heatmap explained above - use it as a guide for choosing the penalty for additional clones in complex cases.


### CNAs visualization plots

_This part is fun to read (and important)_

As the result, `ClinCNV` produces several types of plots. At first, this plot shows the CNAs at different chromosomes. Some researchers prefer this visualization.

![Chromosomal plot][chrom_plot]

Red color denotes deletion, blue color denotes duplication, the color intensity shows the sub-clonal fraction of tumor cells harboring this variant. We can see that at chr11 the deletion was fragmented (same cancer cell fraction but multiple deletions occur) - this is not usual for `ClinCNV` and usually it shows that something real happened there (e.g., microduplications, so the non-deleted regions are, technically, LOH), however, it may still be a technical artifact. Increase the minimum length of the detected variant to get rid of this if you don't like it (step 5 from the list of commands above).

We also provide a plot of allele-specific copy-number changes. The vertical lines show the borders of chromosomes (and chromosome arms), horizontal purple lines on two bottom panels show integer copy-numbers.

![Allele specific plot][allele_cna_plot]

This is the same sample as before, but in another representation. The top panel shows B-allele frequency track (more accurately, the maximum likelihood of the closest possible BAF), the middle panel shows coverage, the bottom panel shows allele decomposition. Let's discuss this plot a bit more.

**The colors.** __Black points__ denotes data points (it does not distinguish between off-target and on-target data points). It may be BAF (top panel) or log-ratio of coverages (middle panel). Different colors at 2 tops panels denote different type of CNAs: brown colors indicates the presence of LOH, blue color shows presence of duplication, red color represents deletion. Lines are drawn "where it is expected". If you see a CNA when the lines do not intersect with the black dots - don't trust these calls, something different happened there. What could happen is the 2nd round of CNAs - when one CNA happens and then some other happen within the borders of the previously found CNA - but we will discuss it later.

**The bottom panel.** For duplications when one of the alleles remained normal (one copy), we chose to show it with green color. Everything that was duplicated is depicted with blue color. Red color shows the deletion. May be it would be more natural to show it with green color (since one allele remained "normal"), but we have tried that color scheme and decided to show it with red - in order to get the attention. **Horizontal black bars** show allelic copy-numbers. Let's look at chr5. There is a horizontal bar at copy-number 1 and the same bar at copy-number 3 - which means that total copy-number is equal to 3 and the minor allele has copy-number 1. Chr9? One allele was deleted so the black bar dropped to 0 copies, but another bar is a copy-number 2 - it indicates LOH event. Intensity of the color shows the cancer cell fraction of the variant.

_Easy, yes?_

Let's move to more complicated examples.

![Allele specific plot][somatic_2nd_round]

I recommend to open this file in a separate window and inspect the variants manually. It is a panel sequenced sample, almost all these black dots of coverage are off-target data. Let's check the chr2. It looks inconsistent! Coverage profile at chr2q shows that there are several deletions, but the allelic plot shows there is only one! And BAF lines also did not change!

The answer is: **several rounds of CNAs happened in chr2q**. So there was a duplication __after__ a deletion -- that's why the coverage profile changed, but allelic decomposition did not. But why it is not a separate CNA you'll ask - and you'll be right! I suggest you to look at chr5q. Coverage drop here indicates a small deletion. But BAF strongly suggests: it is A BIG DELETION in a large sub-clone. **Whom will we trust?** They both are right (at least be believe so). Probably LOH happened first and then - deletion of this part. **You can turn a second round of CNAs calling off**. Step 13. Doctors don't like this 2nd round of CNAs -- but otherwise for some samples calls will be segmented: it will be LOH for SNV-rich regions and deletion for SNV-free regions. So may be better turn it off until you will see something really weird. `ClinCNV` does not reconstruct evolutionary history of changes and it may miss the nested structure of CNAs - the tool is created for accurate segmentation and calling, other tasks should be performed during post-processing. Moreover, 2nd round of CNAs is downweighted comparing to more simplistic explanation - but sometimes data says this.

Tetraploidy and something could also be an explanation - however, we can not be sure. We prefer to assume that the tumor is mostly diploid and recover the clonal structure with this complex modelling rather than increasing ploidy. `ClinCNV` works with tetraploid tumors under condition that at least part of the sample remained diploid.

Even though, `ClinCNV` may be terribly wrong. E.g. in this situation (pictured below) due to not so big number of SNVs `ClinCNV` allowed extremely long homozygous deletions. This sample has to be recalled with step 12: `--guideBaseline chrN:X-Y`, where baseline should be at the same level as homozygous deletions. We can see that chr4 is "homozygously deleted" in this call - so put `--guideBaseline chr4:54506451-190696238` -- approximate coordinate of the homozygous deletion. Of course you may use homozygous deletions from chroms 10 or 18 as the baseline - but better choose it as long as possible.

![Allele specific plot - wrong calls!][clincnv_wrong]

For correction of mistakes not only the baseline can be changes, but also filtering step (step 9). In general, QC filtering of variants improves calling of variants, but if the tumor has at least one small cancer cell fraction sub-clone - this filtering may distort the calls. Try to turn it off (set to 0).


### Results table

Results of somatic CNA calling are summarised in a table:

![Table of results][results_table]

Top three lines just shows type of the analysis, the version used and when the analysis was finished and may be ignored.

The fourth line is informative. At first, _QC_ shows predicted False Discovery Rate. It is measured as the percentage of CNVs that were detected as allele imbalanced, but show no deviation in BAFs. _gender of sample_ is, actually, gender of sample. Then there is estimated ploidy (rarely equal exactly to to 2 even if no CNAs were detected since it is estimated from coverage and it is, after all the transformations, continuous-like). Then clonality -- it shows different sub-clones' cancer cell fraction. The biggest number is considered as _purity_ of the tumor. Other numbers are given relatively to all the DNA sequenced, thus, not normalized by purity.

Next line has column names. What they denote:


First three columns of the table show coordinates of the variant.

4th column - _"tumor_CN_change"_ - shows a copy number after correction on percentage of this clone in sequenced tumor mass. Usually this number is the most interesting for the researchers.

5th and 6th - _"major_CN"_ and _"minor_CN"_ shows allelic decomposition of copy-number alteration.

_"tumor_clonality"_ shows to which subclone this particular variant belongs.

_"CN_change"_ is an absolute copy-number change. At first `ClinCNV` finds an absolute copy-number change and then it infers "tumor_CN_change" and "tumor_clonality". This value may also be reported if you don't want to interpret the clonal structure.

_"loglikelihood"_ shows the quality score of particular variant, bigger the score is - more credible is the variant.

_"median_loglikelihood"_ shows the median of log-likelihood. If it is negative - then this variant was caused by several outliers and is not real. It is 0 for LOH since the BAF signal is sparse.

_"number_of_regions"_ is the number of datapoints within the variant (on- and off-target regions). More datapoints usually means more credible variant.

_"major_CN_allele2"_ and _"minor_CN_allele2"_ denotes "the second round of CNAs that probably happened". Always empty when you turn this mode off (step 13). 

_"tumor_clonality_2"_ denotes the percentage of cells that experienced the second round of CNAs. It is **always** smaller than the first Cancer Cell Fraction (well, it is expected).

_"genes"_ field is empty if you have not annotated your `bed` file with genes. 

Next 4 fields (_"ontarget_RD_CI_lower"_ etc) shows the bounds of 95% confidence interval for read depth.

_"Lowmed_tumor_BAF"_ and _"Highmed"_ shows medians of tumor BAFs which are below the median of normal BAFs and above. It is useful for diagnostics of small variants (p-value may be non-significant, but the shift could be seen).

_"BAF_qval_fdr"_ is a q-value for BAFs within the variant. Is expected to be low for allele-imbalanced variants and high for short variants with only couple of SNVs inside or allele-balanced variants (e.g., CN2 of one allele and CN2 of another allele).

_"Overall_qvalue"_ is q-value obtained from merged p-values from coverage and BAFs.


## I do not have time! I have hundreds of samples!



# Still have questions?

Ask them at `german dot demidov at medizin dot uni-tuebingen dot de`. We are glad to help you with everything connected to `ClinCNV`!



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
[IGV_tracks_all]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_IGV_tracks_all.png "IGV track of copy number"
[IGV_tracks_chrom]: https://github.com/imgag/ClinCNV/raw/master/doc/images/somatic_IGV_tracks_chrom.png "IGV track of copy number"

[stripe_low_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/stripe_low_clone.png "Heatmap of likelihood landscape"
[stripe_high_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/stripe_high_clone.png "Heatmap of likelihood landscape"
[ellipsoid_low_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/ellipsoid_low_clone.png "Heatmap of likelihood landscape"
[ellipsoid_high_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/ellipsoid_high_clone.png "Heatmap of likelihood landscape"
[strange_pattern]: https://github.com/imgag/ClinCNV/raw/master/doc/images/strange_pattern.png "Heatmap of likelihood landscape"
[strange_pattern2]: https://github.com/imgag/ClinCNV/raw/master/doc/images/strange_pattern2.png "Heatmap of likelihood landscape"
[barplot_one_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/barplot_one_clone.png "Barplot of clonal landscape"
[barplot_two_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/barplot_two_clone.png "Barplot of clonal landscape"
[barplot_multiple_clone]: https://github.com/imgag/ClinCNV/raw/master/doc/images/barplot_multiple_clone.png "Barplot of clonal landscape"

[chrom_plot]: ./images/CNAs_plot_other_form.png "Chromosomal plot of CNAs"
[allele_cna_plot]: ./images/CNAs_plot.png "Allele-specific plot of CNAs"
[somatic_2nd_round]: ./images/somatic_2nd_round.png "2nd round of somatic CNAs"
[clincnv_wrong]: ./images/clincnv_wrong.png "ClinCNV wrong"

[results_table]: https://github.com/imgag/ClinCNV/raw/master/doc/images/results_table.png "Table of results"