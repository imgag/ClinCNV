# How to run `ClinCNV` for trios analysis

The first question to answer is __why to run `ClinCNV` for trios__? Why not calling each sample separately? Calling variants separately, we use a fixed threshold usually. This threshold may allow a detection of variant in one of the three related samples analysed, but it may be too strict for detection of the same variant in another sample, even if it is presented there. Calling in trios does not improve the de novo variants calling, but it improves the concordance between the calls.

So here is the pipeline which is pretty similar to germline one.

## Pipeline for trios calling

1. Simple analysis, having only on-target coverage from the cohort of samples (generation of such files described in the documentation):
`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed --triosFile trios.txt`

2. Specifying output folder **(here and until the end of the numbered list instruction you should add the line to the line you got on the previous step)**:
`--out /your/folder`. Folder `normal` will be created in your output folders, results will be put into subfolders, one sample per subfolder.

3. Adding offtarget coverages (does not increase _Sensitivity_ a lot since germline CNVs are rarely as long as off-target window, but improves _Specificity_ efficiently removing CNVs formed by probes standing far away from each other, off-target coverage has to be extracted for *sufficiently large number of samples*, but not necessarily for all the samples from `normal.cov` file):
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

`Rscript clinCNV.R --normal normal.cov --bed bedFile.bed --triosFile trios.txt --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 1 --maxNumGermCNVs 2000 --maxNumIter 5 --fdrGermline 20 --minimumNumOfElemsInCluster 50 --numberOfThreads 4`



# How to interpret results

At the end you will get a table like this:

![Table with results][table_of_results]

The QC is still constant and equal to 1, in future it will be modified to something more meaningful.

The columns are:

`#chr	start	end` - are self explainable - they contain coordinates of detected variants.

`kid_CN_change` denotes a copy number of detected variant in the proband.

`mother_CN_change` and `father_CN_change` denote copy-numbers in parents.

`priority` is something we still work on. Priority 2 shows that the child has different copy-numbers from mother and father. Priority 1: at least one parent has the copy-number equal to the copy-number of the child. Priority 0: there is no copy-number change in the child.

`loglikelihood` denotes the score of the variant. Bigger is better.

`no_of_regions` shows how many data points were affected by this variant. Data point = coverage depth from 3 members of the trio for one particular targeted region.

`length_in_KB` shows the length of variant.

`loglik_kid`,  `loglik_mother`, `loglik_father` shows the individual contributions of each sample into the CNV log-likelihood score. If it is equal to 0, it means that the copy-number of the sample is unchanged. If it is smaller than 0 (it may happen, shown in the red cell in the table), it indicates that the data is contradictive to the inheritance model used. `ClinCNV` detects variants according to the inheritance patterns (e.g., father had copy-numbers 2 and 1 in his alleles, mother had 1 and 0, the kid may inherit one of the alleles from mother, one from father, which gives possible copy-numbers as 3, 2, 1 or 0, or have a de novo variant - but for this case we assume that copy-numbers of mother and father are equal to 2). Negative likelihood shows that the sample is more likely to be diploid, but the maximum likelihood was achieved with another model, taking into account the inheritance schemes. In the particular case, it looks like "homozygous deletion" in mother sample was detected due to lack of hybridization in this part in mother's sample (the kid's most probable state is diploid for this region).


[table_of_results]: ./images/trios_report.png "Table with results"