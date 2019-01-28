# How to run `ClinCNV` for germline analysis

## Tagret sequencing

* Simple analysis, having only on-target coverage:

`R clinCNV.R --normal normal.cov --bed bedFile.bed`

* Specifying output folder:

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder`

* Addinog offtarget coverages (does not increase _Sensitivity_ a lot since germline CNVs are rarely as long as off-target window, but improves _Specificity_ efficiently removing CNVs formed by probes standing far away from each other):

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed`

* Playing with _Sensitivity_ and _Specificity_ balance (increase of the threshold `--scoreG` leads to higher _Specificity_ and lower _Sensitivity_ and vice versa, default value is 40):

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60`

* Increasing or decreasing minimum length of detected variant (default = 3 data points or bigger):

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 1`

* Increasing the number of CNVs you'd expect to get from that sample (default = 100, but for WGS samples better make it 1000 or even bigger, when such threshold is exceeded, quality threshold is increased and sample is re-analysed as many times as specified with `--maxNumIter` flag):

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 1 --maxNumGermCNVs 2000 --maxNumIter 5`

* Running the FDR correction (you specify number of iterations, higher the number is - more accurate FDR correction will be, but we do not recommend to make it bigger than 20 since it will slow down the calling):

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 1 --maxNumGermCNVs 2000 --maxNumIter 5 --fdrGermline 20`

* Dividing your large cohort into several clusters of similar samples for more accurate parameter estimation and speeding the tool up. You specify the minimum size of the cluster (we recommend to keep it in between 20 and 100):

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 1 --maxNumGermCNVs 2000 --maxNumIter 5 --fdrGermline 20 --minimumNumOfElemsInCluster 50`

* Speeding up the tool by using several cores (only some parts of the pipeline are parallelised):

`R clinCNV.R --normal normal.cov --bed bedFile.bed --out /your/folder --normalOfftarget normalOff.cov --bedOfftarget bedFileOff.bed --scoreG 60 --lengthG 1 --maxNumGermCNVs 2000 --maxNumIter 5 --fdrGermline 20 --minimumNumOfElemsInCluster 50 --numberOfThreads 4`



## WGS

Same simple scenario:

`R clinCNV.R --normal normal.cov --bed bedFile.bed`

# How to interpret results