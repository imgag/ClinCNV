# Installation of `ClinCNV`

Open `R` and install the following packages:

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
Rscript $fold"/clinCNV.R" --bed $fold"/samples/bed_file.bed" --normal $fold"/samples/coverages_normal.cov"
```

`ClinCNV` is a set of `R` scripts. To install all the necessary packages execute [the following script](https://github.com/imgag/megSAP/blob/master/data/install_deps_clincnv.R).

For files' preparation we recommend to use [`ngs-bits`](https://github.com/imgag/ngs-bits). We recommend to install `ngs-bits` using [`bioconda`](https://github.com/imgag/ngs-bits/blob/master/doc/install_bioconda.md). 

