#!/usr/bin/env bash



folderWithNormals="/mnt/share/data/coverage/ssSC_v2"
folderWithTumors="/mnt/share/data/coverage/ssSC_v2-tumor"
folderWithEnrichmentKits="/mnt/share/data/enrichment/"
bedFile="ssSC_v2_2015_01_26.bed"
columnWhereCoverageStarts=4
outputFolder="/mnt/share/opt/clincnv-somatic-1.0/CNV_results/"
fileWithPairs="pairs.txt"
folderWithScript=$PWD
reanalyseCohort="F"

cd $folderWithScript

# prepare bed file
BedAnnotateGC -in $folderWithEnrichmentKits$bedFile -out $bedFile
BedAnnotateGenes -in $bedFile -out "annotated."$bedFile

# merge normal files 
/mnt/share/opt/R-3.4.0/bin/Rscript --vanilla mergeFilesFromFolder.R -i $folderWithNormals -o normal.txt -n $columnWhereCoverageStarts 
/mnt/share/opt/R-3.4.0/bin/Rscript --vanilla mergeFilesFromFolder.R -i $folderWithTumors -o tumor.txt -n $columnWhereCoverageStarts 

# run calling
/mnt/share/opt/R-3.4.0/bin/Rscript --vanilla firstStep.R --normal normal.txt --tumor tumor.txt --out $outputFolder --pair $fileWithPairs --bed "annotated."$bedFile --colNum $columnWhereCoverageStarts --folderWithScript $folderWithScript --reanalyseCohort TRUE

