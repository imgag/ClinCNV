#!/usr/bin/env bash

set -e

folderWithNormals="/mnt/share/data/coverage/ssSC_v2"
folderWithTumors="/mnt/share/data/coverage/ssSC_v2-tumor"
folderWithEnrichmentKits="/mnt/share/data/enrichment/"
bedFile="ssSC_v2_2015_01_26.bed"
columnWhereCoverageStarts=4
outputFolder="/mnt/users/ahdemig1/clinCNV_development/ClinCNV_results/"
fileWithPairs="pairs.txt"
folderWithScript=$PWD
reanalyseCohort="F"
# New parameters (inroduced on 9th of August, 2018)
sampleNamesListGermline=""
sampleNamesListSomatic=""
nameOfTheAnalysis="exomes." # names of output files are changed accordingly
typeOfAnalysis="somatic" # two types possible - germline and somatic
scoreGermline=40 # threshold for calling Germline CNAs
lengthGermline=1 # minimum number of regions that forms a Germlinve CNV
scoreSomatic=60 # threshold for calling Somatic CNAs
lengthSomatic=5 # minimum number of regions that forms a Somatic CNA
maximumNumberOfGermlineCNVs=100 # this is a maximum amount of CNVs >=3KBps length expected in WGS 40x sample of European population
maximumNumberOfIterations=3 # tool increases thresholds if the number of CNVs exceeds max amount and re-analyse sample N times (specified)
maximumNumberOfSomaticCNAs=100 # has to be tuned


cd $folderWithScript

# prepare bed file
BedAnnotateGC -in $folderWithEnrichmentKits$bedFile -out $bedFile
BedAnnotateGenes -in $bedFile -out "annotated."$bedFile

# merge normal files 
if [[ ! -f $nameOfTheAnalysis"normal.txt" ]]; then
  if [[ $sampleNamesListGermline = "" ]]; then
    /mnt/share/opt/R-3.4.0/bin/Rscript --vanilla mergeFilesFromFolder.R -i $folderWithNormals -o $nameOfTheAnalysis"normal.txt" -n $columnWhereCoverageStarts 
  else
    /mnt/share/opt/R-3.4.0/bin/Rscript --vanilla mergeFilesFromFolder.R -i $sampleNamesListGermline -o $nameOfTheAnalysis"normal.txt" -n $columnWhereCoverageStarts 
  fi
else
  echo "File "$nameOfTheAnalysis"tumor.txt exists. We do not recalculate it. WARNING: if .bed file do not match coverage file, the pipeline may crash." 
fi 

if [[ $typeOfAnalysis = "somatic" ]]; then
  echo "Somatic framework is used. We create file with coverages for Tumors too."
  if [[ ! -f $nameOfTheAnalysis"tumor.txt" ]]; then
    if [[ $sampleNamesListSomatic = "" ]]; then
  	  /mnt/share/opt/R-3.4.0/bin/Rscript --vanilla mergeFilesFromFolder.R -i $folderWithTumors -o $nameOfTheAnalysis"tumor.txt" -n $columnWhereCoverageStarts 
    else
      /mnt/share/opt/R-3.4.0/bin/Rscript --vanilla mergeFilesFromFolder.R -i $sampleNamesListSomatic -o $nameOfTheAnalysis"tumor.txt" -n $columnWhereCoverageStarts \
      -p $fileWithPairs
    fi
  else
  	echo "File "$nameOfTheAnalysis"tumor.txt exists. We do not recalculate it. WARNING: if .bed file do not match coverage file, the pipeline may crash." 
  fi 
fi



# check if parameters are correctly specified
if [[ $typeOfAnalysis = "germline" ]]; then
	echo "Germline framework is used."
fi

if [[ $typeOfAnalysis != "somatic" && $typeOfAnalysis != "germline" ]]; then
	echo "Framework is not specified as somatic or germline. Quit."
fi


echo "Calling"
# run calling
if [[ $typeOfAnalysis = "somatic" ]]
then
	/mnt/share/opt/R-3.4.0/bin/Rscript --vanilla firstStep.R --normal $nameOfTheAnalysis"normal.txt" --tumor $nameOfTheAnalysis"tumor.txt" \
  --out $outputFolder --pair $fileWithPairs --bed "annotated."$bedFile --colNum $columnWhereCoverageStarts --folderWithScript $folderWithScript --reanalyseCohort TRUE \
  --scoreG $scoreGermline --lengthG $lengthGermline \
  --scoreS $scoreSomatic --lengthS $lengthSomatic \
  --maxNumGermCNVs $maximumNumberOfGermlineCNVs --maxNumIter $maximumNumberOfIterations --maxNumSomCNAs $maximumNumberOfSomaticCNAs
else
	/mnt/share/opt/R-3.4.0/bin/Rscript --vanilla firstStep.R --normal $nameOfTheAnalysis"normal.txt" --out $outputFolder --bed "annotated."$bedFile \
  --colNum $columnWhereCoverageStarts --folderWithScript $folderWithScript --reanalyseCohort TRUE \
  --scoreG $scoreGermline --lengthG $lengthGermline \
  --maxNumGermCNVs $maximumNumberOfGermlineCNVs --maxNumIter $maximumNumberOfIterations
fi
