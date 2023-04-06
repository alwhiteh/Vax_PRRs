#!/bin/bash

# sct_markers.sh
start=`date +%s`

##################################################
# Move the RDS files from staging
ls
pwd
mkdir files
ls
cp /staging/ajwhitehead/vax_prr/h_lung/h_lung_combined.RDS.gz files/
cp /staging/ajwhitehead/vax_prr/m_lung/m_lung_combined.RDS.gz files/
# Unzip the files
gunzip files/*.gz

##################################################

# Copy the R script
cp /staging/ajwhitehead/vax_prr/SCT_Find_Markers.R ./

# Run the R script
Rscript SCT_Find_Markers.R

# Move the output files into staging 
cp mlung_markers.csv /staging/ajwhitehead/vax_prr/m_lung/
cp hlung_markers.csv /staging/ajwhitehead/vax_prr/h_lung/
ls

# Remove the large files so they don't get transferred back 
rm -r files
rm *.csv 
rm SCT_Find_Markers.R
ls
end=`date +%s`

runtime=$((end-start))
echo "Runtime was $runtime"
