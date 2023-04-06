#!/bin/bash

# sct_markers.sh
start=`date +%s`

##################################################
# Move the RDS files from staging
ls
pwd
mkdir files
ls
cp /staging/ajwhitehead/vax_prr/m_lung_as_h/mlungh_ripped_sct.RDS.gz files/
# Unzip the files
gunzip files/*.gz

##################################################

# Copy the R script
cp /staging/ajwhitehead/vax_prr/SCT_and_Markers.R ./

# Run the R script
Rscript SCT_and_Markers.R

# Move the output files into staging 
cp files/mlungh_markers.csv /staging/ajwhitehead/vax_prr/m_lung_as_h/
cp files/mlungh_ripped_sct.RDS /staging/ajwhitehead/vax_prr/m_lung_as_h/
ls

# Remove the large files so they don't get transferred back 
rm -r files
rm *.csv 
rm SCT_and_Markers.R
ls
end=`date +%s`

runtime=$((end-start))
echo "Runtime was $runtime"
