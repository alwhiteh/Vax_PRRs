#!/bin/bash

# run_seurat.sh

##################################################
# Move the RDS files from staging
ls
pwd
mkdir files
ls
cp /staging/ajwhitehead/vax_prr/*.RDS files/
##################################################
# Copy the R script
cp /staging/ajwhitehead/vax_prr/make_h_lung_nostrip.R ./

# Run the R script
Rscript make_h_lung_nostrip.R

# Generate a directory for the output
mkdir /staging/ajwhitehead/vax_prr/h_lung

# Gzip the large RDS file that was outputted by the R script
gzip *.RDS

# Move the output files into staging 
cp h_lung_combined.RDS.gz /staging/ajwhitehead/vax_prr/h_lung/
cp hlung.combined.5markers.csv /staging/ajwhitehead/vax_prr/h_lung/
cp hlung.combined.markers.csv /staging/ajwhitehead/vax_prr/h_lung/
ls

# Remove the large files so they don't get transferred back 
rm -r files
rm h_lung_combined.RDS.gz
rm hlung.combined.5markers.csv
rm hlung.combined.markers.csv 
rm *.R
ls
