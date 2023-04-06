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
# Copy the R script from staging
cp /staging/ajwhitehead/vax_prr/make_m_lung.R ./

# Run the R script
Rscript make_m_lung.R

# Generate a directory for the output
mkdir /staging/ajwhitehead/vax_prr/m_lung

# Gzip the large RDS file that was outputted by the R script
gzip *.RDS

# Move the output files into staging 
cp m_lung_combined.RDS.gz /staging/ajwhitehead/vax_prr/m_lung/
cp mlung.combined.5markers.csv /staging/ajwhitehead/vax_prr/m_lung/
cp mlung.combined.markers.csv /staging/ajwhitehead/vax_prr/m_lung/
ls

# Remove the large files so they don't get transferred back 
rm -r files
rm m_lung_combined.RDS.gz
rm mlung.combined.5markers.csv
rm mlung.combined.markers.csv 
ls
