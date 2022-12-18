library(Seurat)
# Read in the Seurat Objects
hlung.combined <- readRDS("files/h_lung_combined.RDS") #150k cells over 6 datasets
mlung.combined <- readRDS("files/m_lung_combined.RDS") 

# Process the human samples
DefaultAssay(hlung.combined) <- "SCT"
hlung.combined <- PrepSCTFindMarkers(hlung.combined)
hlung.markers <- FindAllMarkers(hlung.combined)
# Process mouse
DefaultAssay(mlung.combined) <- "SCT"
mlung.combined <- PrepSCTFindMarkers(mlung.combined)
mlung.markers <- FindAllMarkers(mlung.combined)

# Write markers to output
write.csv(hlung.markers, "hlung_markers.csv")
write.csv(mlung.markers, "mlung_markers.csv")

