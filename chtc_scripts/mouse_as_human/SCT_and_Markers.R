library(Seurat)
# Read in the Seurat Objects
hlung.combined <- readRDS("files/mlungh_ripped_sct.RDS")

# The active assay default is RNA
DefaultAssay(hlung.combined) <- "integrated"
all_genes <- rownames(hlung.combined)
hlung.combined <- ScaleData(hlung.combined, 
                            features = all_genes, 
                            verbose = FALSE)
hlung.combined <- RunPCA(hlung.combined, 
                         features = all_genes, 
                         npcs = 30)
hlung.combined <- RunUMAP(hlung.combined, 
                          reduction = "pca", 
                          dims = 1:30)
hlung.combined <- FindNeighbors(hlung.combined, 
                                reduction = "pca", 
                                dims = 1:30)
hlung.combined <- FindClusters(hlung.combined, 
                               resolution = 0.3)
hlung.combined <- SCTransform(hlung.combined, 
                              assay = "RNA",
                              verbose = FALSE)

# Find Markers
DefaultAssay(hlung.combined) <- "SCT"
hlung.combined <- PrepSCTFindMarkers(hlung.combined)
hlung.markers <- FindAllMarkers(hlung.combined)

# Write outputs
write.csv(hlung.markers, "files/mlungh_markers.csv")
saveRDS(hlung.combined, "files/mlungh_ripped_sct.RDS")
