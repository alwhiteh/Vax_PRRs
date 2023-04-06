# Merge a human and 'humanized' mouse dataset
library(Seurat)
# Read in data
mouse <- readRDS("files/mlungh_anno.RDS")
mouse$organism <- "mouse"
human <- readRDS("files/hlung_anno.RDS")
human$organism <- "human"

## Generate Test data - delete later
# mouse <- readRDS("mouse_as_human/mlungh_tiny.RDS")
# human <- readRDS("human/hlung_tiny.RDS")
# mouse$organism <- "mouse"
# human$organism <- "human"


# DefaultAssay(mouse) <- "integrated"
# DefaultAssay(human) <- "integrated"
# lung_list <- list(mouse, human)
# lung_list <- lapply(X = lung_list, FUN = SCTransform)

# Set assay to SCT 
DefaultAssay(mouse) <- "SCT"
DefaultAssay(human) <- "SCT"

lung_list <- list(mouse, human)

# for (i in 1:length(lung_list)) {
#   lung_list[[i]] <- FindVariableFeatures(lung_list[[i]], 
#                                            selection.method = "vst", 
#                                            nfeatures = 2000,
#                                            verbose = FALSE)
# }

features <- SelectIntegrationFeatures(object.list = lung_list)

lung_list <- PrepSCTIntegration(object.list = lung_list, 
                                  anchor.features = features)



lung_list <- lapply(X = lung_list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find the integration features using rPCA (recommended for large datasets)
lung.anchors <- FindIntegrationAnchors(object.list = lung_list, 
                                        anchor.features = features, 
                                        reduction = "cca", 
                                        dims = 1:50, 
                                        k.anchor = 20, 
                                        normalization.method = "SCT")#,
                                       #reference = 2) #human ref

# create list of common genes to keep (and integrate)
to_integrate <- Reduce(c, lapply(lung.anchors@object.list, rownames))

# integrate data and keep full geneset
lung.combined <- IntegrateData(anchorset = lung.anchors, 
                                features.to.integrate = to_integrate,
                               normalization.method = "SCT",
                               dims = 1:50)

# The active assay default is RNA
DefaultAssay(lung.combined) <- "SCT"
all_genes <- rownames(lung.combined)
DefaultAssay(lung.combined) <- "integrated"
lung.combined <- ScaleData(lung.combined,
                           features = all_genes,
                           verbose = FALSE)
lung.combined <- RunPCA(lung.combined, 
                         features = all_genes, 
                         npcs = 30)
lung.combined <- RunUMAP(lung.combined, 
                          reduction = "pca", 
                          dims = 1:30)
lung.combined <- FindNeighbors(lung.combined, 
                                reduction = "pca", 
                                dims = 1:30)
lung.combined <- FindClusters(lung.combined, 
                               resolution = 0.3)
#lung.combined <- SCTransform(lung.combined, 
#                              assay = "RNA",
#                              verbose = FALSE)

# Find Markers
DefaultAssay(lung.combined) <- "SCT"
lung.combined <- PrepSCTFindMarkers(lung.combined)
lung.markers <- FindAllMarkers(lung.combined)

# Write outputs
write.csv(lung.markers, "files/lung_markers.csv")
saveRDS(lung.combined, "files/lung_sct.RDS")

