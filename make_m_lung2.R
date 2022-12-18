# Load the packages
getwd()
library(dplyr)
library(stringr)
library(Seurat)

#Load the RDS files
m_lung_angel <- readRDS("files/m_lung_angelidis.RDS") #4872
m_lung_sam2 <- readRDS("files/m_lung_sam2.RDS") #4621
m_lung_sam3 <- readRDS("files/m_lung_sam3.RDS") #4621
m_lung_sam4 <- readRDS("files/m_lung_sam4.RDS") #4621
m_lung_kras <- readRDS("files/m_lung_kras.RDS") #E4344
m_lung_rare_f <- readRDS("files/m_lung_rare_f.RDS") #3747
m_lung_rare_m <- readRDS("files/m_lung_rare_m.RDS") #3747
m_lung_koenitzer1 <- readRDS("files/m_lung_koenitzer1.RDS")
m_lung_koenitzer2 <- readRDS("files/m_lung_koenitzer2.RDS")
m_lung_koenitzer3 <- readRDS("files/m_lung_koenitzer3.RDS")
m_lung_list <- list(m_lung_angel,
                    m_lung_sam2,
                    m_lung_sam3,
                    m_lung_sam4,
                    m_lung_kras, 
                    m_lung_rare_f,
                    m_lung_rare_m,
                    m_lung_koenitzer1,
                    m_lung_koenitzer2,
                    m_lung_koenitzer3)


# remove the low RNA count cells, high mitochondrial (dead) cells
#DefaultAssay(mlung.combined) <- "RNA"

for (i in 1:length(m_lung_list)) {
  m_lung_list[[i]] <- PercentageFeatureSet(m_lung_list[[i]], 
                                           pattern = "^MT-",
                                           col.name = "percent.mt")
  m_lung_list[[i]] <- subset(m_lung_list[[i]], 
                             subset = nFeature_RNA > 200  & percent.mt < 18)
}

# Normalize the data and find variable features for each dataset
for (i in 1:length(m_lung_list)) {
  m_lung_list[[i]] <- SCTransform(m_lung_list[[i]], 
                                  method = "glmGamPoi", 
                                  verbose = FALSE, 
                                  vars.to.regress = "percent.mt")
  m_lung_list[[i]] <- FindVariableFeatures(m_lung_list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = 2000,
                                           verbose = FALSE)
}
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = m_lung_list)
m_lung_list <- PrepSCTIntegration(object.list = m_lung_list, anchor.features = features)
m_lung_list <- lapply(X = m_lung_list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# Find the integration features using rPCA (recommended for large datasets)
mlung.anchors <- FindIntegrationAnchors(object.list = m_lung_list, 
                                        anchor.features = features, 
                                        reduction = "rpca", 
                                        dims = 1:30, 
                                        k.anchor = 20, 
                                        normalization.method = "SCT",
                                        reference = 1)

# create list of common genes to keep (and integrate)
to_integrate <- Reduce(c, lapply(mlung.anchors@object.list, rownames))

# integrate data and keep full geneset
mlung.combined <- IntegrateData(anchorset = mlung.anchors, 
                                features.to.integrate = to_integrate)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(mlung.combined) <- "integrated"
mlung.combined <- ScaleData(mlung.combined, verbose = FALSE)
mlung.combined <- RunPCA(mlung.combined, npcs = 30, verbose = FALSE)
mlung.combined <- RunUMAP(mlung.combined, reduction = "pca", dims = 1:30)
mlung.combined <- FindNeighbors(mlung.combined, reduction = "pca", dims = 1:30)
mlung.combined <- FindClusters(mlung.combined, resolution = 0.3)

#Find Markers for each cluster
DefaultAssay(mlung.combined) <- "SCT"
mlung.combined <- PrepSCTFindMarkers(mlung.combined)
mlung.combined.markers <- FindAllMarkers(mlung.combined, # this takes a while
                                         only.pos = TRUE, 
                                         min.pct = 0.25, 
                                         logfc.threshold = 0.25) 

write.csv(mlung.combined.markers,"mlung.combined.markers.csv")
mlung.combined.markers %>% 
  group_by(cluster) %>% 
  slice_head(n = 5) -> Top5Markers
write.csv(Top5Markers,"mlung.combined.5markers.csv")

# Save the combined dataset as an RDS file
saveRDS(mlung.combined, "m_lung_combined.RDS")
