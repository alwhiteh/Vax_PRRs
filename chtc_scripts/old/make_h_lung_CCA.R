# Load the packages
getwd()
library(dplyr)
library(stringr)
library(Seurat)

#Load the RDS files
h_lung_endo_tata <- readRDS("files/h_lung_endo_tata.RDS") #4872
h_lung_epi_tat <- readRDS("files/h_lung_epi_tat.RDS") #4621
h_lung_haber <- readRDS("files/h_lung_haber.RDS") #4621
h_lung_imm_tata <- readRDS("files/h_lung_immune_tata.RDS") #4621
h_lung_kras <- readRDS("files/h_lung_kras.RDS") #E4344 # Use as ref
h_lung_okuda <- readRDS("files/h_lung_okuda.RDS") #3747
h_lung_strip <- readRDS("files/h_lung_strip.RDS") #3747
h_nasal_teich <- readRDS("files/h_nasal_teich.RDS")

#We need to filter out the pathological samples from haberman
h_lung_haber$celltype <- Idents(h_lung_haber)
Idents(h_lung_haber) <- "Diagnosis"
health_haber_cells <- WhichCells(h_lung_haber, idents = "Control")
h_lung_haber <- subset(h_lung_haber, cells =  health_haber_cells)

m_lung_list <- list(h_lung_endo_tata,
                    h_lung_epi_tat,
                    h_lung_haber,
                    h_lung_imm_tata,
                    h_lung_kras, 
                    h_lung_okuda,
                    h_lung_strip,
                    h_nasal_teich)

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
                                        reduction = "cca", 
                                        dims = 1:30, 
                                        k.anchor = 20, 
                                        normalization.method = "SCT",
					reference = 5)

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

write.csv(mlung.combined.markers,"hlung.combined.markers.csv")
mlung.combined.markers %>% 
  group_by(cluster) %>% 
  slice_head(n = 5) -> Top5Markers
write.csv(Top5Markers,"hlung.combined.5markers.csv")

# Save the combined dataset as an RDS file
saveRDS(mlung.combined, "h_lung_combined.RDS")
