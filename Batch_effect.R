##-----------------------------
## BATCH EFFECT CORRECTION: scRNA-seq Pipeline
##-----------------------------

# Set working directory
setwd("~/Documents/Drive G/Personal Project/scRNA_tutorial/GSE132771_RAW")

# Load Seurat
library(Seurat)

#-----------------------------
# Load data
#-----------------------------

NML1 <- Read10X(data.dir = "./NML_I")
NML2 <- Read10X(data.dir = "./NML_II")
NML3 <- Read10X(data.dir = "./NML_III")

# Create Seurat objects
NML1 <- CreateSeuratObject(counts = NML1, min.cells = 3, min.features = 200, project = "NML1")
NML2 <- CreateSeuratObject(counts = NML2, min.cells = 3, min.features = 200, project = "NML2")
NML3 <- CreateSeuratObject(counts = NML3, min.cells = 3, min.features = 200, project = "NML3")

# Merge datasets
Merged_NML <- merge(x = NML1, y = c(NML2, NML3), add.cell.ids = c("rep1", "rep2", "rep3"), project = "Merged_NML")

#-----------------------------
# Quality Control
#-----------------------------

# Calculate % mitochondrial genes
Merged_NML <- PercentageFeatureSet(Merged_NML, pattern = "^MT-", col.name = "percent.mt")

# Visualize QC metrics
VlnPlot(Merged_NML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells based on QC thresholds
Merged_NML <- subset(Merged_NML, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 &
                                   nCount_RNA < 30000 & percent.mt < 10)

# Visualize after filtering
VlnPlot(Merged_NML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Save QC-processed object
saveRDS(Merged_NML, file = "./Merged_NML.RDS")

#-----------------------------
# Normalization, Variable Features, Scaling
#-----------------------------

Merged_NML <- NormalizeData(Merged_NML)
Merged_NML <- FindVariableFeatures(Merged_NML, selection.method = "vst", nfeatures = 2000)
Merged_NML <- ScaleData(Merged_NML)

# Visualize variable features
top10 <- head(VariableFeatures(Merged_NML), 10)
VariableFeaturePlot(Merged_NML) + 
  LabelPoints(points = top10, repel = TRUE)

#-----------------------------
# Dimensionality Reduction
#-----------------------------

Merged_NML <- RunPCA(Merged_NML)
DimPlot(Merged_NML, reduction = "pca")
ElbowPlot(Merged_NML)

#-----------------------------
# Clustering & UMAP/tSNE
#-----------------------------

Merged_NML <- FindNeighbors(Merged_NML, dims = 1:20)
Merged_NML <- FindClusters(Merged_NML, resolution = 0.3)

Merged_NML <- RunUMAP(Merged_NML, dims = 1:20)
DimPlot(Merged_NML, reduction = "umap", label = TRUE, repel = TRUE)

Merged_NML <- RunTSNE(Merged_NML, dims = 1:20)
DimPlot(Merged_NML, reduction = "tsne", label = TRUE, repel = TRUE)

# Marker gene visualization
FeaturePlot(Merged_NML, features = c("EPCAM", "CLDN5", "COL1A2", "PTPRC"), cols = c('lightgrey', 'blue'))

# UMAP grouped by sample
DimPlot(Merged_NML, reduction = "umap", group.by = "orig.ident")

#-----------------------------
# Integration / Batch Correction
#-----------------------------

# Load pre-QC merged data
Merged_NML <- readRDS("./Merged_NML.RDS")

Merged_NML <- NormalizeData(Merged_NML)
Merged_NML <- FindVariableFeatures(Merged_NML)
Merged_NML <- ScaleData(Merged_NML)
Merged_NML <- RunPCA(Merged_NML)

# CCA integration
Merged_NML <- IntegrateLayers(Merged_NML, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# Re-join layers
Merged_NML[["RNA"]] <- JoinLayers(Merged_NML[["RNA"]])

# Clustering and UMAP on integrated data
Merged_NML <- FindNeighbors(Merged_NML, reduction = "integrated.cca", dims = 1:20)
Merged_NML <- FindClusters(Merged_NML, resolution = 0.3)
Merged_NML <- RunUMAP(Merged_NML, dims = 1:20, reduction = "integrated.cca")

# Plot UMAP before and after integration
p1 <- DimPlot(Merged_NML, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Merged_NML, reduction = "umap", group.by = "orig.ident")
p1 + p2


