## ### BATCH effect
setwd("~/Documents/Drive G/Personal Project/scRNA_tutorial/GSE132771_RAW")

# Loading seurat

library(Seurat)

#Loading the scRNA datafiles

NML1 <- Read10X(data.dir = "./NML_I") 
NML2 <- Read10X(data.dir = "./NML_II")
NML3 <- Read10X(data.dir = "./NML_III")

#to remove a file use the rm function
rm (IPF)
#creating a seurat object from the dgmatrix

NML1 <- CreateSeuratObject(counts = NML1, min.cells = 3, min.features = 200, project = "NML1")
NML2 <- CreateSeuratObject(counts = NML2, min.cells= 3, min.features = 200, project ="NML2")
NML3 <- CreateSeuratObject(counts = NML3, min.cells = 3, min.features = 200, project = "NML3")

#Merge all the three datasets

Merged_NML <- merge(x= NML1, y = c(NML2, NML3), add.cell.ids= c("rep1", "rep2", "rep3"), project = "Merged_NML")

# The next step after merging the dataset would be to perform quality control

# getting the percentage mitochondrial genes, name the col.name as percent.mt

Merged_NML <- PercentageFeatureSet(Merged_NML, pattern = "^MT-", col.name = "percent.mt")
View(Merged_NML@meta.data)

# visualize the variability in the merged data set among the replicates to start with data preprocessing
# we use violin plot to determine the cutoff values for QC

VlnPlot(Merged_NML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Quality control
#remove cells with low feature counts or too high feature counts, too much of mitochondrial genes

Merged_NML<- subset(Merged_NML, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & 
                      nCount_RNA <30000 & percent.mt < 10)

# view the data distribution again using violin plot

VlnPlot(Merged_NML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(object = Merged_NML, file = "./Merged_NML.RDS")

Merged_NML <- readRDS("./Merged_NML.RDS")
# Normalize the data

Merged_NML_QC <- NormalizeData(Merged_NML, normalization.method = "LogNormalize", scale.factor = 10000)

# the normalized count matrix is scored in pbmc[["RNA"]]$data
View(Merged_NML_QC[["RNA"]]$data)

#Find the features top 2000 genes with high variablity, this is done to reduce the dimentionality of the data
# The function FindVariableFetures stores these features in the VariableFeature slot in the seurat object
Merged_NML_QC

Merged_NML_QC <- FindVariableFeatures(Merged_NML_QC, selection.method = "vst", nfeatures = 2000)

Merged_NML_QC
# To view the VariableFeatures
VariableFeatures(Merged_NML_QC)
VariableFeaturePlot(Merged_NML_QC)

#plot gene variance vs gene expression
top10 <- head(VariableFeatures(Merged_NML_QC),10)
LabelPoints(plot = VariableFeaturePlot(Merged_NML_QC),points = top10, repel = TRUE)
# we could use these 2000 genes to perform the downstream analysis.

#Scale the data using ScaleData
Merged_NML_QC <- ScaleData(Merged_NML_QC) #2000 variable features.

View (Merged_NML_QC[["RNA"]]$scale.data)

# We can perform scaling using all genes if the dataset is small.
#all.genes <- rownames(Merged_NML_QC)
#Merged_NML_QC <- ScaleData(Merged_NML_QC, features = all.genes)
# Done with standrad pre-processing workfolw
#DImentioanlity reduction
# Principle component analysis: propcess of computing prinicple component, it's a linear dimentionality reduction
#technique. It projects each data point onto first few first few PC.
# By default it runs firts 50 PC

Merged_NML_QC <- RunPCA(Merged_NML_QC)
# Examine and visialize PCA using
# DimPlot and DimHeatmap
DimPlot(Merged_NML_QC, reduction = "pca") #+ NoLegend()
DimPlot(Merged_NML_QC, reduction = "pca", dim = c(1,5))
DimPlot(Merged_NML_QC, reduction = "pca", dim = c(2,3)) 

DimHeatmap(Merged_NML_QC, dims = 1, cells = 500, nfeatures = 30, balanced = TRUE)

# Two ways to determine the dimentionality of teh data set:
# Jackstraw function and Elbow Plot function
#Elbow lot to identify pc dims with significant variance

ElbowPlot(Merged_NML_QC)

#### Cell Clustering

# first calculate the distances between the pca vectors for each pair of cells
Merged_NML_QC<- FindNeighbors(Merged_NML_QC, dims = 1:20)
#The cluster the cells
Merged_NML_QC <- FindClusters(Merged_NML_QC, resolution = 0.3)

#### Non-linear dimentional reduction technique (UMAP/tSNE)

Merged_NML_QC <- RunUMAP( Merged_NML_QC, dims = 1:20) # use the same pc that we used for find neighbors

DimPlot(Merged_NML_QC, reduction = "umap", label = TRUE,repel = TRUE)

Merged_NML_QC <- RunTSNE(object = Merged_NML_QC)

DimPlot(Merged_NML_QC, reduction = "tsne", label = TRUE, repel = TRUE)

# Epithelial, Endothelial, mesnchymal, immune cells repsectively.
FeaturePlot(Merged_NML_QC, features = c("EPCAM", "CLDN5", "COL1A2", "PTPRC"), 
            cols = c('lightgrey','blue'))

#UMAP grouped by the orig.ident

plot1 <- DimPlot(Merged_NML_QC, reduction ="umap", group.by = "orig.ident")

# Reading in the RDS dataset without pre-processing step. QC done

Merged_NML <- readRDS("./Merged_NML.RDS")

## Data Integration
#ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
Merged_NML <- NormalizeData(Merged_NML)
Merged_NML <- FindVariableFeatures(Merged_NML)
Merged_NML <- ScaleData(Merged_NML)
Merged_NML <- RunPCA(Merged_NML)

Merged_NML <- IntegrateLayers(object = Merged_NML, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
Merged_NML[["RNA"]] <- JoinLayers(Merged_NML[["RNA"]])

Merged_NML <- FindNeighbors(Merged_NML, reduction = "integrated.cca", dims = 1:20)
Merged_NML <- FindClusters(Merged_NML, resolution = 0.3)
Merged_NML <- RunUMAP(Merged_NML, dims = 1:20, reduction = "integrated.cca")
plot2 <- DimPlot(Merged_NML, reduction ="umap", group.by = "orig.ident")
plot1 + plot2


#### CCAIntegration######

saveRDS(Merged_IPF, file = "./Merged_IPF.rds")
temp <- Merged_IPF
temp <- IntegrateLayers(object = temp, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)
temp[["RNA"]] <- JoinLayers(temp[["RNA"]])

temp <- FindNeighbors(temp, reduction = "integrated.cca", dims = 1:20)
temp <- FindClusters(temp, resolution = 0.3)
temp <- RunUMAP(temp, dims = 1:20, reduction = "integrated.cca", reduction.name = "umap.cca")
DimPlot(temp, reduction = "umap", group.by = "orig.ident")
Merged_IPF_ori<- readRDS("./Merged_IPF.rds")
temp<- readRDS("../Merged_IPF")
plot1<- DimPlot(Merged_IPF_ori, reduction ="umap", group.by = "orig.ident")
plot2<- DimPlot(temp, reduction ="umap", group.by = "orig.ident")
plot1 + plot2

Merged_IPF <- readRDS("./Merged_IPF.rds")
View(Merged_IPF@meta.data)
Merged_IPF[["RNA"]] <- split (Merged_IPF[["RNA"]], f = Merged_IPF$orig.ident)
# Cluster biomarkers
cluster2.markers <- FindMarkers(Merged_IPF, ident.1 = 2)
head(cluster2.markers, n = 5)
