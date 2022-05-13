## Dataset build for scRNA-seq analysis of organoid timecourse experiment
## Joe Raymond

#install.packages('Seurat')
library(Seurat)

## Data Loading
# Data is loaded from the CellRanger output into a Seurat object.
# This chunk of code will only function properly if data has been downloaded
# with its original filenames and stored by itself in a directory (path.to.h5s).
path.to.h5s <- 'your_local_path_to_datasets'
samples <- list.files(path.to.h5s)
samp.names <- gsub(x = samples, pattern = '.h5', replacement = '')
paths <- list.files(path.to.h5s, full.names = T)
obj.save.path <- 'where_you_want_to_save_seurat_object'
project.name <- 'Organoids'
# Reads through each of the files and creates and merges Seurat objects
dataset.data <- Read10X_h5(paths[1])
dataset.x <- CreateSeuratObject(dataset.data, project = project.name)
dataset.list <- c()
for(i in paths[2:length(paths)]){
  dataset.data <- Read10X(i)
  dataset <- CreateSeuratObject(dataset.data, project = project.name)
  dataset.list <- c(dataset.list, dataset)
}
dataset <- merge(x = dataset.x, y = dataset.list, add.cell.ids = samp.names)
rm(path.to.h5s, samples, samp.names, paths, project.name,
   dataset.data, dataset.x, dataset.list)

## Adding metadata
# pulling the age information from the cell names
d166 <- colnames(dataset)[grepl(pattern = '_166', x = colnames(dataset))]
d90 <- colnames(dataset)[grepl(pattern = '_90', x = colnames(dataset))]
d30 <- colnames(dataset)[grepl(pattern = '_30', x = colnames(dataset))]
dataset$age <- '0'
dataset@meta.data[d166,]$age <- '166'
dataset@meta.data[d90,]$age <- '90'
dataset@meta.data[d30,]$age <- '30'
# calculating percent of mitochondrial genes captured
dataset$percent.mt <- PercentageFeatureSet(dataset, pattern = "^MT-")

## Filtering Cells
#Data was filtered to exclude cells with more than 6000 
#counts per cell or more than 10% mitochondrial reads.
dataset <- subset(dataset,
            subset = nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA < 35000)

## Normalization, scaling, and dimensionality reduction
# Data is normalized and scaled using SCTransform and principal component 
# analysis is run on the scaled/normalized data.
dataset <- SCTransform(dataset, assay = 'RNA', verbose = F)
dataset <- RunPCA(dataset, assay = 'SCT')
ElbowPlot(dataset) # used to help select number of PCs to use

## Clustering and 2D visualization
# We determine k-nearest neighbors for each cell and use this graph to calculate
# a shared nearest neighbors (SNN) graph. A smart local moving algorithm uses
# this graph to identify clusters. Uniform Manifold Approximation and Projection
# (UMAP) is used to reduce the high dimensional PC data into a 2D visualization.
dims <- 12
dataset <- FindNeighbors(dataset, reduction = 'pca',
                         assay = 'SCT', dims = 1:dims)
dataset <- FindClusters(dataset, resolution = .2, dims = 1:dims)
dataset <- RunUMAP(dataset, dims = 1:dims)
save(dataset, file = obj.save.path)
rm(dims)

## Transferring simplified cell type annotations from Polioudakis et al
# The protocol below was used to annotate the organoid cells using the 
# transcriptional signatures from an annotated fetal brain data. This process is
# computationally intensive and the metadata gained by running this code
# can be imported directly from the metadata provided in the data submission.
dataset.ref <- load('path/to/Polioudakis/fetalbrain/data')
Sys.setenv('R_MAX_VSIZE'=30*1000*1024^2) #requesting 30Gb memory for calculation
anchors <- FindTransferAnchors(reference = dataset.ref, query = dataset,
                               normalization.method = 'SCT', verbose = T)
predictions <- TransferData(anchorset = anchors, refdata = dataset.ref$celltype)
dataset <- AddMetaData(object = dataset, metadata = predictions)
rm(dataset.ref, cells.to.keep)
save(dataset, file = obj.save.path)

### Integration of Polioudakis dataset with organoid dataset ###
# Prior to the addition of the anchor transfer workflow within Seurat, we
# simply integrated our data with our reference set and used the nearest
# neighbor graph to assign cell types. We kept the UMAP plot from this analysis
# for publication because we thought it did a better job projecting the 
# natural developmental trajectory of the brain than the projection obtained
# when looking only at the organoid data. For this reason, we have provided the
# code for the integration of the data, which is ultimately only used for the
# UMAP projection in the publication.

## Experiment specific variables
# These names will be stored as metadata category "exp" for use in
# downstream analysis. The variable reference_dataset refers to which 
# sample is the reference that other datasets are being mapped to.
dataset1.cellids <- "Timecourse_Organoids"
dataset2.cellids <- "Fetal_Brain"
reference_dataset <- 2
dataset1 <- load('/path/to/organoids')
dataset2 <- load('/path/to/braindata')
obj.save.path.int <- '/save/path/for/data'
## Integration 
# Anchor features are found using CCA between each cell and every other cell not
# from the same dataset. These anchors are ranked based on variation across 
# datasets and the least variable anchors are used to integrate data.
dataset1$exp <- dataset1.cellids
dataset2$exp <- dataset2.cellids
dataset.list <- list(dataset1, dataset2)
rm(dataset1, dataset2, dataset)
dataset.features <- SelectIntegrationFeatures(object.list = dataset.list,
                                              nfeatures = 3000)
options(future.globals.maxSize= 10*1000*1024^2) #10GB
dataset.list <- PrepSCTIntegration(object.list = dataset.list,
                                   anchor.features = dataset.features, 
                                   verbose = FALSE)
dataset.anchors <- FindIntegrationAnchors(object.list = dataset.list,
                                          normalization.method = "SCT",
                                          anchor.features = dataset.features,
                                          verbose = FALSE,
                                          reference = reference_dataset)
dataset.integrated <- IntegrateData(anchorset = dataset.anchors,
                                    normalization.method = "SCT",
                                    verbose = FALSE)
dataset.integrated$original.clusters <- dataset.integrated$seurat_clusters
dataset <- dataset.integrated
rm(dataset.list, dataset.anchors, dataset.features, dataset.integrated)
## Running analysis on integrated data 
# Integrated data goes through standard scRNA-seq analysis after integration.
dataset <- RunPCA(dataset, verbose = FALSE)
ElbowPlot(dataset)
dims <- 14
dataset <- FindNeighbors(dataset, dims = 1:dims)
dataset <- FindClusters(dataset, resolution = .2)
dataset <- RunUMAP(dataset, dims = 1:dims)
save(dataset, file = obj.save.path.int)
## Transferring projection to original organoid dataset
dataset.int <- dataset
dataset <- load(obj.save.path)
emb <- dataset.int@reductions$umap@cell.embeddings
emb <- emb[colnames(dataset),]
dataset@reductions$umap@cell.embeddings <- as.array(emb)
save(dataset, file = obj.save.path)
rm(emb, dataset.int)
