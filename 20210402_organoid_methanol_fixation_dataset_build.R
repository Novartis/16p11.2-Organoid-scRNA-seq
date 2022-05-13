## Dataset build for scRNA-seq analysis of organoid methanol fixation experiment
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
paths <- paste0(path.to.h5s, samples)
obj.save.path <- 'where_you_want_to_save_seurat_object'
project.name <- 'Organoids'
# Reads each file and create seurat objects to merge
dataset.data <- Read10X_h5(paths[1])
dataset.x <- CreateSeuratObject(dataset.data, project = project.name)
dataset.data <- Read10X_h5(paths[2])
dataset.y <- CreateSeuratObject(dataset.data, project = project.name)
dataset <- merge(x = dataset.x, y = dataset.y, add.cell.ids = samp.names)
rm(path.to.h5s, samples, samp.names, paths, project.name,
   dataset.data, dataset.x, dataset.y)

## Filtering cells, setting percent mitochondrial genes, and adding metadata
dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset <- subset(dataset,
            subset = nFeature_RNA < 5000 & nCount_RNA < 30000 & percent.mt < 5)
# simple function to name samples
JR.name.samples <- function(dataset,
                            samples = c("sample1","sample2","sample3")){
  require(Seurat)
  lane.ids <- as.numeric(gsub("[^0-9]","",colnames(dataset)))
  condition <- as.character(factor(lane.ids, labels = samples))
  sample.labels <- data.frame(row.names = colnames(dataset), lane.ids, condition)
  dataset <- AddMetaData(dataset, sample.labels)
  return(dataset)
}
dataset <- JR.name.samples(dataset, samples = c('Fresh', 'Fixed'))


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
dims <- 14
dataset <- FindNeighbors(dataset, reduction = 'pca', assay = 'SCT', dims = 1:dims)
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
# In addition to transferring the cell type annotations from the fetal brain
# data, we integrated the datasets using Seurat's integration workflow. The
# workflow used to perform this integration in shown below. The UMAP from this
# integration was used in lieu of the UMAP produced from the organoid data by 
# itself for the sake of keeping things consistent between different figures.

## Experiment specific variables
# These names will be stored as metadata category "exp" for use in
# downstream analysis. The variable reference_dataset refers to which 
# sample is the reference that other datasets are being mapped to.
dataset1.cellids <- "Organoids"
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
## Transferring cell type annotations to integrated data
cts <- dataset$predicted.id #list of cell types for organoids
dataset.int$predicted.id <- dataset.int$celltype #original cell type annotations
dataset.int@meta.data[colnames(dataset),]$predicted.id <- cts
save(dataset.int, file = obj.save.path.int)
rm(cts)



