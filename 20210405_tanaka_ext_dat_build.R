## Dataset build for scRNA-seq analysis of select external publications
## Joe Raymond

#install.packages('Seurat')
library(Seurat)

# Prior to running this script, data must be downloaded from Tanaka et al's 
# 2020 publication: https://doi.org/10.1016/j.celrep.2020.01.038, which contains
# three datasets of interest as well as Trujillo et al's 2019 publication:
# https://doi.org/10.1016/j.stem.2019.08.002
path.dat <- '/path/to/Tanaka/data/matrix'
path.meta <- '/path/to/Tanaka/metadata'
path.truj <- '/path/to/Trujillo/data/matrix'
obj.save.path <- '/where/to/save/seurat/object'
meta <- read.delim(file = path.meta, header = T)
unique(meta$Age)
#We want Velasco, Birey, Quadrato, and Trujillo (not included)
meta <- dplyr::filter(meta, Age == '3m')
unique(meta$Dataset)
dat.to.keep <- c("Birey 3m", "Quadrato 3m rep1", "Quadrato 3m rep2",
                 "Velasco 3m rep1", "Velasco 3m rep2")
meta <- dplyr::filter(meta, Dataset %in% dat.to.keep == T)
dataset.data <- readMM(paste0(path.dat, 'matrix.mtx'))
cols <- read.delim(paste0(path.dat, 'barcodes.tsv'), header = F)
rows <- read.delim(paste0(path.dat, 'features.tsv'), header =F)
rownames(dataset.data) <- rows[,1]
colnames(dataset.data) <- cols[,1]
dataset.data <- dataset.data[,meta$cellId]
dataset <- CreateSeuratObject(counts = dataset.data, project = 'Tanaka')
dataset@meta.data <- meta
rownames(dataset@meta.data) <- colnames(dataset)
#Adding missing Trujillo data
dataset.data <- Read10X(data.dir = path.truj)
dataset.truj <- CreateSeuratObject(counts = dataset.data, project = 'Trujillo')
dataset.truj$Dataset <- 'Trujillo 3m'
dataset.truj$Protocol <- 'Trujillo'
dataset.truj$Age <- '3m'
dataset.truj$cellID <- colnames(dataset.truj)
dataset.all <- merge(dataset, dataset.truj)
dataset <- dataset.all
save(dataset, file = obj.save.path)
rm(dataset.all, dataset.data, dataset.truj, cols, rows, meta)

## Normalization, scaling, and dimensionality reduction
# Data is normalized and scaled using SCTransform and principal component 
# analysis is run on the scaled/normalized data.
dataset <- SCTransform(dataset, assay = 'RNA', verbose = F)
dataset <- RunPCA(dataset, assay = 'SCT')
ElbowPlot(dataset)

## Clustering and 2D visualization
# We determine k-nearest neighbors for each cell and use this graph to calculate
# a shared nearest neighbors (SNN) graph. A smart local moving algorithm uses
# this graph to identify clusters. Uniform Manifold Approximation and Projection
# (UMAP) is used to reduce the high dimensional PC data into a 2D visualization.
dims <- 10
dataset <- FindNeighbors(dataset, reduction = 'pca',
                         assay = 'SCT', dims = 1:dims)
dataset <- FindClusters(dataset, resolution = .2, dims = 1:dims)
dataset <- RunUMAP(dataset, dims = 1:dims)
save(dataset, file = obj.save.path)
rm(dims)

## Transferring Anchors for each dataset individually
# Using the same protocol as in our internal datasets to identify putative
# cell types based on fetal brain data.
datasets <- SplitObject(dataset, split.by = 'Protocol')
load('path/to/Polioudakis/fetalbrain/data')
dataset.ref <- dataset
meta <- c()
# iterating through each dataset
for(i in datasets){
  dataset <- i
  Sys.setenv('R_MAX_VSIZE'=30*1000*1024^2)
  anchors <- FindTransferAnchors(reference = dataset.ref, query = dataset,
                                 normalization.method = 'SCT', verbose = T)
  predictions <- TransferData(anchorset = anchors,
                              refdata = dataset.ref$celltype)
  dataset <- AddMetaData(object = dataset, metadata = predictions)
  meta <- rbind(meta, dataset@meta.data)
}
load(obj.save.path)
dataset@meta.data <- meta
save(dataset, file = obj.save.path)
rm(dataset.ref, datasets, i, meta, predictions, cells.to.keep, anchors)