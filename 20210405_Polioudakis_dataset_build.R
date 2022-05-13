## Dataset build for fetal brain data from Polioudakis et al (2019)
## Joe Raymond

#install.packages('Seurat')
library(Seurat)

## Loading Data
# Prior to running this script, data must be downloaded from Polioudakis et al's
# 2019 publication: 10.1016/j.neuron.2019.06.011
load('/path/to/data/matrix')
metadata <- read.csv('path/to/metadata', row.names = 1)
obj.save.path <- '/where/to/save/seurat/object'
dataset <- CreateSeuratObject(counts = raw_counts_mat, meta.data = metadata)
# ensuring that all cells in metadata exist in data matrix and vice versa
meta.rows.to.keep <- intersect(rownames(dataset@meta.data), colnames(dataset))
dataset@meta.data <- dataset@meta.data[meta.rows.to.keep,]
rm(raw_counts_mat, metadata, meta.rows.to.keep)

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
dims <- 14
dataset <- FindNeighbors(dataset, reduction = 'pca',
                         assay = 'SCT', dims = 1:dims)
dataset <- FindClusters(dataset, resolution = .2, dims = 1:dims)
dataset <- RunUMAP(dataset, dims = 1:dims)
save(dataset, file = obj.save.path)
rm(dims)

# Simple function to improve readability of cell selection
JR.which.cells <- function(dataset = dataset,
                           meta.col = 'genotype',
                           which = 'WT'){
  return(rownames(dataset@meta.data[dataset@meta.data[,meta.col] == which,]))
}

## Cell type annotation
# For the sake of simplicity (and readability) we have collapsed some of the
# cell type categories into broader categories. The code below should make it 
# relatively easy to see how these broader categories are defined in relation
# to the original publication.
dataset$celltype <- as.character(dataset$Cluster)
meta <- dataset@meta.data
maturing.neurons <- JR.which.cells(dataset, meta.col = 'Cluster',
                                   which = 'ExM')
migrating.neurons <- JR.which.cells(dataset, meta.col = 'Cluster',
                                    which = 'ExN')
maturing.upper <- JR.which.cells(dataset, meta.col = 'Cluster',
                                 which = 'ExM-U')
oRG.cells <- JR.which.cells(dataset, meta.col = 'Cluster', which = 'oRG')
vRG.cells <- JR.which.cells(dataset, meta.col = 'Cluster', which = 'vRG')
mit.cells <- c(JR.which.cells(dataset, meta.col = 'Cluster', which = 'PgG2M'),
               JR.which.cells(dataset, meta.col = 'Cluster', which = 'PgS'))
inh.cells <- c(JR.which.cells(dataset, meta.col = 'Cluster', which = 'InCGE'),
               JR.which.cells(dataset, meta.col = 'Cluster', which = 'InMGE'))
dp.cells <- c(JR.which.cells(dataset, meta.col = 'Cluster', which = 'ExDp1'),
              JR.which.cells(dataset, meta.col = 'Cluster', which = 'ExDp2'))
dataset@meta.data[maturing.neurons,]$celltype <- 'Maturing Neurons'
dataset@meta.data[maturing.upper,]$celltype <- 'Maturing Upper Layer Neurons'
dataset@meta.data[migrating.neurons,]$celltype <- 'Migrating Neurons'
dataset@meta.data[oRG.cells,]$celltype <- 'Outer Radial Glia'
dataset@meta.data[vRG.cells,]$celltype <- 'Ventral Radial Glia'
dataset@meta.data[mit.cells,]$celltype <- 'Cycling Cells'
dataset@meta.data[inh.cells,]$celltype <- 'Inhibitory Neurons'
dataset@meta.data[dp.cells,]$celltype <- 'Deep Layer Excitatory Neurons'
dataset@meta.data[JR.which.cells(dataset,
            meta.col = 'Cluster', which = 'Mic'),]$celltype <- 'Microglia'
dataset@meta.data[JR.which.cells(dataset,
            meta.col = 'Cluster', which = 'End'),]$celltype <- 'Endothelial'
dataset@meta.data[JR.which.cells(dataset,
            meta.col = 'Cluster', which = 'Per'),]$celltype <- 'Pericytes'
dataset@meta.data[JR.which.cells(dataset,
            meta.col = 'Cluster', which = 'OPC'),]$celltype <- 'OPCs'
dataset@meta.data[JR.which.cells(dataset,
    meta.col = 'Cluster', which = 'IP'),]$celltype <- 'Intermediate Progenitors'
save(dataset, file = obj.save.path)