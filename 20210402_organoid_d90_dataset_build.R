## Dataset build for scRNA-seq analysis of d90 organoids
## Joe Raymond

#install.packages('Seurat')
library(Seurat)

# Helper function to run scds on dataset
# mostly just splits the object by sample and runs the analysis on each
# individual sample before recombining metadata
JR.scds.helper <- function(dataset = dataset, sample.meta.col.name = 'sample'){
  require(scds)
  require(Seurat)
  require(scater)
  datasets.list <- SplitObject(dataset, split.by = sample.meta.col.name)
  meta <- c()
  for(i in datasets.list){
    sce <- as.SingleCellExperiment(i)
    sce <- cxds_bcds_hybrid(sce) #main function
    if(length(meta > 0)){
      meta <- rbind(meta, as.data.frame(colData(sce)))
    } else{
      meta <- as.data.frame(colData(sce))
    }
  }
  dataset@meta.data <- meta
  return(dataset)
}

## Data Loading
# Data is loaded from the CellRanger output into a Seurat object.
# This chunk of code will only function properly if data has been downloaded
# with its original filenames and stored by itself in a directory (path.to.h5s).
path.to.h5s <- 'your_local_path_to_datasets'
samples <- list.files(path.to.h5s)
samp.names <- gsub(x = samples, pattern = '.h5', replacement = '')
paths <- list.files(path.to.h5s, full.names = T)
obj.save.path <- 'where_you_want_to_save_seurat_object'
dataset.data <- Read10X_h5(paths[1])
dataset.x <- CreateSeuratObject(dataset.data, project = 'd90_All')
dataset.list <- c()
# iterate through list of samples 2-26 for input for data merge
for(i in paths[2:length(paths)]){
  dataset.data <- Read10X_h5(i)
  dataset <- CreateSeuratObject(dataset.data, project = 'd90_all')
  dataset.list <- c(dataset.list, dataset)
}
dataset <- merge(x = dataset.x, y = dataset.list, add.cell.ids = samp.names)
rm(dataset.data, paths, start.path, end.path, i,
   dataset.x, dataset.list, samp.ids)

## Adding some metadata
# Basic metadata is added prior to further analysis to ensure variables
# that might need to be regressed can be easily addressed.
meta <- dataset@meta.data
meta <- cbind(meta, 'CB' = rownames(meta))
meta <- separate(data = meta, col = 'CB',
        into = c('cell_line', 'genotype', 'age', 'batch', 'experiment', 'CB'),
        sep = '_')
meta$sample <- paste(meta$cell_line, meta$genotype, meta$age,
                     meta$batch, meta$experiment, sep = "_")
meta$type <- 'iPSC'
meta[meta$cell_line == 'H1',]$type <- 'ESC'
meta[meta$cell_line == 'H9',]$type <- 'ESC'
dataset@meta.data <- meta
rm(meta)

## Filtering out cells
# Cells with an abnormally high number of reads or
# genes are filtered out as suspected doublets.
summary(dataset@meta.data$nCount_RNA > 15000)
#setting filter at 15,000 removes 3114 cells
dataset <- subset(dataset, subset = nCount_RNA < 15000)

##Identifying Doublets
# Using [SCDS](https://github.com/kostkalab/scds) to detect doublets
# We expect .8% doublets per 1000 cells and will use this as a basis to choose
# the correct score to filter by.
dataset <- JR.scds.helper(dataset, sample.meta.col.name = 'sample')
dim(dataset)
summary(as.factor(dataset$sample))
#doublet rate expected based on 10X documentation
doublet.rate.exp <- (summary(as.factor(dataset$sample))/1000*.8)*summary(as.factor(dataset$sample))/100
#df will contain data frame of the number of doublets detected 
# per sample at a given hybrid score
df <- data.frame(row.names = unique(dataset$sample))
df.names <- c()
for(i in seq(.1, 2, by = .1)){
  name <- paste0('hybridscore_', i)
  col <- c()
  for(j in rownames(df)){
    meta <- dataset@meta.data
    meta <- meta[meta$sample == j,]
    x <- length(rownames(meta[meta$hybrid_score > i,]))
    col <- c(col, x)
  }
  df <- cbind(df, col)
  df.names <- c(df.names, name)
}
colnames(df) <- df.names
norm <- df-doublet.rate.exp #subtracting expected doublets from found
# doublets at each score
#looking for lowest colMeans and colSums in norm to show number of 
# doublets detected is closest to the expected number from 10X
colMeans(norm)
colSums(norm)
dataset <- subset(dataset, subset = hybrid_score < 1.4)
dim(dataset)
rm(df, df.names, norm, x, i , j, col, name, meta, doublet.rate.exp)

## Normalization, scaling, and dimensionality reduction
# Data is normalized and scaled using SCTransform and principal component 
# analysis is run on the scaled/normalized data.
options(future.globals.maxSize= 10*1000*1024^2) #10GB
dataset <- SCTransform(dataset, assay = 'RNA',
                    vars.to.regress = c('nFeature_RNA', 'batch', 'experiment'),
                    verbose = F)
dataset <- RunPCA(dataset, assay = 'SCT')
ElbowPlot(dataset)

## Clustering and 2D visualization
# We determine k-nearest neighbors for each cell and use this graph to calculate
# a shared nearest neighbors (SNN) graph. A smart local moving algorithm uses
# this graph to identify clusters. Uniform Manifold Approximation and Projection
# (UMAP) is used to reduce the high dimensional PC data into a 2D visualization.
dims <- 15 
dataset <- FindNeighbors(dataset, reduction = 'pca', assay = 'SCT',
                         dims = 1:dims)
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

