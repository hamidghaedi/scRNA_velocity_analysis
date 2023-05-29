# Seurat + velocuto R + scvelo 

# libs (ran on CC)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

# Reading the data and doing pre-processing with Seurat
## example data 
curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile = 'SCG71.loom')
#The following step needs velocytoR to be installed  
ldat <- ReadVelocity(file = "SCG71.loom")
# convert object to Seurat obj
bm <- as.Seurat(x = ldat)

# pre-processing in Seurat obj
bm[["RNA"]] <- bm[["spliced"]]
bm <- SCTransform(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
DefaultAssay(bm) <- "RNA"
# saving and converting seurat obj inot other format
SaveH5Seurat(bm, filename = "mouseBM.h5Seurat")
Convert("mouseBM.h5Seurat", dest = "h5ad")