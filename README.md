# scRNA_velocity_analysis
scRNA velocity analysis using Seurat, velocyto and scvelo 

The goal is  to gain insights into the dynamics of gene expression changes in individual cells over time. This analysis provides information about the rate and direction of changes in gene expression, allowing researchers to infer the future states of individual cells and identify cellular transitions and trajectories.

ScRNA velocity analysis relies on the concept of RNA velocity, which is the measurement of the ratio between unspliced and spliced mRNA molecules within a cell. By quantifying the abundance of these RNA species, it is possible to infer the directionality and speed of gene expression changes, providing valuable information about cell fate decisions, differentiation processes, and developmental trajectories.

## 1 ) Running `velocyto` on cellranger outputs

The command is simple and is like :

` velocyto run10x -@ 16 --samtools-memory 80000 -m /home/ghaedi/projects/def-gooding-ab/ghaedi/sc/hg38_rmsk.gtf $d /home/ghaedi/projects/def-gooding-ab/ghaedi/sc/refdata-gex-GRCh38-2020-A/genes/genes.gtf` 

It relies on `$d/out` subdirectory in cellranger outputs. Also it requires bam file sorting using samtools. However it is possible to run `velocyto` to get the `loom` outputs - RNA velocity outputs- but it did not throw an error and generated no output!It seemed there is an error with samtools and it turns out that the program wont generate an error message when something is wrong with samtools. So let's sort the bam files beforehand and then proceed with the `velocyto` :

```bash
#!/bin/bash
#SBATCH --account=#
#SBATCH --job-name=samtools_sort
#SBATCH --qos=privileged
#SBATCH --nodes=1                
#SBATCH --tasks-per-node=16       
#SBATCH --mem 80g
#SBATCH --time 72:00:00
#SBATCH --output=samtool_sort.%J.out
#SBATCH --error=samtool_sort.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#


module load samtools


base_directory="/home/ghaedi/projects/def-gooding-ab/ghaedi/sc/raw/cellranger_outs"

directories=(
  "SRR12603789"
  "SRR12603788"
  "SRR12603786"
  "SRR12603790"
  "SRR12603787"
  "SRR12603780"
  "SRR12603781"
  "SRR12603782"
  "SRR12603783"
  "SRR12603784"
  "SRR12603785"
)

# loop
for dir in "${directories[@]}"
do
  echo $dir
  input_bam="$base_directory/$dir/outs/possorted_genome_bam.bam"
  output_bam="$base_directory/$dir/outs/cellsorted_possorted_genome_bam.bam"

  samtools sort -l 7 -m 7647M -t CB -O BAM  -o "$output_bam" "$input_bam"
done
```
Now we can  continue with `velocyto` : 

```shell
velocyto run10x -@ 16 --samtools-memory 80000 -m /home/ghaedi/projects/def-gooding-ab/ghaedi/sc/hg38_rmsk.gtf $d /home/ghaedi/projects/def-gooding-ab/ghaedi/sc/refdata-gex-GRCh38-2020-A/genes/genes.gtf
```

## 2) pre-processing in R 

```r
# Seurat + velocuto R + scvelo 

# libs (ran on CC)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(harmony)
library(dplyr)
library(ggplot2)
library(SCP)

# preparing epi-seurat cell ids to be used with loom files
# subset su file to only include epi cells:

## on CC
epi_seurat <- readRDS("~/projects/def-gooding-ab/ghaedi/sc/epi_seurat.RDS")

## on local
epi_seurat <- readRDS("~/scRNA/github/epi_seurat.RDS")

# defining cells 
epi_cells <- rownames(epi_seurat@meta.data)

# Renaming idents
epi_seurat <- RenameIdents(object = epi_seurat, 
                           "0" = "basal_cell",
                           "1" = "cancer_associated_luminal_cell",
                           "2" = "differentiated_luminal_cell",
                           "3" = "unique_luminal_cell",
                           "4" = "immunomodulatory_luminal_cell",
                           "5" = "adhesion_signaling_luminal_cell")


# changing cell names to match with epi_cells
substr(epi_cells,12,12) <- ":"
epi_cells <- substr(epi_cells, 1, nchar(epi_cells)-2)
epi_cells <- paste0(epi_cells, "x")

# adding modified name to the metadatad
epi_seurat@meta.data$cells <- epi_cells
epi_seurat$clusters <- Idents(epi_seurat)

# saving metadata 
epiMet <- epi_seurat@meta.data
epiMet <- epiMet[, colnames(epiMet) %in% c("cells", "gender", "age", "Grade", "Invasiveness", "clusters")]
rownames(epiMet) <- epiMet$cells

#______ read files into Seurat objects on CC

# Reading loom files
base_dir <- "/home/ghaedi/projects/def-gooding-ab/ghaedi/sc/raw/cellranger_outs/"
# create list of samples
samples <- list.files(base_dir)


# read files inot Seurat objects on CC
for (sample in samples){
  print(paste0(sample))
  loom_file <- paste0(base_dir,sample, "/velocyto/", sample, ".loom")
  if (file.exists(loom_file)) {
    ldat <- ReadVelocity(file = loom_file)
    # convert object to Seurat obj
    su <- as.Seurat(x = ldat)
    # define a cells column 
    su$cells <- rownames(su@meta.data)
    # Subset seurat obj to include only epi cells
    su <- subset(su, subset = cells %in% epi_cells)
    #
    suMet <- su@meta.data
    # merging meta data
    suMet <- dplyr::left_join(suMet, epiMet, by = 'cells')
    rownames(suMet) <- suMet$cells
    #
    su@meta.data <- suMet
    # pre-processing in Seurat obj
    su[["RNA"]] <- su[["spliced"]]
    # 
    assign(sample, su)
    } else {
    print(paste0("Skipping file: ", sample))
  }
}

#______ read files into Seurat objects on local machine
# Reading loom files
base_dir <- "C:/Users/qaedi/OneDrive - Queen's University/Documents/scRNA_velocity/loom_files"
# create list of samples
samples <- list.files(base_dir)


# read files inot Seurat objects on CC
for (sample in samples){
  print(paste0(sample))
  loom_file <- sample
  if (file.exists(loom_file)) {
    ldat <- ReadVelocity(file = loom_file)
    # convert object to Seurat obj
    su <- as.Seurat(x = ldat)
    # define a cells column 
    su$cells <- rownames(su@meta.data)
    # Subset seurat obj to include only epi cells
    su <- subset(su, subset = cells %in% epi_cells)
    #
    suMet <- su@meta.data
    # merging meta data
    suMet <- dplyr::left_join(suMet, epiMet, by = 'cells')
    rownames(suMet) <- suMet$cells
    #
    su@meta.data <- suMet
    su[["RNA"]] <- su[["spliced"]]
    assign(sample, su)
  } else {
    print(paste0("Skipping file: ", sample))
  }
}



# now merging all objects into one Seurat obj

merged_seurat <- merge(x = SRR12603782,
                       y = c(SRR12603783,
                             SRR12603784,
                             SRR12603785))

# setting active assay
#DefaultAssay(merged_seurat) <- "RNA"
# adding sample name
merged_seurat$orig.ident <- substr(merged_seurat$cells,1,11)

# Integration: Preprocessing
# Perform log-normalization and feature selection, as well as SCT normalization on global object
merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData() %>%
  SCTransform(vars.to.regress = c("orig.ident"))


# Calculate PCs using variable features determined by SCTransform (3000 by default)
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)


# Integration
harmonized_seurat <- RunHarmony(merged_seurat, 
                                group.by.vars = c("orig.ident", "gender"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:40)
# to set reduction to harmony and finding the clusters
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8))
#
DefaultAssay(harmonized_seurat) <- "RNA"


# saving and converting seurat obj inot other format
SaveH5Seurat(harmonized_seurat, filename = "./scvelo_files/harmonized_seurat.h5Seurat")
Convert("./scvelo_files/harmonized_seurat.h5Seurat", dest = "h5ad")
```
RNA velocity analysis can accept either a single sample or an integrated object as input. When examining the cell transition between sub-clusters of cells, our objective is to observe the transitions between groups of cells present in the dataset. Integrating samples into one integrated object allows us to combine cells from different samples into a single entity. This approach has both advantages and disadvantages. The advantage is that we have a more diverse set of cells to analyze. However, the disadvantage is that each sample is unique due to tumor heterogeneity, and studying them separately provides valuable insights. In this case, we will not be analyzing the integrated object, although we have generated one using the code mentioned above.

The rest of analysis is performed in python environment (notebook ....):

```python
#import libs
import os
import scvelo as scv
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import statsmodels.graphics.mosaicplot as mplt
#import pyreader

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

# setting working directory
os.chdir("C:\\Users\\qaedi\\OneDrive - Queen's University\\Documents\\scRNA_velocity\\scvelo_files")
# setting working directory
os.chdir("C:\\Users\\qaedi\\OneDrive - Queen's University\\Documents\\scRNA_velocity\\scvelo_files")
```

## SRR12603783 HG MIBC
```python
# reading file
adata = scv.read("SRR12603783.h5ad")
adata

#AnnData object with n_obs × n_vars = 6545 × 36601
#    obs: 'orig.ident', 'nCount_spliced', 'nFeature_spliced', 'nCount_unspliced', #'nFeature_unspliced', 'nCount_ambiguous', 'nFeature_ambiguous', 'cells', 'gender', 'age', 'Grade', #'Invasiveness', 'clusters', 'nCount_RNA', 'nFeature_RNA', 'nCount_SCT', 'nFeature_SCT', #'SCT_snn_res.0.1', 'SCT_snn_res.0.2', 'SCT_snn_res.0.4', 'SCT_snn_res.0.6', 'SCT_snn_res.0.8', #'seurat_clusters', 'epi_cluster'
#    var: 'features', 'ambiguous_features', 'spliced_features', 'unspliced_features'
#    obsm: 'X_umap'
#    layers: 'ambiguous', 'spliced', 'unspliced'


# To see cells in each clusters
adata.obs.clusters.value_counts()
#1    2626
#4    2127
#0     795
#2     703
#5     152
#3     142
#Name: clusters, dtype: int64

# Define a dictionary mapping the numbers to labels
label_mapping = {
    0: 'basal_cell',
    1: 'cancer_associated_luminal_cell',
    2: 'differentiated_luminal_cell',
    3: 'unique_luminal_cell',
    4: 'immunomodulatory_luminal_cell',
    5: 'adhesion_signaling_luminal_cell'
}

# Replace the numbers with their corresponding labels
adata.obs['clusters'] =adata.obs['clusters'].replace(label_mapping)
adata.obs['clusters'] = adata.obs['clusters'].astype('category')

# Define a dictionary mapping the numbers to labels
label_mapping = {
    0: 'basal_cell',
    1: 'cancer_associated_luminal_cell',
    2: 'differentiated_luminal_cell',
    3: 'unique_luminal_cell',
    4: 'immunomodulatory_luminal_cell',
    5: 'adhesion_signaling_luminal_cell'
}

# Replace the numbers with their corresponding labels
adata.obs['clusters'] =adata.obs['clusters'].replace(label_mapping)
adata.obs['clusters'] = adata.obs['clusters'].astype('category')
```

### Basic pre-processing
```python
scv.pl.proportions(adata, figsize=(18,2), save= 'basic_features_SRR12603783.png')
```
![basic_features_SRR12603783.png](https://github.com/hamidghaedi/scRNA_velocity_analysis/blob/main/image/scvelo_proportions_basic_features_SRR12603783.png)

After basic preprocessing (gene selection and normalization), we compute the first- and second-order moments (means and uncentered variances) for velocity estimation:

```python
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
#Filtered out 31130 genes that are detected 20 counts (shared).
#Normalized count data: X, spliced, unspliced.
#Extracted 2000 highly variable genes.
#Logarithmized X.
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
#computing neighbors
#    finished (0:00:06) --> added 
#    'distances' and 'connectivities', weighted adjacency matrices (adata.obsp)
#computing moments based on connectivities
#    finished (0:00:01) --> added 
#    'Ms' and 'Mu', moments of un/spliced abundances (adata.layers)
```

### Velocity Tools _ dynamic model
The core of the software is the efficient and robust estimation of velocities, obtained with:

```python
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

#recovering dynamics (using 1/16 cores)
#Error displaying widget: model not found
#    finished (0:03:31) --> added 
#    'fit_pars', fitted parameters for splicing dynamics (adata.var)
#computing velocities
#    finished (0:00:01) --> added 
#    'velocity', velocity vectors for each individual cell (adata.layers)
#computing velocity graph (using 1/16 cores)
#Error displaying widget: model not found
#    finished (0:00:11) --> added 
#    'velocity_graph', sparse matrix with cosine correlations (adata.uns)
```
### Velocity stream plot

```python
scv.pl.velocity_embedding_stream(adata, basis="umap", color="SCT_snn_res.0.4", title= 'SRR12603783 [HG, MIBC]', save= 'srr12603783_velocity_stream.png')
#computing velocity embedding
#    finished (0:00:01) --> added
#    'velocity_umap', embedded velocity vectors (adata.obsm)
#saving figure to file ./figures/scvelo_srr12603783_velocity_stream.png
```
![scvelo_srr12603783_velocity_stream.png](https://github.com/hamidghaedi/scRNA_velocity_analysis/blob/main/image/scvelo_srr12603783_velocity_stream.png)