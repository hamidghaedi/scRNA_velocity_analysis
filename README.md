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

## 2) pre-processing in R 

```r
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

```
## 3) velocity analysis in python

```python
import scvelo as scv
adata = scv.read("mouseBM.h5ad")
# to see AnnData structure
adata

## Pre-Processing
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters")
```
![]()

```python
scv.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120)
```
![]()

```python
scv.tl.recover_dynamics(adata)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color="latent_time", color_map="gnuplot")
```
![]()

```python
top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color="seurat_clusters", n_convolve=100)
```
![]()
