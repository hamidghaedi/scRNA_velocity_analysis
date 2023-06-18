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

```python
# To see the proportion of different cell types in sample specefic clsuters:
# Plot the stacked bar chart
plt.figure(figsize=(10, 6))
data.plot(kind='bar', stacked=True)
plt.xlabel('sample-specefic clsuters')
plt.ylabel('Counts')
plt.title('')
plt.legend(loc='upper right')
# Save the plot as an image file
plt.savefig('./figures/srr12603783_frequency_sample_clsuter_in_cell_clusters.png')

plt.show()
```
![srr12603783_frequency_sample_clsuter_in_cell_clusters.png](https://github.com/hamidghaedi/scRNA_velocity_analysis/blob/main/image/srr12603783_frequency_sample_clsuter_in_cell_clusters.png)


The most fine-grained resolution of the velocity vector field we get at single-cell level, with each arrow showing the direction and speed of movement of an individual cell.
```python
scv.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120, title= 'SRR12603783 [HG, MIBC]', save= 'srr12603783_velocity_embedding.png')
```
![srr12603783_velocity_embedding.png](https://github.com/hamidghaedi/scRNA_velocity_analysis/blob/main/image/scvelo_srr12603783_velocity_embedding.png)

### Latent time and top genes

The dynamical model recovers the latent time of the underlying cellular processes. This latent time represents the cell’s internal clock and approximates the real time experienced by cells as they differentiate, based only on its transcriptional dynamics.

```python
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save= "srr12603783_latent_time.png")
```
![srr12603783_latent_time.png](https://github.com/hamidghaedi/scRNA_velocity_analysis/blob/main/image/scvelo_srr12603783_latent_time.png)

### Kinetic rate paramters
The rates of RNA transcription, splicing and degradation are estimated without the need of any experimental data.

They can be useful to better understand the cell identity and phenotypic heterogeneity.
```python
df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata, 'fit*', dropna=True).head()
```
The above code generates a plot , which I askip that for visualization here.But lt's take a look at the result table:

|        |  fit_r2 | fit_alpha | fit_beta | fit_gamma |   fit_t_ | fit_scaling | fit_std_u | fit_std_s | fit_likelihood |   fit_u0 | fit_s0 | fit_pval_steady | fit_steady_u | fit_steady_s | fit_variance | fit_alignment_scaling |
|--------:|----------:|---------:|----------:|---------:|------------:|----------:|----------:|---------------:|---------:|-------:|----------------:|-------------:|-------------:|-------------:|----------------------:|----------|
|    ENO1 |  0.170526 | 3.296458 | 13.728501 | 0.156885 |    7.864689 |  0.026324 |  0.081372 |       3.353864 | 0.219171 |    0.0 |             0.0 |     0.496356 |     0.237947 |     9.392603 |              1.058804 | 2.702026 |
| CAMK2N1 |  0.094471 | 1.098705 | 16.866885 | 0.165689 |    7.069970 |  0.012907 |  0.018368 |       0.821829 | 0.181957 |    0.0 |             0.0 |     0.476764 |     0.057117 |     3.062095 |              1.389385 | 2.526427 |
|    GALE |  0.039663 | 0.127367 |  2.923022 | 0.199394 |    7.737945 |  0.268056 |  0.020084 |       0.129809 | 0.138851 |    0.0 |             0.0 |     0.400852 |     0.067361 |     0.422462 |              2.089054 | 2.294053 |
|    EDN2 |  0.110540 | 0.805942 | 16.574118 | 0.347824 |    3.875601 |  0.035237 |  0.010551 |       0.403489 | 0.000004 |    0.0 |             0.0 |     0.468513 |     0.031373 |     1.050915 |              1.050940 | 1.977689 |
|  SLC2A1 |  0.548495 | 3.259899 | 15.706002 | 0.256389 |   12.335837 |  0.019384 |  0.046016 |       3.556889 | 0.339018 |    0.0 |             0.0 |     0.497239 |     0.153646 |     9.316565 |              0.391465 | 4.524540 |

In the table:

*fit_alpha (Transcription Rate)*: A high value indicates a high transcription rate, suggesting that the gene is actively transcribed, leading to higher mRNA production. Conversely, a low value suggests a lower transcription rate and potentially lower mRNA levels.

*fit_beta (Splicing Rate)*: A high value indicates a high splicing rate, implying that the gene undergoes frequent splicing, resulting in efficient production of mature mRNA. A low value suggests slower or less frequent splicing.

*fit_gamma (Degradation Rate)*: A high value suggests a rapid degradation rate for the gene's mRNA. This could result in lower mRNA abundance or shorter mRNA half-life. Conversely, a low value indicates a slower degradation rate and potentially higher mRNA stability.

*fit_t_ (Switching Time Point)*: A high value indicates a late switching time point, suggesting that the gene is activated or deactivated at a later stage of cellular differentiation or other biological processes. A low value suggests an early switching time point.

*fit_scaling (Scaling Parameter)*: A high value indicates a larger scaling factor used to adjust for under-represented unspliced reads. This suggests that the gene has a relatively higher number of unspliced transcripts compared to the rest of the dataset. A low value suggests a smaller scaling factor and potentially fewer unspliced transcripts.

*fit_std_u and fit_std_s (Standard Deviations)*: These columns represent the standard deviations of unspliced and spliced reads, respectively. A high value indicates greater variability or dispersion in the abundance of unspliced or spliced transcripts for the gene. Conversely, a low value suggests less variability.

*fit_likelihood*: A high value indicates a better fit of the dynamical model to the data for the gene, suggesting that the model accurately captures the transcriptional dynamics. A low value suggests a poorer fit and potential discrepancies between the model and the actual data.

*fit_steady_u and fit_steady_s (Steady-State Levels)*: A high value in either column indicates a higher abundance of the respective transcript type (unspliced or spliced) at the steady state. This suggests a potentially important role for the gene in maintaining stable expression levels. A low value indicates lower steady-state levels.

*fit_pval_steady_u and fit_pval_steady_s (P-values for Steady-State Levels)*: A high value in either column indicates a significant difference in the steady-state levels compared to a reference or null hypothesis. This suggests that the gene's expression at the steady state is distinct and potentially regulated differently. A low value suggests no significant difference.

*fit_variance*: A high value indicates greater variability or uncertainty in the estimated parameters and transcriptional dynamics for the gene. This could be due to complex regulation or heterogeneity. A low value suggests lower variability and more consistent behavior.

*fit_alignment_scaling*: A high value indicates a larger scaling factor used to align the gene-wise latent times to a universal, gene-shared latent time. This suggests a more substantial deviation from the shared latent time in the dynamical model. A low value indicates a smaller deviation and closer alignment to the shared latent time.

In our analysis of bladder cancer, we focused on five specific genes: EDN2, SLC2A1, GBP1, S100A2, and MUC1. These genes play important roles in the transcriptional dynamics and splicing kinetics of bladder cancer cells.

*EDN2*, while showing moderate values in most parameters, may still have relevance in bladder cancer development. Its transcription rate (fit_alpha) and splicing rate (fit_beta) indicate that it contributes to the regulation of gene expression in the context of bladder cancer. Although its impact may not be as strong as other genes, it could still play a role in the disease.

*SLC2A1*, on the other hand, exhibits higher values in multiple parameters. Its elevated transcription rate and splicing rate suggest that SLC2A1 is actively transcribed and spliced in bladder cancer cells. This gene's expression levels may vary significantly, indicating its involvement in the disease process.

*GBP1*, with low to moderate values across most parameters, may have a less prominent role in bladder cancer development compared to other genes. Its relatively low transcription rate and splicing rate suggest a lesser impact on the transcriptional dynamics. However, further studies are needed to determine its specific role in bladder cancer.

*S100A2* emerges as a gene with high values in multiple parameters, indicating its potential significance in bladder cancer. The elevated transcription rate, splicing rate, and standard deviation of spliced reads suggest that S100A2 plays a crucial role in the transcriptional dynamics and splicing kinetics of bladder cancer cells. Its high likelihood and steady-state levels further support its involvement in the disease development.

*MUC1*, similar to EDN2, exhibits low to moderate values in most parameters. Its relatively low transcription rate and splicing rate suggest a lesser impact on the transcriptional dynamics in bladder cancer. While MUC1 may not be as pronounced as other genes in the list, it could still have some influence on the disease.

Overall, our analysis suggests that S100A2 is a key gene in bladder cancer development due to its high values in various parameters related to transcription and splicing. The other genes, including EDN2, SLC2A1, GBP1, and MUC1, may also contribute to the disease but to a lesser extent. Further investigation and experimental validation are necessary to fully understand the roles of these genes in bladder cancer and their potential as therapeutic targets or biomarkers.

Let's look at the genes for each cluster of cells:
```Python
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color="SCT_snn_res.0.4", n_convolve=100 ,colorbar=True, save= "srr12603783_top_genes_kinetic_parameters.png")
```
![srr12603783_top_genes_kinetic_parameters.png](https://github.com/hamidghaedi/scRNA_velocity_analysis/blob/main/image/scvelo_heatmap_srr12603783_top_genes_kinetic_parameters.png)

The same analyses for other HG MIBC samples and also NMIBC samples can be find in the notebook.
