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

epi_seurat$clusters <- Idents(epi_seurat)
 markers <- c("RPSA", "KRT17", "CEACAM6")
FeaturePlot(object = epi_seurat,
            features = markers,
            order = TRUE,
            min.cutoff = "q10",
            reduction = "umap",
            label = TRUE,
            repel = TRUE)
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


# DE and Enrichement analysis

#Differential expression analysis
epi_seurat <- RunDEtest(srt = epi_seurat, group_by = "clusters", fc.threshold = 1.5, only.pos = FALSE)

#
png(filename = "vplcano_DE.png", width = 16, height = 8.135, units = "in", res = 300)
VolcanoPlot(srt = epi_seurat, group_by = "clusters")
dev.off()

#
DEGs <- epi_seurat@tools$DEtest_clusters$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1.5 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
epi_seurat <- AnnotateFeatures(epi_seurat, species = "Homo_sapiens", db = c("TF", "SP"))

png(filename = "heatmap_DE.png", width = 16, height = 8.135, units = "in", res = 300)
ht <- FeatureHeatmap(
  srt = epi_seurat, group.by = "clusters", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Homo_sapiens", db = c("GO_BP", "KEGG", "WikiPathway"), anno_terms = TRUE)
print(ht$plot)
dev.off()

# Enrichement analysis

epi_seurat <- RunEnrichment(
  srt = epi_seurat, group_by = "clusters", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05"
)

png(filename = "enrichment_1.png", width = 16, height = 8.135, units = "in", res = 300)
EnrichmentPlot(
  srt = epi_seurat, group_by = "clusters", group_use = c("basal_cell", "differentiated_luminal_cell"),
  plot_type = "bar"
)
dev.off()

png(filename = "enrichment_2.png", width = 16, height = 8.135, units = "in", res = 300)
EnrichmentPlot(srt = epi_seurat, group_by = "clusters", plot_type = "comparison")
dev.off()

# Enrichment analysis(GSEA)
epi_seurat <- RunGSEA(
  srt = epi_seurat, group_by = "clusters", db = "GO_BP", species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)

png(filename = "enrichment_3.png", width = 16, height = 8.135, units = "in", res = 300)
GSEAPlot(srt = epi_seurat, group_by = "clusters", plot_type = "comparison")
dev.off()


#Trajectory inference
png(filename = "trajectory.png", width = 16, height = 8.135, units = "in", res = 300)
epi_seurat <- RunSlingshot(srt = epi_seurat, group.by = "clusters", reduction = "UMAP")
dev.off()





# reading each sample data

obj <- LoadH5Seurat("./scvelo_files/SRR12603783.h5Seurat")

Idents(obj) <- obj$SCT_snn_res.0.4

# tarjectory
obj <- RunSlingshot(srt = obj, group.by = "SCT_snn_res.0.4", reduction = "UMAP")
