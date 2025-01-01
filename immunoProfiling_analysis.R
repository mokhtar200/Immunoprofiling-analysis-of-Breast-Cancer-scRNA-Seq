# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# Step 1: Load and preprocess data
# Load the single-cell RNA-seq data
# Ensure your FASTQ files are in the correct directory
sc_data <- Read10X(data.dir = "path_to_your_data/")  # Update with actual path

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, project = "BreastCancer_ImmunoProfiling")

# Step 2: Quality Control and Filtering
# Filtering cells based on thresholds for genes detected and percent of mitochondrial genes
seurat_obj <- seurat_obj %>%
  subset(nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10)

# Step 3: Normalization and scaling
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Step 4: Dimensionality reduction
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Step 5: Clustering of cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Step 6: Immunoprofiling Analysis
# Identifying immune cell populations based on marker genes
immune_markers <- c("CD3D", "CD19", "CD14", "CD8A", "PDCD1")  # Example immune markers
seurat_obj <- AddModuleScore(seurat_obj, features = list(immune_markers), name = "Immune_Cell_Signature")

# Step 7: Visualize immune cell populations
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "tsne", group.by = "ident") + ggtitle("Breast Cancer Immune Profiling")

# Step 8: Differential Expression Analysis
# Find differentially expressed genes between immune cell populations
de_genes <- FindMarkers(seurat_obj, ident.1 = "1", ident.2 = "2")

# Step 9: Save results and Seurat object for future use
write.csv(de_genes, "immune_marker_differential_expression.csv")
saveRDS(seurat_obj, "breast_cancer_immunoprofiling_seurat_obj.rds")
