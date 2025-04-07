install.packages("BiocManager")  # Ensure BiocManager is installed
BiocManager::install("glmGamPoi")

install.packages("magrittr")
install.packages("pROC")
# Load Required Libraries
library(Seurat)
library(Matrix)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)
library(Harmony)
library(readxl)
library(glmGamPoi) 
library(magrittr)


# Set Working Directory (Modify this)
setwd("/Users/pradeepchowdary/Desktop/project")

# Load Gene Count Matrix (FeatureCounts Output)
counts <- read.table("gene_counts.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
counts <- counts[, grep("SRR", colnames(counts))]  # Keep only sample columns


# Ensure counts are numeric
counts_matrix <- as.matrix(counts)
mode(counts_matrix) <- "numeric"

# Convert to sparse matrix
counts_sparse <- Matrix(counts_matrix, sparse = TRUE)


# Load Metadata (Contains Tumor/Normal Labels)
metadata <- read_excel("metadata.xlsx")  # Ensure it has SampleID & Condition (Tumor/Normal)

# Create Seurat Object
seurat_obj <- CreateSeuratObject(counts = counts, project = "Tumor_vs_Normal", min.cells = 3, min.features = 200)
seurat_obj <- AddMetaData(seurat_obj, metadata)

# Rename Row Names in Metadata to Match Seurat Object
rownames(seurat_obj@meta.data) <- gsub("Aligned.sortedByCoord.out.bam", "", rownames(seurat_obj@meta.data))

# Load required libraries
library(Seurat)
library(biomaRt)

# Step 1: Connect to Ensembl and Retrieve Gene Symbol Mapping
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")  

# Get Ensembl-to-Gene Symbol Mapping
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      mart = ensembl)

# Convert to named vector for easy lookup
ensembl_to_symbol <- setNames(gene_mapping$hgnc_symbol, gene_mapping$ensembl_gene_id)

# Step 2: Convert Ensembl IDs in Seurat Object
rownames(seurat_obj) <- ensembl_to_symbol[rownames(seurat_obj)]

# Step 3: Remove Genes Without Valid Symbols
valid_genes <- !is.na(rownames(seurat_obj))
seurat_obj <- subset(seurat_obj, features = rownames(seurat_obj)[valid_genes])

# Step 4: Recreate Seurat Object to Ensure Correct Filtering
counts_filtered <- GetAssayData(seurat_obj, layer = "counts")[valid_genes, ]
meta_data <- seurat_obj@meta.data
seurat_obj <- CreateSeuratObject(counts = counts_filtered, meta.data = meta_data)

# Step 5: Compute Mitochondrial Percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Verify final object dimensions
dim(seurat_obj)  # Should return (filtered genes, number of cells)

# Visualize QC metrics: nFeature_RNA, nCount_RNA, and percent.mt
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Apply filtering: Remove cells with extreme values
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)

# Check dimensions after filtering
dim(seurat_obj)  # Should return (genes, high-quality cells)


# Normalization 
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
dim(seurat_obj[["SCT"]])  # Should return (number of genes, number of cells)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Visualize top variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

dim(seurat_obj)  # Should return (number of genes, number of cells)


# PCA & Dimensionality Reduction
# Dynamically adjust npcs to avoid the error
max_pcs <- min(nrow(seurat_obj) - 1, ncol(seurat_obj) - 1)  # Ensure valid PCA dimensions
seurat_obj <- RunPCA(seurat_obj, npcs = max_pcs)

# Visualize PCA
ElbowPlot(seurat_obj)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.2)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10, perplexity = 5)

# Plot clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(seurat_obj, reduction = "tsne", label = TRUE, pt.size = 0.5)

table(Idents(seurat_obj))

cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View top marker genes per cluster
head(cluster_markers)

FeaturePlot(seurat_obj, features = c("ARHGAP9", "DDX5", "GPBP1"), reduction = "umap")



# top 5 marker genes per cluster
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

VlnPlot(seurat_obj, features = c("ARHGAP9", "DDX5", "GPBP1"), group.by = "seurat_clusters")


VlnPlot(seurat_obj, features = c("ARHGAP9", "DDX5", "GPBP1"), group.by = "Condition")


# DEG

head(seurat_obj@meta.data)
Idents(seurat_obj) <- seurat_obj$Condition
table(Idents(seurat_obj)) 


deg_results <- FindMarkers(seurat_obj, ident.1 = "Tumor", ident.2 = "Normal",
                           min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox")

# Apply FDR correction
deg_results$p_val_adj <- p.adjust(deg_results$p_val, method = "fdr")

head(deg_results)

# Save results
write.csv(deg_results, "DEGs_Tumor_vs_Normal.csv")


# Enhanbcedvolcano
EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "DEGs: Tumor vs Normal",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 5.0)

top_degs <- rownames(deg_results[order(deg_results$p_val_adj), ])[1:20]
print(top_degs)


# pheatmap

valid_genes <- intersect(top_degs, rownames(GetAssayData(seurat_obj, slot = "scale.data")))
length(valid_genes)  # Should be > 0

# Extract only valid genes for heatmap
expr_data <- GetAssayData(seurat_obj, slot = "scale.data")[valid_genes, ]

pheatmap(expr_data, 
         annotation_col = seurat_obj@meta.data["Condition"], 
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = TRUE,
         main = "Top DEGs Heatmap (Valid Genes)")


# to validate 
VlnPlot(seurat_obj, features = c("ACTL6A", "TSPAN8", "ID4"), group.by = "Condition")

# biomarker identification
biomarkers <- deg_results %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>%
  arrange(p_val_adj)

head(biomarkers)


library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs (if needed)
# gene ontology
entrez_ids <- mapIds(org.Hs.eg.db, keys = rownames(biomarkers),
                     column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
ego <- enrichGO(gene = entrez_ids, 
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)
dotplot(ego, showCategory = 15, title = "GO Enrichment")

# kegg pathway
ekegg <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
barplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment")

library(pheatmap)
scaled_data <- GetAssayData(seurat_obj, assay = "SCT", layer = "scale.data")
common_genes <- intersect(rownames(biomarkers), rownames(scaled_data))
expr_data <- scaled_data[common_genes, ]
pheatmap(expr_data, annotation_col = seurat_obj@meta.data["Condition"],
         scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Biomarker Expression Heatmap")


library(pROC)

# Example for one biomarker
gene_expr <- FetchData(seurat_obj, vars = "FOXA1")
roc_obj <- roc(seurat_obj$Condition, gene_expr$FOXA1)
plot(roc_obj, main = "ROC Curve for GENE_OF_INTEREST")
auc(roc_obj)#





