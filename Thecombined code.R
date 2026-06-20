# ============================================================
# Purpose:
#   1. Load dataset using Read10X
#   2. Create Seurat object and add metadata
#   3. Quality control and filtering
#   4. Normalization and highly variable genes
#   5. PCA
#   6. Batch correction using Harmony
#   7. UMAP and clustering BEFORE annotation
#   8. FindAllMarkers BEFORE annotation
#   9. Differential expression TNBC vs ER+
#   10. Marker-based annotation
#   11. Reference mapping validation
#   12. CAF extraction and reclustering
# ============================================================

# -----------------------------
# Step 0: Clean session
# -----------------------------
rm(list = ls())
gc()

# -----------------------------
# Step 1: Load packages
# -----------------------------
required_packages <- c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "harmony",
  "dplyr",
  "patchwork"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", ")
    )
  )
}

library(Seurat)
library(Matrix)
library(ggplot2)
library(harmony)
library(dplyr)
library(patchwork)

set.seed(42)
options(stringsAsFactors = FALSE)


# -----------------------------
# Step 2: Define project directories
# -----------------------------
base_dir <- "C:/Users/A/Desktop/Bioinfo project/scRNA_data-project"

data_dir <- base_dir

results_dir <- file.path(
  base_dir,
  "results"
)

processed_dir <- file.path(
  base_dir,
  "processed"
)

qc_dir <- file.path(results_dir, "QC")
pca_dir <- file.path(results_dir, "PCA")
marker_dir <- file.path(results_dir, "Markers")
annotation_dir <- file.path(results_dir, "Annotation")
caf_dir <- file.path(results_dir, "CAF_analysis")
umap_dir <- file.path(results_dir, "UMAP")
de_dir <- file.path(results_dir, "Differential_expression")

for (folder in c(
  results_dir,
  processed_dir,
  qc_dir,
  pca_dir,
  marker_dir,
  annotation_dir,
  caf_dir,
  umap_dir,
  de_dir
)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
# -----------------------------
# Step 3: Prepare Read10X files
# -----------------------------
tenx_dir <- file.path(data_dir, "tenx")

dir.create(
  tenx_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

raw_genes <- read.delim(
  file.path(data_dir, "genes.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(raw_genes) == 1) {
  
  features_fixed <- data.frame(
    V1 = raw_genes[[1]],
    V2 = raw_genes[[1]],
    V3 = "Gene Expression"
  )
  
} else if (ncol(raw_genes) == 2) {
  
  features_fixed <- data.frame(
    V1 = raw_genes[[1]],
    V2 = raw_genes[[2]],
    V3 = "Gene Expression"
  )
  
} else {
  
  features_fixed <- raw_genes[, 1:3]
}

write.table(
  features_fixed,
  file = file.path(tenx_dir, "features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

file.copy(
  file.path(data_dir, "barcodes.tsv"),
  file.path(tenx_dir, "barcodes.tsv"),
  overwrite = TRUE
)

file.copy(
  file.path(data_dir, "matrix.mtx"),
  file.path(tenx_dir, "matrix.mtx"),
  overwrite = TRUE
)

gzip_copy <- function(input_file, output_file) {
  
  input_con <- file(input_file, "rb")
  output_con <- gzfile(output_file, "wb")
  
  repeat {
    
    bytes <- readBin(
      input_con,
      what = "raw",
      n = 1e6
    )
    
    if (length(bytes) == 0) break
    
    writeBin(bytes, output_con)
  }
  
  close(input_con)
  close(output_con)
}

gzip_copy(
  file.path(tenx_dir, "features.tsv"),
  file.path(tenx_dir, "features.tsv.gz")
)

gzip_copy(
  file.path(tenx_dir, "barcodes.tsv"),
  file.path(tenx_dir, "barcodes.tsv.gz")
)

gzip_copy(
  file.path(tenx_dir, "matrix.mtx"),
  file.path(tenx_dir, "matrix.mtx.gz")
)

# -----------------------------
# Step 4: Read dataset
# -----------------------------
counts_raw <- Read10X(
  data.dir = tenx_dir,
  gene.column = 2,
  cell.column = 1
)

if (is.list(counts_raw)) {
  
  counts <- counts_raw[[1]]
  
} else {
  
  counts <- counts_raw
}

print(dim(counts))

# -----------------------------
# Step 5: Create Seurat object
# -----------------------------
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "BRCA_scRNA",
  min.cells = 3,
  min.features = 200
)

print(seurat_obj)

# -----------------------------
# Step 6: Add metadata
# -----------------------------
metadata <- read.csv(
  file.path(data_dir, "metadata.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

metadata <- metadata[
  colnames(seurat_obj),
  ,
  drop = FALSE
]

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata
)

# -----------------------------
# Step 7: Mito percentage
# -----------------------------
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)


# -----------------------------
# Step 8: QC before filtering
# -----------------------------
qc_vln_before <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  group.by = ifelse("subtype" %in% colnames(seurat_obj@meta.data), "subtype", NULL),
  ncol = 3,
  pt.size = 0
)

print(qc_vln_before)
ggsave(file.path(qc_dir, "01_QC_before_filtering.png"), qc_vln_before, width = 12, height = 5, dpi = 300)

qc_summary_before <- summary(seurat_obj@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mito")])
writeLines(capture.output(qc_summary_before), file.path(qc_dir, "01_QC_summary_before_filtering.txt"))

# -----------------------------
# Step 8.1: QC filtering
# -----------------------------
cells_before_filtering <- ncol(seurat_obj)

seurat_obj_filtered <- subset(
  seurat_obj,
  subset = nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mito < 10
)

cells_after_filtering <- ncol(seurat_obj_filtered)

filtering_summary <- data.frame(
  Step = c("Before filtering", "After filtering", "Removed cells"),
  Number_of_cells = c(
    cells_before_filtering,
    cells_after_filtering,
    cells_before_filtering - cells_after_filtering
  )
)

print(filtering_summary)
write.csv(filtering_summary, file.path(qc_dir, "02_filtering_summary.csv"), row.names = FALSE)

# -----------------------------
# Step 8.2: QC after filtering
# -----------------------------
qc_vln_after <- VlnPlot(
  seurat_obj_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
  group.by = ifelse("subtype" %in% colnames(seurat_obj_filtered@meta.data), "subtype", NULL),
  ncol = 3,
  pt.size = 0
)

print(qc_vln_after)
ggsave(file.path(qc_dir, "03_QC_after_filtering.png"), qc_vln_after, width = 12, height = 5, dpi = 300)

qc_summary_after <- summary(seurat_obj_filtered@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mito")])
writeLines(capture.output(qc_summary_after), file.path(qc_dir, "03_QC_summary_after_filtering.txt"))

saveRDS(seurat_obj_filtered, file = file.path(processed_dir, "02_filtered_Read10X.rds"))

# -----------------------------
# Step 9: Normalize
# -----------------------------
seurat_obj_filtered <- NormalizeData(
  seurat_obj_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# -----------------------------
# Step 10: HVGs
# -----------------------------
seurat_obj_filtered <- FindVariableFeatures(
  seurat_obj_filtered,
  selection.method = "vst",
  nfeatures = 2000
)

VariableFeaturePlot(seurat_obj_filtered)

# -----------------------------
# Step 11: Scale
# -----------------------------
variable_genes <- VariableFeatures(
  seurat_obj_filtered
)

seurat_obj_filtered <- ScaleData(
  seurat_obj_filtered,
  features = variable_genes
)

# -----------------------------
# Step 12: PCA
# -----------------------------
seurat_obj_filtered <- RunPCA(
  seurat_obj_filtered,
  features = variable_genes
)

ElbowPlot(
  seurat_obj_filtered,
  ndims = 50
)

# -----------------------------
# Step 13: Batch correction
# -----------------------------
seurat_obj_filtered <- harmony::RunHarmony(
  object = seurat_obj_filtered,
  group.by.vars = "orig.ident"
)

# -----------------------------
# Step 14: UMAP using Harmony
# -----------------------------
seurat_obj_filtered <- RunUMAP(
  seurat_obj_filtered,
  reduction = "harmony",
  dims = 1:30,
  seed.use = 42
)

# -----------------------------
# Step 15: Find neighbors
# -----------------------------
seurat_obj_filtered <- FindNeighbors(
  seurat_obj_filtered,
  reduction = "harmony",
  dims = 1:30
)

# -----------------------------
# Step 16: Clustering
# -----------------------------
seurat_obj_filtered <- FindClusters(
  seurat_obj_filtered,
  resolution = 0.1
)

# -----------------------------
# Step 17: UMAP clusters
# -----------------------------
DimPlot(
  seurat_obj_filtered,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
)

ggsave(
  file.path(umap_dir, "UMAP_clusters_resolution_0.1.png"),
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# Step 18: FindAllMarkers
# -----------------------------
cluster_markers <- FindAllMarkers(
  object = seurat_obj_filtered,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  cluster_markers,
  file.path(
    marker_dir,
    "all_cluster_markers.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# Step 19: Top markers
# -----------------------------
top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    n = 10,
    order_by = avg_log2FC
  )

write.csv(
  top10_markers,
  file.path(
    marker_dir,
    "top10_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

# -----------------------------
# Step 20: Visualize markers
# -----------------------------
FeaturePlot(
  seurat_obj_filtered,
  features = c(
    "COL1A1",
    "EPCAM",
    "CD3D",
    "MS4A1"
  )
)

DotPlot(
  seurat_obj_filtered,
  features = c(
    "COL1A1",
    "DCN",
    "EPCAM",
    "KRT18",
    "CD3D",
    "IL7R",
    "MS4A1",
    "CD79A"
  )
) + RotatedAxis()


# -----------------------------
# Step 21: UMAP by subtype
# -----------------------------
possible_subtype_columns <- c(
  "subtype",
  "Subtype",
  "cancer_subtype",
  "Cancer_Subtype",
  "condition",
  "Condition",
  "diagnosis",
  "Diagnosis",
  "group",
  "Group"
)

subtype_column <- NULL

for (col in possible_subtype_columns) {
  
  if (col %in% colnames(seurat_obj_filtered@meta.data)) {
    
    subtype_values <- seurat_obj_filtered@meta.data[[col]]
    subtype_values <- subtype_values[!is.na(subtype_values)]
    
    if (length(unique(subtype_values)) > 1) {
      subtype_column <- col
      break
    }
  }
}

if (!is.null(subtype_column)) {
  
  subtype_umap <- DimPlot(
    seurat_obj_filtered,
    reduction = "umap",
    group.by = subtype_column
  )
  
  ggsave(
    file.path(
      umap_dir,
      paste0("UMAP_by_", make.names(subtype_column), ".png")
    ),
    subtype_umap,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# -----------------------------
# Step 22: Differential expression TNBC vs ER+
# -----------------------------
if (!is.null(subtype_column)) {
  
  subtype_values <- unique(
    as.character(seurat_obj_filtered@meta.data[[subtype_column]])
  )
  
  subtype_values <- subtype_values[!is.na(subtype_values)]
  
  group_1 <- NA
  group_2 <- NA
  
  if ("TNBC" %in% subtype_values) {
    group_1 <- "TNBC"
  } else if ("Triple-negative" %in% subtype_values) {
    group_1 <- "Triple-negative"
  } else if ("Triple Negative" %in% subtype_values) {
    group_1 <- "Triple Negative"
  }
  
  if ("ER+" %in% subtype_values) {
    group_2 <- "ER+"
  } else if ("ER-positive" %in% subtype_values) {
    group_2 <- "ER-positive"
  } else if ("ER positive" %in% subtype_values) {
    group_2 <- "ER positive"
  }
  
  if (!is.na(group_1) && !is.na(group_2)) {
    
    Idents(seurat_obj_filtered) <- seurat_obj_filtered@meta.data[[subtype_column]]
    
    de_results <- FindMarkers(
      object = seurat_obj_filtered,
      ident.1 = group_1,
      ident.2 = group_2,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "wilcox"
    )
    
    de_results$gene <- rownames(de_results)
    
    de_fc_col <- if ("avg_log2FC" %in% colnames(de_results)) {
      "avg_log2FC"
    } else {
      "avg_logFC"
    }
    
    de_results$significance <- "Not significant"
    
    de_results$significance[
      de_results[[de_fc_col]] > 0.25 &
        de_results$p_val_adj < 0.05
    ] <- paste0("Higher in ", group_1)
    
    de_results$significance[
      de_results[[de_fc_col]] < -0.25 &
        de_results$p_val_adj < 0.05
    ] <- paste0("Higher in ", group_2)
    
    write.csv(
      de_results,
      file.path(
        de_dir,
        paste0(
          "DE_",
          make.names(group_1),
          "_vs_",
          make.names(group_2),
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    de_results$p_val_adj_plot <- ifelse(
      de_results$p_val_adj == 0,
      .Machine$double.xmin,
      de_results$p_val_adj
    )
    
    volcano_plot <- ggplot(
      de_results,
      aes(
        x = .data[[de_fc_col]],
        y = -log10(p_val_adj_plot),
        color = significance
      )
    ) +
      geom_point(alpha = 0.7, size = 1.5) +
      geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_minimal() +
      labs(
        title = paste0(group_1, " vs ", group_2),
        x = "Average log2 fold change",
        y = "-log10 adjusted p-value",
        color = "Category"
      )
    
    ggsave(
      file.path(
        de_dir,
        paste0(
          "Volcano_",
          make.names(group_1),
          "_vs_",
          make.names(group_2),
          ".png"
        )
      ),
      volcano_plot,
      width = 9,
      height = 7,
      dpi = 300
    )
    
    top_de_genes <- de_results %>%
      filter(p_val_adj < 0.05) %>%
      arrange(desc(abs(.data[[de_fc_col]]))) %>%
      slice_head(n = 20) %>%
      pull(gene)
    
    top_de_genes <- top_de_genes[
      top_de_genes %in% rownames(seurat_obj_filtered)
    ]
    
    if (length(top_de_genes) > 2) {
      
      seurat_obj_filtered <- ScaleData(
        seurat_obj_filtered,
        features = unique(c(variable_genes, top_de_genes)),
        verbose = FALSE
      )
      
      de_heatmap <- DoHeatmap(
        seurat_obj_filtered,
        features = top_de_genes,
        group.by = subtype_column,
        size = 3
      )
      
      ggsave(
        file.path(
          de_dir,
          paste0(
            "Heatmap_top_DE_",
            make.names(group_1),
            "_vs_",
            make.names(group_2),
            ".png"
          )
        ),
        de_heatmap,
        width = 12,
        height = 9,
        dpi = 300
      )
    }
  }
}

# ============================================================
# Step 23: MANUAL ANNOTATION
# ============================================================

# Use marker genes + CellMarker database
# https://bio-bigdata.hrbmu.edu.cn/CellMarker/

# Example annotation
new_cluster_ids_1 <- c(
  "T_cells",
  "Epithelial",
  "Myeloid",
  "Endothelial",
  "CAFs",
  "Pericytes",
  "Plasma",
  "B_cells",
  "Proliferating",
  "Basal_epithelial",
  "pDCs",
  "Unknown"
)

names(new_cluster_ids_1) <- levels(seurat_obj_filtered)

seurat_obj_filtered <- RenameIdents(
  seurat_obj_filtered,
  new_cluster_ids_1
)
names(new_cluster_ids_1) <- levels(
  seurat_obj_filtered
)

seurat_obj_filtered <- RenameIdents(
  seurat_obj_filtered,
  new_cluster_ids_1
)

# Save annotation
seurat_obj_filtered$annotated_celltype <- Idents(
  seurat_obj_filtered
)

# -----------------------------
# Step 24: Annotated UMAP
# -----------------------------
DimPlot(
  seurat_obj_filtered,
  reduction = "umap",
  group.by = "annotated_celltype",
  label = TRUE,
  repel = TRUE
)

# ============================================================
# Step 25: Reference mapping validation
# ============================================================

# Compare your annotation with metadata annotation

table(
  Your_Annotation =
    seurat_obj_filtered$annotated_celltype,
  
  Metadata_Annotation =
    seurat_obj_filtered$celltype_major
)

# ============================================================
# Step 26: CAF extraction
# ============================================================

caf_clusters <- c("CAFs")

caf_obj <- subset(
  seurat_obj_filtered,
  idents = caf_clusters
)

print(caf_obj)

# ============================================================
# Step 27: CAF reclustering
# ============================================================

caf_obj <- NormalizeData(caf_obj)

caf_obj <- FindVariableFeatures(
  caf_obj,
  selection.method = "vst",
  nfeatures = 2000
)

caf_variable_genes <- VariableFeatures(caf_obj)

caf_obj <- ScaleData(
  caf_obj,
  features = caf_variable_genes
)

caf_obj <- RunPCA(
  caf_obj,
  features = caf_variable_genes
)

caf_obj <- RunHarmony(
  object = caf_obj,
  group.by.vars = "orig.ident",
)

caf_obj <- RunUMAP(
  caf_obj,
  reduction = "harmony",
  dims = 1:20
)

caf_obj <- FindNeighbors(
  caf_obj,
  reduction = "harmony",
  dims = 1:20
)

caf_obj <- FindClusters(
  caf_obj,
  resolution = 0.2
)

DimPlot(
  caf_obj,
  reduction = "umap",
  label = TRUE
)

# ============================================================
# Step 28: CAF markers
# ============================================================

caf_markers <- FindAllMarkers(
  caf_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  caf_markers,
  file.path(
    caf_dir,
    "CAF_markers.csv"
  ),
  row.names = FALSE
)

# ============================================================
# Step 29: CAF subtype annotation
# ============================================================

# myCAF:
# ACTA2 TAGLN MYL9

# iCAF:
# CXCL12 CXCL14 IL6

# matrixCAF:
# COL1A1 COL3A1 FN1

# apCAF:
# HLA-DRA CD74

# proliferativeCAF:
# MKI67 TOP2A

FeaturePlot(
  caf_obj,
  features = c(
    "ACTA2",
    "CXCL12",
    "COL1A1",
    "HLA-DRA",
    "MKI67"
  )
)

DotPlot(
  caf_obj,
  features = c(
    "ACTA2",
    "TAGLN",
    "CXCL12",
    "IL6",
    "COL1A1",
    "FN1",
    "HLA-DRA",
    "CD74",
    "MKI67",
    "TOP2A"
  )
) + RotatedAxis()

# ============================================================
# Step 30: Manual CAF subtype annotation
# ============================================================

caf_new_ids <- c(
  "myCAF",
  "iCAF",
  "matrixCAF",
  "apCAF",
  "proliferativeCAF"
)

names(caf_new_ids) <- levels(caf_obj)

caf_obj <- RenameIdents(
  caf_obj,
  caf_new_ids
)

caf_obj$CAF_subtype <- Idents(caf_obj)

DimPlot(
  caf_obj,
  reduction = "umap",
  group.by = "CAF_subtype",
  label = TRUE
)

# ============================================================
# Step 31: Save objects
# ============================================================

saveRDS(
  seurat_obj_filtered,
  file.path(
    processed_dir,
    "final_annotated_seurat.rds"
  )
)

saveRDS(
  caf_obj,
  file.path(
    processed_dir,
    "final_CAF_object.rds"
  )
)
# ============================================================
# Step 32: Cell chat
# ============================================================


library(CellChat)

library(dplyr)

library(Seurat)

caf_obj <- readRDS(
  "C:/Users/A/Desktop/Bioinfo project/scRNA_data-project/processed/final_CAF_object.rds"
)
cellchat_most_accu <- createCellChat(
  object = GetAssayData(caf_obj, layer = "data")
  meta = caf_obj@meta.data,
  group.by = "CAF_subtype")

CellChat_DB <- CellChatDB.human
cellchat_most_accu@DB <- CellChat_DB

# Give me the genes that characterize each cluster.
# The first line in this step will be done depending on Expression level, Comparison between clusters, and Statistical test
# In the second code I am trying to link:ligand (from one cluster), and receptor (in another cluster)


# ------------------> Cellchat process:
#Just give me the genes related to communication. Befor:(All genes (thousands))
#After:(only the important genes (ligand + receptor))

cellchat_most_accu <- subsetData(cellchat_most_accu)

# Give me the genes that characterize each cluster.
# The first line in this step will be done depending on Expression level, Comparison between clusters, and Statistical test
# In the second code I am trying to link:ligand (from one cluster), and receptor (in another cluster)
cellchat_most_accu <- identifyOverExpressedGenes(cellchat_most_accu)
cellchat_most_accu <- identifyOverExpressedInteractions(cellchat_most_accu)

# There are 3 steps :
# 1- calculating the probability that one cell communicates with another cell:
#Ligand in cluster A, Receptor in cluster B
# and Level of expression
# 2- Gathering these signals at the pathways level
# 3- here the code builds the final network, It gathers all the signals and works

cellchat_most_accu <- computeCommunProb(cellchat_most_accu) 
cellchat_most_accu <- computeCommunProbPathway(cellchat_most_accu) 
cellchat_most_accu <- aggregateNet(cellchat_most_accu)



# ------------------> visualization step 

# circle plot
png("circle_plot_2.png", width = 2000, height = 2000, res = 300)
netVisual_circle(cellchat_most_accu@net$count, weight.scale = TRUE, label.edge = FALSE)
dev.off()

# heatmap
png("Heatmap_2.png", width = 2000, height = 2000, res = 300)


netVisual_heatmap(cellchat_most_accu, measure = "count")
dev.off()

# bubble plot (IMPORTANT FIX)
png("Bubble_Blot_2.png", width = 2000, height = 2000, res = 300)
netVisual_bubble(
  cellchat_most_accu,
  sources.use = c("myCAF","iCAF"),
  targets.use = c("matrixCAF","apCAF"),
  remove.isolate = TRUE
)
dev.off()


# PATHWAYS:

cellchat_most_accu@netP$pathways

png("sRanking Pathways_2.png", width = 2000, height = 2000, res = 300)

Rank<- rankNet(cellchat_most_accu, mode = "single")

dev.off()

# Heatmap showing the pathwaysrepresented for each type of CAFs



cellchat_most_accu <- netAnalysis_computeCentrality(cellchat_most_accu)

png("signaling_role_heatmap_2.png", width = 2000, height = 2000, res = 300)
netAnalysis_signalingRole_heatmap(cellchat_most_accu)
dev.off()

