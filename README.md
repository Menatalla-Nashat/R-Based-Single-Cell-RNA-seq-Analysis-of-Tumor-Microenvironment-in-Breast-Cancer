Single-cell transcriptomic analysis reveals CAF heterogeneity and signaling programs in breast cancer

1. Project Description

This project explored the diversity of cancer-associated fibroblasts (CAFs) within the breast cancer tumor microenvironment through scRNA-seq analysis.

The analysis was performed in R using Seurat and CellChat packages, and the entire process contains quality control, normalization, clustering, cell type annotation, CAF subtype characterization, differential expression, and cell-cell communication analysis.

To understand the communication pattern and biological functions of different CAF subtypes, cell-cell communication networks of ligands and receptors, and signal pathways between CAF subtypes were explored using CellChat.


2. How to Run the Tool

2.1. Open RStudio.

2.2. Set the working directory to the project folder.

2.3. Make sure the following input files are available in the project directory:
   - matrix.mtx
   - barcodes.tsv
   - genes.tsv
   - metadata.csv

2.4. Install the required R packages:

install.packages(c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "harmony",
  "dplyr",
  "patchwork",
  "CellChat"))

2.5. Run the preprocessing and CAF analysis script:
Data loading using Read10X
Quality control and filtering
Normalization and scaling
PCA and Harmony batch correction
UMAP visualization and clustering
Cell type annotation
CAF extraction and subclustering
CAF subtype annotation

2.6. Run the CellChat analysis script:
Create CellChat object
Load CellChat human database
Detect overexpressed genes and ligand-receptor interactions
Compute communication probabilities
Generate signaling pathway analysis
Create visualization plots
Output files will be automatically generated and saved in the results folders and working directory.

2.7. Generated outputs include:

UMAP plots
Marker analysis files
CAF subtype plots
Circle plot
Heatmaps
Bubble plot
Pathway ranking plot
Signaling role heatmap

3. Required Packages

  "Seurat","Matrix","ggplot2","harmony","dplyr","patchwork","CellChat"

4. Dataset Information

The project uses GSE176078 human breast cancer single-cell RNA sequencing (scRNA-seq) data.

The data consists of:
- A count matrix of gene expression levels
- A barcode file
- A genes/features file
- Annotations related to metadata

The input files needed to run the analysis are:
- matrix.mtx
- barcodes.tsv
- genes.tsv
- metadata.csv

This data contains ~100,064 cells from 26 primary breast tumors, including:
- ER+
- HER2 Positive
- TNBC

All the data files must be located inside the main directory of the project before you run the scripts.

For CAF-specific analysis, we use the processed Seurat object:

- final_CAF_object.rds

which results from CAF extraction and subclustering.

5. Parameter Values Used

Quality Control Parameters
- Min detected genes per cell: nFeature_RNA > 300
- Max detected genes per cell: nFeature_RNA < 6000
- Max mitochondrial percentage: percent.mt < 10


Normalization Parameters
- Normalization method: LogNormalize
- Scale factor: 10000


Highly Variable Genes (HVGs)
- Selection method: vst
- Num of variable genes: 2000



PCA and Dim Reduction
- PCA dimensions used: 1:30
- UMAP dimensions used: 1:30
- Random seed: 42


Clustering Parameters
- Main clustering resolution: 0.1
- CAF reclustering resolution: 0.2


Differential Expression Parameters
- Minimum percentage expression: min.pct = 0.25
- Log fold-change threshold: logfc.threshold = 0.25


CellChat Parameters
- Database used: CellChatDB.human
- Source CAF subtypes:
 - myCAF
 - iCAF

- Target CAF subtypes:
 - matrixCAF
 - apCAF


Visualization Parameters
- Output format: PNG
- Width: 2000
- Height: 2000
- Resolution: 300 dpi


6. Screenshots of Results

The figures which represents the output of the analyses have been collected in the folder named screenshots.

Included images:
- 01_Cbefore_filtering.png : QC metrics before filtering.
- 03_Cafter_filtering.png : QC metrics after filtering.
- UMAP_CLUSTERS.png : UMAP visualization of identified clusters.
- Volcano Plot_TNBC_vs ER plus_DEGs.png : Volcano plot showing differentially expressed genes between TNBC and ER+ breast cancer cells.
- Major cell type annotation.png : Major cell populations annotated from scRNA-seq analysis.
- Annotation for CAF subtypes.png : Major cell types annotated.
- CAF subtypes Markers Feature Plot.png : Marker gene expression profiles for CAF subtypes.
- circle_plot.png : CAF subtype intercellular signaling communication network visualized by CellChat.
- signaling_rle_eatmap.png: signaling pathways' activity heatmap.
- sRanking_pahways2.png: ranking of signaling pathways in terms of information transfer.


7. How to Reproduce the Work

1. Download all dataset files and place in the root directory of the project.
- matrix.mtx
- barcodes.tsv
- genes.tsv
- metadata.csv
2. Open the project in R Studio.
3. Install all necessary R libraries:
- Seurat
- Matrix
- ggplot2
- harmony
- dplyr
- patchwork
- CellChat
4. Execute sequentially the preprocessing and clustering script:
- Quality control and filtering
- Normalization
- PCA and Harmony batch correction
- UMAP clustering
- Cell type annotation
- Extract and recluster CAF populations
5. Execute the CellChat analysis script:
- Inference of ligand-receptor interactions
- Calculation of communication probabilities
- Signaling pathway analysis
- CellChat visualizations
6. Output files and images will be saved into the respective results folders or the working directory automatically.
7. Some representative output images are found within the screenshots folder.
