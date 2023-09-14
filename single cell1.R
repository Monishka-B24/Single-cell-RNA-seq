library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data = Read10X(data.dir = "C:/Users/Monishka/Downloads/cancer1/GSE143868_RAW/GSM4275611_D322_Biop_Int1.tar/GSM4275611_D322_Biop_Int1/D322_Biop_Int1")
pbmc = CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features = 200)
pbmc.data[1:50, 1:10]

pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2200 & percent.mt < 15)

pbmc = NormalizeData(pbmc)
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 = head(VariableFeatures(pbmc), 10)

plot1 = VariableFeaturePlot(pbmc)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)
pbmc@assays$RNA@scale.data[1:50, 1:5]

pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(pbmc)

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 1)
head(pbmc@meta.data)

pbmc = RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap", label = T)
pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(pbmc.markers)

filtered_markers = pbmc.markers %>% filter(p_val_adj < 0.01, avg_log2FC > 1)

cluster_genes_list <- vector("list", length(unique(filtered_markers$cluster)))

# Populate the list with gene names for each cluster
for (i in seq_along(cluster_genes_list)) {
  cluster <- unique(filtered_markers$cluster)[i]
  genes <- filtered_markers$gene[filtered_markers$cluster == cluster]
  cluster_genes_list[[i]] <- genes
}
# Determine the maximum number of genes in a cluster
max_genes <- max(sapply(cluster_genes_list, length))

# Fill in shorter lists with NA to match max_genes
cluster_genes_list <- lapply(cluster_genes_list, function(x) {
  length_diff <- max_genes - length(x)
  c(x, rep(NA, length_diff))
})

# Create a data frame with gene names for each cluster
cluster_genes_df <- as.data.frame(do.call(cbind, cluster_genes_list))

cluster_id_to_cell_type <- c("0" ="Basal squamous epithelial cells",
                             "1"= "Trophoblast cells",
                             
                             "2"= "Glandular epithelial cells",
                             "3"="Specialized epithelial cells",
                             "4"="Glandular epithelial cells",
                             "5"="Endothelial cells",
                             "6"="miscellaneous",
                             "7"="miscellaneous",
                             "8"="mix",
                             "9"="mix",
                             "10"="mixture",
                             "11"="Blood & immune cells")

# Assuming you have a Seurat object named 'pbmc' and a mapping of clusters to cell types named 'cluster_to_celltype'

# Assign cell types to clusters and add a new column 'cell_type' to metadata
pbmc@meta.data$cell_type <-cluster_id_to_cell_type[as.character(pbmc@meta.data$seurat_clusters)]

# Now the 'cell_type' column is added to the metadata with assigned cell types
print(pbmc@meta.data)

DimPlot(pbmc,group.by = "cell_type")
