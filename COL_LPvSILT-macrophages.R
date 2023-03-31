library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(xlsx)
library(monocle3)

#Load Datasets
WT_CP.data<-Read10X(data.dir="~/CellRanger/Mortha_Pailin__WT_CP_3pr_v3_1/outs/filtered_feature_bc_matrix")
KO_CP.data<-Read10X(data.dir="~/CellRanger/Mortha_Pailin__KO_CP_3pr_v3_1/outs/filtered_feature_bc_matrix")
WT_LP.data<-Read10X(data.dir="~/CellRanger/Mortha_Pailin__WT_LP_3pr_v3_1/outs/filtered_feature_bc_matrix")
KO_LP.data<-Read10X(data.dir="~/CellRanger/Mortha_Pailin__KO_LP_3pr_v3_1/outs/filtered_feature_bc_matrix")

# Initialize the Seurat object with the raw (non-normalized data).
WT_CP <- CreateSeuratObject(counts = WT_CP.data, project = "WT_CP", min.cells = 3, min.features = 200)
KO_CP <- CreateSeuratObject(counts = KO_CP.data, project = "KO_CP", min.cells = 3, min.features = 200)
WT_LP <- CreateSeuratObject(counts = WT_LP.data, project = "WT_LP", min.cells = 3, min.features = 200)
KO_LP <- CreateSeuratObject(counts = KO_LP.data, project = "KO_LP", min.cells = 3, min.features = 200)

# Calculate percentage of mitochondrial genes and dissociation associated genes
WT_CP[["percent.mt"]] <- PercentageFeatureSet(WT_CP, pattern = "^mt-")
KO_CP[["percent.mt"]] <- PercentageFeatureSet(KO_CP, pattern = "^mt-")
WT_LP[["percent.mt"]] <- PercentageFeatureSet(WT_LP, pattern = "^mt-")
KO_LP[["percent.mt"]] <- PercentageFeatureSet(KO_LP, pattern = "^mt-")

dag <- read.csv(file = "~/dag_genes.csv")
dag_WT_CP <- intersect(rownames(WT_CP), dag$gene2)
dag_KO_CP <- intersect(rownames(KO_CP), dag$gene2)
dag_WT_LP <- intersect(rownames(WT_LP), dag$gene2)
dag_KO_LP <- intersect(rownames(KO_LP), dag$gene2)
WT_CP[["percent.dag"]] <- PercentageFeatureSet(WT_CP, features = dag_WT_CP)
KO_CP[["percent.dag"]] <- PercentageFeatureSet(KO_CP, features = dag_KO_CP)
WT_LP[["percent.dag"]] <- PercentageFeatureSet(WT_LP, features = dag_WT_LP)
KO_LP[["percent.dag"]] <- PercentageFeatureSet(KO_LP, features = dag_KO_LP)

#Filter low quality, dead or stressed cells and putative doublets
plotQC <- function(data) {
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(data, feature1 = "percent.dag", feature2 = "percent.mt")
  plot4 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.dag")
  plot1 + plot2 + plot3 + plot4
}

plotQC(WT_CP)
VlnPlot(WT_CP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.dag"), ncol = 4)

WT_CP <- subset(WT_CP, subset = nFeature_RNA > 200 & nFeature_RNA < 5800 & percent.mt < 6 & percent.dag < 21)

plotQC(KO_CP)
VlnPlot(KO_CP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.dag"), ncol = 4)

KO_CP <- subset(KO_CP, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 7 & percent.dag < 19)

plotQC(WT_LP)
VlnPlot(WT_LP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.dag"), ncol = 4)

WT_LP <- subset(WT_LP, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 9 & percent.dag < 20)

plotQC(KO_LP)
VlnPlot(KO_LP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.dag"), ncol = 4)

KO_LP <- subset(KO_LP, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 9 & percent.dag < 21)

#specify strain/region in a metadata column for each dataset
add_metadata_cols <- function(dataset, strain, region) {
  dataset[["strain"]] <- strain
  dataset[["region"]] <- region
  return(dataset)
}

WT_CP <- add_metadata_cols(WT_CP, "WT", "CP")
WT_LP <- add_metadata_cols(WT_LP, "WT", "LP")
KO_CP <- add_metadata_cols(KO_CP, "KO", "CP")
KO_LP <- add_metadata_cols(KO_LP, "KO", "LP")

#MERGE ALL FOUR DATASETS TOGETHER
all_merged <- merge(x = WT_LP, y = c(WT_CP, KO_CP, KO_LP))
all_merged[["all_merged.percent.mt"]] <- PercentageFeatureSet(all_merged, pattern = "^mt-") 
all_merged <- SCTransform(all_merged, vars.to.regress = c("nCount_RNA", "all_merged.percent.mt"), verbose = TRUE, return.only.var.genes = F)
all_merged <- RunPCA(all_merged, verbose = FALSE)

ElbowPlot(all_merged, ndims = 40)
DimHeatmap(all_merged, dims = 25:42, cells = 500, balanced = TRUE)

#retrieve cell counts
table(Idents(all_merged), all_merged$orig.ident)
# cell counts: KO_CP=2107, KO_LP=1840, WT_CP=7256, WT_LP=4166

#Clustering, no batching needed
all_merged <- RunUMAP(all_merged, dims = 1:30, reduction = "pca")
all_merged <- FindNeighbors(object = all_merged, dims = 1:30, reduction = "pca")
all_merged <- FindClusters(object = all_merged, resolution = 0.2)
DimPlot(all_merged, group.by = "orig.ident", label = FALSE) + DimPlot(all_merged, label = FALSE)

#Finding differentially expressed genes
all_merged.markers <- FindAllMarkers(all_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_merged_sig <- all_merged.markers[which(all_merged.markers$p_val_adj < 0.05), ]
all_merged_sig <- all_merged_sig[!grepl("^mt-", rownames(all_merged_sig)), ]
all_merged_sig <- all_merged_sig[!grepl("^Rp", rownames(all_merged_sig)), ]
write.csv(all_merged_sig, file = "~/all_merged_sig_dims30res02.csv")

#Example of making a heatmap
top30 <- all_merged_sig %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
DoHeatmap(all_merged, cells = WhichCells(all_merged, downsample = 50, seed = 111), features=top30$gene) + NoLegend()

#Example of making featureplots
FeaturePlot(all_merged, features = c("Cebpb","Csf1r","Adgre1","Fcgr1","Cx3cr1","Ccr2","Itgam"), order=TRUE, label=TRUE)

#Example of making violin plots
VlnPlot(all_merged, features = c("C1qa","C1qc","Mafb","Maf","Fcgr1","Adgre1","Cd68","Cd14"), pt.size = F, ncol = 4)

#Subsetting and re-clustering for macrophages and monocytes
MPMo <- subset(all_merged, idents = c(1,2,4,8))
MPMo <- SCTransform(MPMo, vars.to.regress = c("nCount_RNA", "all_merged.percent.mt"), verbose = TRUE, return.only.var.genes = F)
MPMo <- RunPCA(MPMo, verbose = FALSE)
ElbowPlot(MPMo, ndims = 25)
DimHeatmap(MPMo, dims = 1:24, cells = 500, balanced = TRUE)
MPMo <- FindNeighbors(MPMo, dims = 1:10, verbose = FALSE, reduction = "pca")
MPMo <- FindClusters(MPMo, verbose = FALSE, resolution = 0.4)
MPMo <- RunUMAP(MPMo, dims = 1:10, reduction = "pca")
DimPlot(MPMo, label = TRUE)

#Example of retrieving cluster proportions
prop.table(table(Idents(MPMo)))
prop.table(table(Idents(MPMo), MPMo$orig.ident))
table(Idents(MPMo), MPMo$orig.ident)
table(Idents(MPMo))

##Monocle trajectory test
#MANUALLY IMPORTING SEURAT DATA INTO MONOCLE 3
#cell metadata
cells <- MPMo@meta.data[,1:13]

#gene df
genes <- MPMo@assays$RNA@counts@Dimnames[[1]]
gene_df <- data.frame(ids = genes, row.names = genes)
colnames(gene_df) <- "gene_short_name"

#expression matrix
expression_matrix <- MPMo@assays[["RNA"]]@counts

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cells,
                         gene_metadata = gene_df)

#extract highly variable genes from seurat analysis for monocle analysis
HVG <- VariableFeatures(MPMo)

cds <- preprocess_cds(cds, num_dim = 12, use_genes = HVG)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA", umap.n_neighbors = 100, umap.min_dist = 0.5)
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)

p1<-plot_cells(cds, color_cells_by = "partition")
p2<-plot_cells(cds, color_cells_by = "seurat_clusters")
p1+p2
p3<-plot_cells(cds)
p2+p3
#using neighbors=100 min.dist=0.2 res=1e-3 generates 4 partitions with similar #clusters to Seurat(except broke Mos to 2)
#using neighbors=100 min.dist=0.5 res=1e-3 generates 1 partition with same #clusters to Seurat

p1<-plot_cells(cds, cell_size = 0.75, label_cell_groups=FALSE, label_branch_points=FALSE, show_trajectory_graph = FALSE)
p2<-plot_cells(cds, color_cells_by = "seurat_clusters", cell_size = 0.75, label_cell_groups=FALSE, label_branch_points=FALSE, show_trajectory_graph = FALSE)
p3<-plot_cells(cds, color_cells_by="partition", group_cells_by="partition", label_cell_groups=FALSE, cell_size = 0.75, show_trajectory_graph = FALSE)
p1+p2

cds <- learn_graph(cds)
p1<-plot_cells(cds, cell_size = 0.75, label_branch_points = FALSE, label_leaves = FALSE)
p2<-plot_cells(cds, color_cells_by = "seurat_clusters",cell_size = 0.75, label_branch_points = FALSE, label_leaves = FALSE)
p1+p2
plot_cells(cds, genes = c("Ly6c2","Ccr2","Timd4","Cd4"), cell_size = 1, label_branch_points = FALSE, label_leaves = FALSE, label_cell_groups = FALSE, label_roots = FALSE)
plot_cells(cds, genes = c("Timd4","Lyve1","Cd4","Ccr2"), cell_size = 1, label_branch_points = FALSE, label_leaves = FALSE)
plot_cells(cds, color_cells_by = "orig.ident", cell_size = 0.8, graph_label_size = 1)

cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_branch_points=FALSE, label_leaves = FALSE, graph_label_size=1.5, cell_size = 0.75)

###Optional: add UMAP from Seurat
cds@int_colData@listData$reducedDims@listData[["UMAP"]] <- MPMo@reductions[["umap"]]@cell.embeddings
plot_cells(cds)

#visualize specific genes as a function of pseudotime
#choose genes, then subset CDS into a CDS with only those genes
genesSub <- c("C1qa", "Ccr2", "Lyz1", "Lyve1", "Cd4", "Il22ra2")
cds_subset <- cds[rowData(cds)$gene_short_name %in% genesSub]
plot_genes_in_pseudotime(cds_subset, 
                         color_cells_by = "seurat_clusters",
                         min_expr=0.5, cell_size = 1.5)