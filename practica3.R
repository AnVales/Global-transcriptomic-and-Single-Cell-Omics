# PRACTICE 3 EXERCISES #
if(!require(matlib)){install.packages("matlib")} 
if(!require(rlang)){install.packages("rlang")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(openxlsx)){install.packages("Seurat")}
if(!require(openxlsx)){install.packages("dplyr")}
library(matlib)
library(rlang)
library(openxlsx)
library(Seurat)
library(Matrix)
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest")

BiocManager::install("igraph")


install.packages("remotes")

remotes::install_version("SDMTools", "1.1-221")

install.packages('devtools')
devtools::install_version(package = 'Seurat', version = package_version('3.1.1'))

library(Seurat)

# colors
newcolors<-colorRampPalette(colors=c("gray48","skyblue","yellow"))(256)


# PART 1 #
sco_data <- read.table("./SCOplanaria.txt", row.names=1) #path
# dim(sco_data) -> 24883  1871
sco_matrix = as.matrix(sco_data)
sco_sparse_matrix = Matrix(data=sco_matrix, sparse=TRUE)

# PART 2 #
# using the following parameters:

# Keep genes expressed in at least 3 cells and cells with at least 200 genes detected.
seurat_sco = CreateSeuratObject(counts=sco_sparse_matrix, project="SCOP", min.cells=3, min.features=200)

# Assay data with 18561 features for 1301 cells + First 10 features (dd-Smed-v6-10001-0, ...)
seurat_sco[["RNA"]] 

# 18561  1301
dim(seurat_sco[["RNA"]])

# shows first 5 rows and columns data
seurat_sco[["RNA"]]@counts[1:5,1:5]

# Keep cells with 200- 2500 number of detected genes
seurat_sco_subset = subset(seurat_sco, subset=nFeature_RNA>200&nFeature_RNA<2500)

# Use LogNormalize method with a scale factor of 10000
seurat_norm = NormalizeData(seurat_sco_subset, normalization.method="LogNormalize", scale.factor=10000)

# Visualize effects of normalizalition in total cell expression
Normalized.Reads <- colSums(seurat_norm[["RNA"]]@data)
seurat_norm <- AddMetaData (seurat_norm,Normalized.Reads,col.name="Normalized.Reads")
VlnPlot(object = seurat_norm, features= c("nCount_RNA", "Normalized.Reads"), ncol = 2)

# Select 3000 variable features
seurat_norm = FindVariableFeatures(seurat_norm, selection.method="vst", nfeatures=3000)
p2 = LabelPoints(plot=VariableFeaturePlot(seurat_norm), 
                 points=head(VariableFeatures(seurat_norm),2), repel=TRUE) + NoLegend()
p2

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_norm), 10)

# plot variable features with and without labels for 10 most variable genes
plot1 <- VariableFeaturePlot(seurat_norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# Regress out variability for number of genes in each cell.
seurat_norm = ScaleData(seurat_norm, features=rownames(seurat_norm))

# Use PCs 1-5 and a resolution value of 0.6, t-SNE 
seurat_pca <- RunPCA(seurat_norm, features = VariableFeatures(object = seurat_norm))

# Genes are grouped based on behaviour and can be visualize in a heatmap (maximum 15 genes per group)
DimHeatmap(seurat_pca, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(seurat_pca)
DimHeatmap(seurat_pca, dims = 1:5, cells = 500, balanced = TRUE)

# Visualize PCA scores for components 1 and 2
VizDimLoadings(seurat_pca, dims = 1:5, reduction = "pca")

# Plot Cells Based on PCA1 and 2 components
DimPlot(seurat_pca, reduction = "pca")


# Seurat cluster cells based on PCA scores
JackStraw_seurat <- JackStraw(seurat_pca, num.replicate = 100)
JackStraw_seurat <- ScoreJackStraw(JackStraw_seurat, dims = 1:15)

# Significant components are those who have a  srong enrichment of low p-values
JackStrawPlot(JackStraw_seurat, dims = 1:15)
windows()
ElbowPlot(JackStraw_seurat)

# t-SNE plot
seurat_pca_plot = FindNeighbors(seurat_pca, dims=1:5)
seurat_pca_plot = FindClusters(seurat_pca_plot, resolution=0.75)
seurat_pca_plot <- RunTSNE(seurat_pca_plot, dims = 1:5)
windows(1)
DimPlot(seurat_pca_plot, reduction = "tsne")

# PART 3 #
# Extract a table of the 5 top biomarkers for each of your clusters.
# Use these biomarkers to generate a heatmap.
# Comment the results

# find markers for every cluster compared to all remaining cells, report only the positive ones
biomarkers = FindAllMarkers(seurat_pca_plot, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
biomarkers_top5 = biomarkers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)

#Visualize as heatmap
DoHeatmap(seurat_pca_plot, features=biomarkers_top5$gene) + NoLegend()

# PART 4 #

# Here we show a table with biomarkers for several cell identities
# Using VlnPlot() and FeaturePlot() functions, show the expression distribution of these markers along your clusters
#  With this information, identify and rename your clusters.

# Gene name / Cell identity
# dd-Smed-v6-61-0 / Early epidermal progenitors
# dd-Smed-v6-2178-0 / Late epidermal progenitors
# dd-Smed-v6-298-0 / Epidermis
# dd-Smed-v6-1410-0 / Muscle progenitors
# dd-Smed-v6-702-0 / Muscle body
# dd-Smed-v6-2548-0 / Neural progenitors
# dd-Smed-v6-9977-0 / GABA neurons
# dd-Smed-v6-48-0 / Phagocytes
# dd-Smed-v6-175-0 / Parenchymal cells
# dd-Smed-v6-1161-1 / Pigment

markers_id = list("dd-Smed-v6-61-0"="Early epidermal progenitors", "dd-Smed-v6-2178-0"="Late epidermal progenitors",
                  "dd-Smed-v6-298-0"="Epidermis", "dd-Smed-v6-1410-0"="Muscle progenitors",
                  "dd-Smed-v6-702-0"="Muscle body", "dd-Smed-v6-2548-0"="Neural progenitors",
                  "dd-Smed-v6-9977-0"="GABA neurons", "dd-Smed-v6-48-0"="Phagocytes",
                  "dd-Smed-v6-175-0"="Parenchymal cells", "dd-Smed-v6-1161-1"="Pigment")

#seurat_pca_plot = names(markers_id)

VlnPlot(seurat_pca_plot, features = c("dd-Smed-v6-61-0", "dd-Smed-v6-2178-0",
                                      "dd-Smed-v6-298-0","dd-Smed-v6-1410-0",
                                      "dd-Smed-v6-702-0", "dd-Smed-v6-2548-0",
                                      "dd-Smed-v6-9977-0", "dd-Smed-v6-48-0",
                                      "dd-Smed-v6-175-0", "dd-Smed-v6-1161-1",
                                      "dd-Smed-v6-1999-0"))     



FeaturePlot(seurat_pca_plot, features = c("dd-Smed-v6-61-0", "dd-Smed-v6-2178-0",
                                          "dd-Smed-v6-298-0","dd-Smed-v6-1410-0",
                                          "dd-Smed-v6-702-0", "dd-Smed-v6-2548-0",
                                          "dd-Smed-v6-9977-0", "dd-Smed-v6-48-0",
                                          "dd-Smed-v6-175-0", "dd-Smed-v6-1161-1",
                                          "dd-Smed-v6-1999-0"))
DimPlot(seurat_pca_plot, reduction = "tsne")

# dd-Smed-v6-1999-0 is a neoblast (stem) marker gene. 
# Show its expression distribution with VlnPlot() and FeaturePlot() and explain the result.

VlnPlot(seurat_pca_plot, features = c("dd-Smed-v6-1999-0"))                  
FeaturePlot(seurat_pca_plot, features = c("dd-Smed-v6-1999-0"))
DimPlot(seurat_pca_plot, reduction = "tsne")

# Compare your results with this lineage tree reconstruction of planarian cell types (Plass et al., 2018).
# What cell type do you think cluster 0 (the central cluster) is? 
# Can all of the planarian cell types be found in your plot? Why? 
# rename each cluster #

new.cluster.ids <- c("Neoblast", "neural progenitors",
                     "early epidermal progenitors", "epidermis", 
                     "GABA neurons", "parenchymal cells",
                     "muscle body", "late epidermal progenitors",
                     "mucle progenitors", "pigment", 
                     "phagocytes")

names(new.cluster.ids) <- levels(seurat_pca_plot)
seurat_pca_names <- RenameIdents(seurat_pca_plot, new.cluster.ids)
DimPlot(seurat_pca_names, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()


