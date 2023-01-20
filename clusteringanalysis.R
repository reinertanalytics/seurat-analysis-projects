library(BiocManager)
library(Seurat)
library(aricode)
library(dplyr)
library(dbscan)
library(cluster)
library(splatter)
library(ggfortify)
library(fpc)

# read data
data <- Read10X(data.dir = "../documents/hg19")

# read count matrix into Seurat object
seudata <- CreateSeuratObject(counts = data, project = "P2", min.cells = 3, min.features = 200)

# trim data
seudata[["percent.mt"]] <- PercentageFeatureSet(seudata, pattern = "^MT-")
seudata <- subset(seudata, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize data
seudata <- NormalizeData(seudata)

# identify most variable features
seuvars <- FindVariableFeatures(seudata, selection.method = 'vst', nfeatures = 2000)
genes <- rownames(seuvars)

# scale data again
seuvars <- ScaleData(seuvars, features = genes)

# run PCA analysis
seuvars <- RunPCA(seuvars, features = VariableFeatures(object = seuvars))

# review
DimPlot(seuvars, reduction = "pca")

chickenleg <- subset(seuclust, idents = c(0, 2, 4, 6))
chickenleg <- ScaleData(chickenleg, vars.to.regress = c("percent.mt"))
chickenleg <- RunPCA(chickenleg)
DimPlot(chickenleg, reduction = "pca")

# review and identify appropriate dimensions
plot <- ElbowPlot(seuvars)
plot

# run clustering after identifying 9 as an appropriate number of dimensions
seuclust <- FindNeighbors(seuvars, dims = 1:9)
seuclust <- FindClusters(seuclust, resolution = 0.5)

# reduce with UMAP and plot
seuclust <- RunUMAP(seuclust, dims = 1:9)
DimPlot(seuclust, reduction = "umap")

# identify two superclusters and subset the first, which is shaped roughly like a chicken leg
chickenleg <- subset(seuclust, idents = c(0, 2, 4, 6))
chickenleg <- ScaleData(chickenleg, vars.to.regress = c("percent.mt"))
chickenleg <- RunPCA(chickenleg)
DimPlot(chickenleg, reduction = "pca")

# this number of dimensions significantly changes the shape
# after some experimentation, four dimensions fixes the shape and allows deeper analysis

chickenleg4 <- FindNeighbors(chickenleg, dims = 1:4)
chickenleg4 <- FindClusters(chickenleg4, resolution = 0.5)
chickenleg4 <- RunUMAP(chickenleg4, dims = 1:4)
DimPlot(chickenleg4, reduction = "umap")

# here in the further processed data the issue becomes a bit more apparent. 
# the top and the bottom of the image, which looks somewhat like a chicken drumstick now
# represented here as clusters 0, 1, and 5, are fairly similar to the original seurat clusters 0, 2, and 6
# cluster 4 here, originally cluster 4, is significantly reduced at the expense of two new clusters, 2 and 3
# in combination with the previous plot where cluster 4 gained share, this indicates uncertainty in the clustering
# the borders between original clusters 0, 2, and 4, as well as potentially implying that additional clusters were necessary

# for the other supercluster, which somewhat resembles a dolphin
dolphin <- subset(seuclust, idents = c(1, 5, 7))
dolphin <- ScaleData(dolphin, vars.to.regress = c("percent.mt"))
dolphin <- RunPCA(dolphin)
DimPlot(dolphin, reduction = "pca")

# Initial subclustering allows cluster 5 to grow at the expense of cluster 1, but not by much.
dolphin5 <- FindNeighbors(dolphin, dims = 1:5)
dolphin5 <- FindClusters(dolphin3, resolution = 0.5)
dolphin5 <- RunUMAP(dolphin3, dims = 1:5)
DimPlot(dolphin3, reduction = "umap")

# while the dolphin shape has been changed around, the overall structure is still mostly the same. 
# at this level, though, the clustering has really weakened
# portions of clusters are mixed with each other and even completely disconnected from their main groupings
# here is some green near the bottom with the pink, and on the left side the blue is bleeding into the red

# here we will perform cluster quality analysis through Seurat::FindMarkers
# analysing the FindMarkers results and plotting under FeaturePlot to see expression relative to different clusters
# features that are assigned to clusters by Seurat as predictive but are too widely shared likely undermine quality

cluster0 <- FindMarkers(seuclust, ident.1 = 0, min.pct = 0.25)
head(cluster0, n = 10)

cluster1 <- FindMarkers(seuclust, ident.1 = 1, min.pct = 0.25)
head(cluster1, n = 10)

cluster2 <- FindMarkers(seuclust, ident.1 = 2, min.pct = 0.25)
head(cluster2, n = 10)

# FeaturePlot where row 1 represents the top markers for cluster 0, row 2 for cluster 1, and row 3 for cluster 2
FeaturePlot(seuclust, features = c("RPS12", "RPS27", "RPS6", "S100A9", "S100A8", "FCN1", "IL32", "LTB", "IL7R"))

# while Clusters 1 and 2 show reasonable differentiation, cluster 0's top markers appear in basically the entire plot
# this indicates they are a poor basis for a cluster and that cluster 0's quality may be poor and probably meaningless
# given that cluster 0's quality was called into question in the chicken leg analysis, as well, cluster 0's quality probably warrants further investigation
# there is a reasonable chance that this is due to the clustering algorithm
# if the initial point is randomly placed and the algorithm attempts to build a starting cluster around it,
# it may be choosing certain primary markers that are actually shared across all clusters
