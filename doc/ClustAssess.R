## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(1001)

## ----setup, fig.width=7, fig.height=4.5---------------------------------------
library(Seurat)
library(ClustAssess)

# we will use the pbmc_small single-cell RNA-seq dataset available via Seurat
DimPlot(pbmc_small, group.by='letter.idents')

## ----pac, fig.width=7, fig.height=4.5-----------------------------------------
# retrieve scaled data for PAC calculation
pbmc.data = GetAssayData(pbmc_small, assay='RNA', slot='scale.data')

# perform consensus clustering
cc.res = consensus_cluster(t(pbmc.data), 
                           k_max=30, 
                           n_reps=100,
                           p_sample=0.8,
                           p_feature=0.8,
                           verbose=TRUE)

# assess the PAC convergence for a few values of k - each curve should 
# have converged to some value
k.plot = c(4,6,8,10)
pac_convergence(cc.res, k.plot)

# visualize the final PAC across k - there seems to be a local maximum at k=7, 
# indicating that 7 clusters leads to a less stable clustering of the data than 
# nearby values of k
pac_landscape(cc.res)

## ----ecs, fig.width=7, fig.height=4.5-----------------------------------------
# first, cluster with Louvain algorithm
pbmc_small = FindClusters(pbmc_small, resolution=0.8, verbose=FALSE)
DimPlot(pbmc_small, group.by='seurat_clusters')

# also cluster with PCA+k-means
pbmc_pca = Embeddings(pbmc_small, 'pca')
pbmc_small@meta.data$kmeans_clusters = kmeans(pbmc_pca, 
                                              centers=3, 
                                              nstart=10, 
                                              iter.max=1e3)$cluster
DimPlot(pbmc_small, group.by='kmeans_clusters')

# where are the clustering results more similar?
pbmc_small@meta.data$ecs = element_sim_elscore(pbmc_small@meta.data$seurat_clusters, 
                                               pbmc_small@meta.data$kmeans_clusters)
FeaturePlot(pbmc_small, 'ecs')
mean(pbmc_small@meta.data$ecs)

## ----jsi, fig.width=7, fig.height=4.5-----------------------------------------
# first, we calculate the markers on the Louvain clustering
Idents(pbmc_small) = pbmc_small@meta.data$seurat_clusters
louvain.markers = FindAllMarkers(pbmc_small, 
                                 logfc.threshold=1, 
                                 min.pct=0.0, 
                                 test.use='roc', 
                                 verbose=FALSE)

# then we get the markers on the k-means clustering
Idents(pbmc_small) = pbmc_small@meta.data$kmeans_clusters
kmeans.markers = FindAllMarkers(pbmc_small, 
                                logfc.threshold=1, 
                                min.pct=0.0, 
                                test.use='roc', 
                                verbose=FALSE)

# next, compare the top 10 markers per cell
pbmc_small@meta.data$marker.gene.jsi = marker_overlap(louvain.markers, 
                                          kmeans.markers, 
                                          pbmc_small@meta.data$seurat_clusters, 
                                          pbmc_small@meta.data$kmeans_clusters, 
                                          n=10, 
                                          rank_by='power')

# which cells have the same markers, regardless of clustering?
FeaturePlot(pbmc_small, 'marker.gene.jsi')
mean(pbmc_small@meta.data$marker.gene.jsi)

## ----frust, fig.width=7, fig.height=4.5---------------------------------------
clustering.list = list()
n.reps = 20
for (i in 1:n.reps){
  # we set nstart=1 and a fairly high iter.max - this should mean that
  # the algorithm converges, and that the variability in final clusterings
  # depends mainly on the random initial cluster assignments
  km.res = kmeans(pbmc_pca, centers=3, nstart=1, iter.max=1e3)$cluster
  clustering.list[[i]] = km.res
}

# now, we calculate the element-wise consistency (aka frustration), which
# performs pair-wise comparisons between all 20 clusterings; the 
# consistency is the average per-cell ECS across all comparisons. The higher
# the consistency, the more consistently is that cell clustered across
# random seeds.
pbmc_small@meta.data$consistency = element_consistency(clustering.list)

# which cells are clustered more consistently?
FeaturePlot(pbmc_small, 'consistency')
mean(pbmc_small@meta.data$consistency)

## ----sessinf------------------------------------------------------------------
sessionInfo()

