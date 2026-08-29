# Evaluating single-cell clustering with ClustAssess

In this vignette we will illustrate several of the tools available in
ClustAssess on a small single-cell RNA-seq dataset.

``` r

library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 'SeuratObject' was built under R 4.4.0 but the current version is
#> 4.4.3; it is recomended that you reinstall 'SeuratObject' as the ABI
#> for R may have changed
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
library(ClustAssess)
library(ggplot2)
theme_set(theme_classic())

data("pbmc_small")

# we will use the pbmc_small single-cell RNA-seq dataset available via Seurat
DimPlot(pbmc_small, group.by = "letter.idents")
```

![](ClustAssess_files/figure-html/setup-1.png)

## Proportion of ambiguously clustered pairs

The proportion of ambiguously clustered pairs (PAC) uses consensus
clustering to infer the optimal number of clusters for the data, by
observing how variably pairs of elements cluster together. The lower the
PAC, the more stable the clustering. PAC was originally presented in
<https://doi.org/10.1038/srep06207>.

``` r

# retrieve scaled data for PAC calculation
pbmc.data <- GetAssayData(pbmc_small, assay = "RNA", layer = "scale.data")

# perform consensus clustering
cc.res <- consensus_cluster(t(pbmc.data),
    k_max = 30,
    n_reps = 100,
    p_sample = 0.8,
    p_feature = 0.8,
    verbose = TRUE
)
#> Calculating consensus clustering

# assess the PAC convergence for a few values of k - each curve should
# have converged to some value
k.plot <- c(4, 6, 8, 10)
pac_convergence(cc.res, k.plot)
```

![](ClustAssess_files/figure-html/pac-1.png)

``` r


# visualize the final PAC across k - there seems to be a local maximum at k=7,
# indicating that 7 clusters leads to a less stable clustering of the data than
# nearby values of k
pac_landscape(cc.res)
```

![](ClustAssess_files/figure-html/pac-2.png)

## Element-centric clustering comparison

We compare the similarity of clustering results between methods (in this
case, between Louvain community detection and k-means) using
Element-Centric Similarity (ECS), which quantifies the clustering
similarity per cell. Higher ECS indicates that an observation is
clustered similarly across methods. ECS was introduced in
<https://doi.org/10.1038/s41598-019-44892-y>.

``` r

# first, cluster with Louvain algorithm
pbmc_small <- FindClusters(pbmc_small, resolution = 0.8, verbose = FALSE)
DimPlot(pbmc_small, group.by = "seurat_clusters")
```

![](ClustAssess_files/figure-html/ecs-1.png)

``` r


# also cluster with PCA+k-means
pbmc_pca <- Embeddings(pbmc_small, "pca")
pbmc_small@meta.data$kmeans_clusters <- kmeans(pbmc_pca,
    centers = 3,
    nstart = 10,
    iter.max = 1e3
)$cluster
DimPlot(pbmc_small, group.by = "kmeans_clusters")
```

![](ClustAssess_files/figure-html/ecs-2.png)

``` r


# where are the clustering results more similar?
pbmc_small@meta.data$ecs <- element_sim_elscore(
    pbmc_small@meta.data$seurat_clusters,
    pbmc_small@meta.data$kmeans_clusters
)
suppressMessages(FeaturePlot(pbmc_small, "ecs") + scale_colour_viridis_c())
```

![](ClustAssess_files/figure-html/ecs-3.png)

``` r

mean(pbmc_small@meta.data$ecs)
#> [1] 0.6632927
```

## Jaccard similarity of cluster markers

As a common step in computational single-cell RNA-seq analyses,
discriminative marker genes are identified for each cluster. These
markers are then used to infer the cell type of the respective cluster.
Here, we compare the markers obtained for each clustering method to ask:
how similarly would each cell be interpreted across clustering methods?
We compare the markers per cell using the Jaccard similarity (defined as
the size of the intersect divided by the size of the overlap) of the
marker gene lists. The higher the JSI, the more similar are the marker
genes for that cell.

``` r

# first, we calculate the markers on the Louvain clustering
Idents(pbmc_small) <- pbmc_small@meta.data$seurat_clusters
louvain.markers <- FindAllMarkers(pbmc_small,
    logfc.threshold = 1,
    min.pct = 0.0,
    test.use = "roc",
    verbose = FALSE
)

# then we get the markers on the k-means clustering
Idents(pbmc_small) <- pbmc_small@meta.data$kmeans_clusters
kmeans.markers <- FindAllMarkers(pbmc_small,
    logfc.threshold = 1,
    min.pct = 0.0,
    test.use = "roc",
    verbose = FALSE
)

# next, compare the top 10 markers per cell
pbmc_small@meta.data$marker.gene.jsi <- marker_overlap(louvain.markers,
    kmeans.markers,
    pbmc_small@meta.data$seurat_clusters,
    pbmc_small@meta.data$kmeans_clusters,
    n = 10,
    rank_by = "power"
)

# which cells have the same markers, regardless of clustering?
suppressMessages(FeaturePlot(pbmc_small, "marker.gene.jsi") + scale_colour_viridis_c())
```

![](ClustAssess_files/figure-html/jsi-1.png)

``` r

mean(pbmc_small@meta.data$marker.gene.jsi)
#> [1] 0.5738636
```

## Element-wise consistency

How consistently are cells clustered by k-means? We will rerun k-means
clustering 20 times to investigate.

``` r

clustering.list <- list()
n.reps <- 20
for (i in 1:n.reps) {
    # we set nstart=1 and a fairly high iter.max - this should mean that
    # the algorithm converges, and that the variability in final clusterings
    # depends mainly on the random initial cluster assignments
    km.res <- kmeans(pbmc_pca, centers = 3, nstart = 1, iter.max = 1e3)$cluster
    clustering.list[[i]] <- km.res
}

# now, we calculate the element-wise consistency (aka frustration), which
# performs pair-wise comparisons between all 20 clusterings; the
# consistency is the average per-cell ECS across all comparisons. The higher
# the consistency, the more consistently is that cell clustered across
# random seeds.
pbmc_small@meta.data$consistency <- element_consistency(clustering.list)

# which cells are clustered more consistently?
suppressMessages(FeaturePlot(pbmc_small, "consistency") + scale_colour_viridis_c())
```

![](ClustAssess_files/figure-html/frust-1.png)

``` r

mean(pbmc_small@meta.data$consistency)
#> [1] 0.6712922
```

## Session info

``` r

sessionInfo()
#> R version 4.4.3 (2025-02-28)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_4.0.3      ClustAssess_1.2.0  Seurat_5.5.1       SeuratObject_5.4.0
#> [5] sp_2.2-3          
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.5        
#>   [4] spatstat.utils_3.2-4   farver_2.1.2           rmarkdown_2.31        
#>   [7] fs_2.1.0               ragg_1.5.2             vctrs_0.7.3           
#>  [10] ROCR_1.0-12            spatstat.explore_3.8-2 htmltools_0.5.9       
#>  [13] progress_1.2.3         sass_0.4.10            sctransform_0.4.3     
#>  [16] parallelly_1.48.0      KernSmooth_2.23-26     bslib_0.12.0          
#>  [19] htmlwidgets_1.6.4      desc_1.4.3             ica_1.0-3             
#>  [22] plyr_1.8.9             plotly_4.12.1          zoo_1.9-0             
#>  [25] cachem_1.1.0           igraph_2.3.3           mime_0.13             
#>  [28] lifecycle_1.0.5        iterators_1.0.14       pkgconfig_2.0.3       
#>  [31] Matrix_1.7-2           R6_2.6.1               fastmap_1.2.0         
#>  [34] fitdistrplus_1.2-6     future_1.75.0          shiny_1.14.0          
#>  [37] digest_0.6.39          patchwork_1.3.2        tensor_1.5.1          
#>  [40] RSpectra_0.16-2        irlba_2.3.7            textshaping_1.0.5     
#>  [43] labeling_0.4.3         progressr_1.0.0        spatstat.sparse_3.2-0 
#>  [46] httr_1.4.8             polyclip_1.10-7        abind_1.4-8           
#>  [49] compiler_4.4.3         withr_3.0.3            S7_0.2.2              
#>  [52] fastDummies_1.7.6      MASS_7.3-64            tools_4.4.3           
#>  [55] lmtest_0.9-40          otel_0.2.0             httpuv_1.6.17         
#>  [58] future.apply_1.20.2    goftest_1.2-3          glue_1.8.1            
#>  [61] nlme_3.1-167           promises_1.5.0         grid_4.4.3            
#>  [64] Rtsne_0.17             cluster_2.1.8          reshape2_1.4.5        
#>  [67] generics_0.1.4         gtable_0.3.6           spatstat.data_3.1-9   
#>  [70] tidyr_1.3.2            data.table_1.18.6.1    hms_1.1.4             
#>  [73] spatstat.geom_3.8-2    RcppAnnoy_0.0.23       ggrepel_0.9.8         
#>  [76] RANN_2.6.3             foreach_1.5.2          pillar_1.11.1         
#>  [79] stringr_1.6.0          spam_2.11-4            RcppHNSW_0.7.0        
#>  [82] later_1.4.8            splines_4.4.3          dplyr_1.2.1           
#>  [85] lattice_0.22-6         survival_3.8-3         deldir_2.0-4          
#>  [88] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4         
#>  [91] knitr_1.51             gridExtra_2.3.1        scattermore_1.2       
#>  [94] xfun_0.60              matrixStats_1.5.0      stringi_1.8.9         
#>  [97] yaml_2.3.12            evaluate_1.0.5         codetools_0.2-20      
#> [100] tibble_3.3.1           cli_3.6.6              uwot_0.2.4            
#> [103] xtable_1.8-8           reticulate_1.46.0      systemfonts_1.3.2     
#> [106] jquerylib_0.1.4        Rcpp_1.1.2             globals_0.19.1        
#> [109] spatstat.random_3.5-1  png_0.1-9              fastcluster_1.3.0     
#> [112] spatstat.univar_3.2-0  parallel_4.4.3         pkgdown_2.2.1         
#> [115] prettyunits_1.2.0      dotCall64_1.2          listenv_1.0.0         
#> [118] viridisLite_0.4.3      scales_1.4.0           ggridges_0.5.7        
#> [121] purrr_1.2.2            crayon_1.5.3           rlang_1.3.0           
#> [124] cowplot_1.2.0
```
