# Automatic ClustAssess pipeline. Generating the ClustAssess shiny-app

`ClustAssess` provide the option to run the entire stability pipeline
automatically, by choosing the parameters based on their highest
stability. Another useful feature is the option to create, based on this
object, a shiny app that user can interact with and perform the
assessment of the PhenoGraph configuration and of the clusters.

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
library(SeuratData)
library(ggplot2)
packageVersion("ClustAssess") # should be 1.0.0
#> [1] '1.2.0'
```

Process the PBMC 3k Seurat object similarly to the
`Stability-based parameter assessment` vignette.

``` r

InstallData("pbmc3k")
#> Warning: The following packages are already installed and will not be
#> reinstalled: pbmc3k
data("pbmc3k")
pbmc3k <- UpdateSeuratObject(pbmc3k)
#> Validating object structure
#> Updating object slots
#> Ensuring keys are in the proper structure
#> Warning: Assay RNA changing from Assay to Assay
#> Ensuring keys are in the proper structure
#> Ensuring feature names don't have underscores or pipes
#> Updating slots in RNA
#> Validating object structure for Assay 'RNA'
#> Object representation is consistent with the most current Seurat version
pbmc3k <- PercentageFeatureSet(pbmc3k, pattern = "^MT-", col.name = "percent.mito")
pbmc3k <- PercentageFeatureSet(pbmc3k, pattern = "^RP[SL][[:digit:]]", col.name = "percent.rp")
# remove MT and RP genes
all.index <- seq_len(nrow(pbmc3k))
MT.index <- grep(pattern = "^MT-", x = rownames(pbmc3k), value = FALSE)
RP.index <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(pbmc3k), value = FALSE)
pbmc3k <- pbmc3k[!((all.index %in% MT.index) | (all.index %in% RP.index)), ]
pbmc3k <- subset(pbmc3k, nFeature_RNA < 2000 & nCount_RNA < 2500 & percent.mito < 7 & percent.rp > 7)
pbmc3k <- NormalizeData(pbmc3k, verbose = FALSE)
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

features <- dimnames(pbmc3k@assays$RNA)[[1]]
var_features <- pbmc3k@assays[["RNA"]]@var.features
n_abundant <- 3000
most_abundant_genes <- rownames(pbmc3k@assays$RNA)[order(Matrix::rowSums(pbmc3k@assays$RNA),
    decreasing = TRUE
)]

pbmc3k <- ScaleData(pbmc3k, features = features, verbose = FALSE)

pbmc3k <- RunPCA(pbmc3k,
    npcs = 30,
    approx = FALSE,
    verbose = FALSE,
    features = intersect(most_abundant_genes, pbmc3k@assays$RNA@var.features)
)
pbmc3k <- RunUMAP(pbmc3k,
    reduction = "pca",
    dims = 1:30,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    verbose = FALSE
)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
```

Initialise the parallel backend.

``` r

RhpcBLASctl::blas_set_num_threads(1)
ncores <- 1
if (ncores > 1) {
    my_cluster <- parallel::makeCluster(
        ncores,
        type = "PSOCK"
    )

    doParallel::registerDoParallel(cl = my_cluster)
}
```

Define the feature sets of interest.

``` r

features <- dimnames(pbmc3k@assays$RNA)[[1]]
var_features <- pbmc3k@assays[["RNA"]]@var.features
n_abundant <- 3000
most_abundant_genes <- rownames(pbmc3k@assays$RNA)[order(Matrix::rowSums(pbmc3k@assays$RNA),
    decreasing = TRUE
)]

steps <- seq(from = 500, to = 3000, by = 500)
ma_hv_genes_intersection_sets <- sapply(steps, function(x) intersect(most_abundant_genes[1:x], var_features[1:x]))
ma_hv_genes_intersection <- Reduce(union, ma_hv_genes_intersection_sets)
ma_hv_steps <- sapply(ma_hv_genes_intersection_sets, length)
```

Apply the automatic assessment.

``` r

automm_output <- automatic_stability_assessment(
    expression_matrix = pbmc3k@assays$RNA@scale.data,
    n_repetitions = 10,
    n_neigh_sequence = seq(from = 5, to = 50, by = 5),
    resolution_sequence = seq(from = 0.1, to = 1, by = 0.1),
    features_sets = list(
        "HV" = var_features,
        "MA" = most_abundant_genes[1:3000]
    ),
    steps = list(
        "HV" = steps,
        "MA" = steps
    ),
    n_top_configs = 2,
    umap_arguments = list(
        min_dist = 0.3,
        n_neighbors = 30,
        metric = "cosine"
    ),
    save_temp = FALSE,
    verbose = TRUE
)
#> [2026-08-29 16:24:47.040217] Assessing the stability of the dimensionality reduction
#> [1] "HV"
#> [1] "MA"
#> [2026-08-29 16:27:42.804271] Choosing the top 2
#> [2026-08-29 16:27:43.17902] Assessing the stability of the connected components
#> [2026-08-29 16:27:45.209415] Assessing the stability of the graph construction parameters
#> [2026-08-29 16:30:07.169616] Assessing the stability of the graph clustering method
#> [1] "HV 500"
#> [1] "HV 1000"
#> [1] "MA 2000"
#> [1] "MA 3000"
```

Close the connections opened when using multiple cores.

``` r

if (ncores > 1) {
    parallel::stopCluster(cl = my_cluster)
}
```

Create the shiny app based on the ClustAssess output. You should also
specify either a seurat object or a normalized expression matrix.
*Note*: Please make sure that the directory mentioned in the parameter
`project_folder` is empty / doesn’t exist.

``` r

# generate using a seurat object
write_shiny_app(
    object = pbmc3k,
    assay_name = "RNA",
    clustassess_object = automm_output,
    project_folder = "clustassess_app_dir_seurat",
    shiny_app_title = "PBMC 3k dataset"
)

# generate using a normalized expression matrix
write_shiny_app(
    object = pbmc3k@assays$RNA@data,
    metadata = pbmc3k@meta.data,
    clustassess_object = automm_output,
    project_folder = "clustassess_app_dir_expr",
    shiny_app_title = "PBMC 3k dataset"
)
```

The app can be run using the following command.

``` r

shiny::runApp("clustassess_app_dir_seurat")
```

Session info

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
#>  [1] ggplot2_4.0.3                 ClustAssess_1.2.0            
#>  [3] Seurat_5.5.1                  SeuratObject_5.4.0           
#>  [5] sp_2.2-3                      pbmcMultiome.SeuratData_0.1.4
#>  [7] pbmc3k.SeuratData_3.1.4       SeuratData_0.2.2.9002        
#>  [9] devtools_2.5.2                usethis_3.2.1                
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.5        
#>   [4] spatstat.utils_3.2-4   farver_2.1.2           rmarkdown_2.31        
#>   [7] fs_2.1.0               ragg_1.5.2             vctrs_0.7.3           
#>  [10] ROCR_1.0-12            memoise_2.0.1          spatstat.explore_3.8-2
#>  [13] progress_1.2.3         htmltools_0.5.9        sass_0.4.10           
#>  [16] sctransform_0.4.3      parallelly_1.48.0      KernSmooth_2.23-26    
#>  [19] bslib_0.12.0           htmlwidgets_1.6.4      desc_1.4.3            
#>  [22] ica_1.0-3              plyr_1.8.9             plotly_4.12.1         
#>  [25] zoo_1.9-0              cachem_1.1.0           igraph_2.3.3          
#>  [28] iterators_1.0.14       mime_0.13              lifecycle_1.0.5       
#>  [31] pkgconfig_2.0.3        Matrix_1.7-2           R6_2.6.1              
#>  [34] fastmap_1.2.0          fitdistrplus_1.2-6     future_1.75.0         
#>  [37] shiny_1.14.0           digest_0.6.39          patchwork_1.3.2       
#>  [40] tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7           
#>  [43] pkgload_1.5.3          textshaping_1.0.5      progressr_1.0.0       
#>  [46] spatstat.sparse_3.2-0  httr_1.4.8             polyclip_1.10-7       
#>  [49] abind_1.4-8            compiler_4.4.3         withr_3.0.3           
#>  [52] S7_0.2.2               fastDummies_1.7.6      pkgbuild_1.4.8        
#>  [55] MASS_7.3-64            rappdirs_0.3.4         sessioninfo_1.2.4     
#>  [58] tools_4.4.3            lmtest_0.9-40          otel_0.2.0            
#>  [61] httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3         
#>  [64] glue_1.8.1             nlme_3.1-167           promises_1.5.0        
#>  [67] grid_4.4.3             Rtsne_0.17             cluster_2.1.8         
#>  [70] reshape2_1.4.5         generics_0.1.4         gtable_0.3.6          
#>  [73] spatstat.data_3.1-9    tidyr_1.3.2            hms_1.1.4             
#>  [76] data.table_1.18.6.1    spatstat.geom_3.8-2    RcppAnnoy_0.0.23      
#>  [79] foreach_1.5.2          ggrepel_0.9.8          RANN_2.6.3            
#>  [82] pillar_1.11.1          stringr_1.6.0          spam_2.11-4           
#>  [85] RcppHNSW_0.7.0         later_1.4.8            splines_4.4.3         
#>  [88] dplyr_1.2.1            lattice_0.22-6         survival_3.8-3        
#>  [91] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2          
#>  [94] pbapply_1.7-4          knitr_1.51             gridExtra_2.3.1       
#>  [97] scattermore_1.2        RhpcBLASctl_0.23-42    xfun_0.60             
#> [100] matrixStats_1.5.0      stringi_1.8.9          yaml_2.3.12           
#> [103] evaluate_1.0.5         codetools_0.2-20       tibble_3.3.1          
#> [106] cli_3.6.6              uwot_0.2.4             xtable_1.8-8          
#> [109] reticulate_1.46.0      systemfonts_1.3.2      jquerylib_0.1.4       
#> [112] Rcpp_1.1.2             globals_0.19.1         spatstat.random_3.5-1 
#> [115] png_0.1-9              spatstat.univar_3.2-0  parallel_4.4.3        
#> [118] ellipsis_0.3.3         pkgdown_2.2.1          prettyunits_1.2.0     
#> [121] dotCall64_1.2          listenv_1.0.0          viridisLite_0.4.3     
#> [124] scales_1.4.0           ggridges_0.5.7         purrr_1.2.2           
#> [127] crayon_1.5.3           rlang_1.3.0            cowplot_1.2.0
```
