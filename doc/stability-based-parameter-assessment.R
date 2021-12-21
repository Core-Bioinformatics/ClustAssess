## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(1001)

## ----setup, warning = F, echo = F, message = F--------------------------------
library(Seurat)
library(readr)
library(Matrix)
library(ggplot2)
library(ClustAssess)
library(patchwork)

## ----load---------------------------------------------------------------------
meta.path = 'http://bioinf.stemcells.cam.ac.uk:3838/clustassess_1gmg5bq/cuomo_metadata.csv'
counts.path = 'http://bioinf.stemcells.cam.ac.uk:3838/clustassess_1gmg5bq/all_genes_with_RP_MT.csv.gz'
cell_metadata = read.csv(meta.path, sep=',', header = TRUE, row.names = 1)
cell_metadata$donor = as.factor(cell_metadata$donor)
cell_metadata$day = as.factor(cell_metadata$day)

cuomo.counts = read_csv(counts.path)
gene_names = cuomo.counts$gene_name
cuomo.counts = cuomo.counts[, 1:1880]

cuomo =  CreateSeuratObject(cuomo.counts, row.names = gene_names)

cuomo = PercentageFeatureSet(cuomo, pattern = '^MT-', col.name = 'percent.mito')
cuomo = PercentageFeatureSet(cuomo, pattern = "^RP[SL][[:digit:]]", col.name = 'percent.rp')
one_level = as.factor(rep("cuomo", dim(cuomo)[2]))

cuomo@meta.data$level.ident = one_level
cuomo@meta.data$day.ident = cell_metadata$day
cuomo@meta.data$donor.ident = cell_metadata$donor
Idents(cuomo) = "donor.ident"

rm(cuomo.counts)
gc()

## ----remove_mt_rp-------------------------------------------------------------
# remove MT and RP genes
all.index = 1:nrow(cuomo)
MT.index <- grep(pattern = "^MT-", x = rownames(cuomo), value = FALSE)
RP.index = grep(pattern = "^RP[SL][[:digit:]]", x = rownames(cuomo), value = FALSE)
cuomo = cuomo[!((all.index %in% MT.index) | (all.index %in% RP.index)  ), ]

## ----normalize, warning=F-----------------------------------------------------
cuomo = NormalizeData(cuomo, verbose = F)
cuomo = FindVariableFeatures(cuomo, selection.method = "vst", nfeatures = 3000, verbose = F)

features = dimnames(cuomo@assays$RNA)[[1]]
var_features = cuomo@assays[["RNA"]]@var.features
n_abundant = 3000
most_abundant_genes = rownames(cuomo@assays$RNA)[order(Matrix::rowSums(cuomo@assays$RNA),
                                                       decreasing=TRUE)]

cuomo = ScaleData(cuomo, features = features, verbose = F)

cuomo = RunPCA(cuomo,
               npcs = 30,
               approx = F,
               verbose = F,
               features = intersect(most_abundant_genes, cuomo@assays$RNA@var.features))
cuomo = RunUMAP(cuomo,reduction = "pca",
                dims = 1:30,
                n.neighbors = 30,
                min.dist = 0.3,
                metric = "cosine",
                verbose = F)
raw_umap = cuomo@reductions$umap@cell.embeddings

## ----qc, fig.width=7, fig.height=7, include = T, echo = F---------------------
wrap_plots(lapply(c(DimPlot(cuomo, group.by=c('day.ident', 'donor.ident'), combine = F),
                    FeaturePlot(cuomo, features = c('nFeature_RNA', 'nCount_RNA','percent.mito', 'percent.rp'), combine = F)),
                  function(x) { x + labs(title = element_blank()) + theme(legend.key.size = unit(0.3, "cm"),
                                                                          legend.key.width = unit(0.1,"cm"),
                                                                          legend.text = element_text(size = 7),
                                                                          axis.text = element_text(size = 7),
                                                                          axis.title = element_text(size = 9)) })) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(nrow = 3)

## ----nrep---------------------------------------------------------------------
n_repetitions = 30

## ----ncores-------------------------------------------------------------------
n_cores = 1

## ----feature_sets-------------------------------------------------------------
ma_hv_genes_intersection = intersect(most_abundant_genes[1:3000], var_features)
steps = seq(from = 500, to = 3000, by = 500)
ma_hv_steps = sapply(steps, function(x) { length(intersect(most_abundant_genes[1:x], var_features))})

## ---- include = FALSE---------------------------------------------------------
start_time_feature_stability = Sys.time()

## ----feature_stability, warning=FALSE-----------------------------------------
pca_feature_stability_object = c(get_feature_stability(data_matrix = cuomo@assays[["RNA"]]@scale.data,
                                                       feature_set = most_abundant_genes,
                                                       steps = steps,
                                                       n_repetitions = n_repetitions,
                                                       feature_type = "MA",
                                                       graph_reduction_type = "PCA",
                                                       npcs = 30,
                                                       min_dist = 0.3,
                                                       n_neighbors = 30,
                                                       metric = "cosine",
                                                       ncores = n_cores,
                                                       ecs_thresh = 1,
                                                       algorithm = 1),
                                 get_feature_stability(data_matrix = cuomo@assays[["RNA"]]@scale.data,
                                                       feature_set = var_features,
                                                       steps = steps,
                                                       n_repetitions = n_repetitions,
                                                       feature_type = "HV",
                                                       graph_reduction_type = "PCA",
                                                       npcs = 30,
                                                       min_dist = 0.3,
                                                       n_neighbors = 30,
                                                       metric = "cosine",
                                                       ncores = n_cores,
                                                       ecs_thresh = 1,
                                                       algorithm = 1),
                                 get_feature_stability(data_matrix = cuomo@assays[["RNA"]]@scale.data,
                                                       feature_set = ma_hv_genes_intersection,
                                                       steps = ma_hv_steps,
                                                       n_repetitions = n_repetitions,
                                                       feature_type = "MA_HV",
                                                       graph_reduction_type = "PCA",
                                                       npcs = 30,
                                                       min_dist = 0.3,
                                                       n_neighbors = 30,
                                                       metric = "cosine",
                                                       ncores = n_cores,
                                                       ecs_thresh = 1,
                                                       algorithm = 1))

## ----stab_boxplot, fig.width=7, fig.height=5----------------------------------
plot_feature_stability_boxplot(pca_feature_stability_object, text_size  = 2.5) +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0))

## ----stab_ecs_inc, fig.width=7, fig.height=5----------------------------------
plot_feature_stability_ecs_incremental(pca_feature_stability_object, dodge_width = 1, text_size = 2) +
    theme(legend.position = c(1,0),
        legend.justification = c(1,0))

## ----stab_mb, fig.width=7, fig.height = 15------------------------------------
plot_feature_stability_mb_facet(pca_feature_stability_object, text_size = 3)

## ----stab_ecs_fa, fig.width=7, fig.height = 15--------------------------------
plot_feature_stability_ecs_facet(pca_feature_stability_object)

## ---- include = F-------------------------------------------------------------
stop_time_feature_stability = Sys.time()

## ----update_object------------------------------------------------------------
cuomo = cuomo[var_features, ]

cuomo = RunPCA(cuomo,
               npcs = 30,
               approx = F,
               verbose = F)
cuomo = RunUMAP(cuomo,reduction = "pca",
                dims = 1:30,
                n.neighbors = 30,
                min.dist = 0.3,
                metric = "cosine",
                verbose = F)

## ---- include = FALSE---------------------------------------------------------
start_time_nn_conn = Sys.time()

## ----nn_conn_comps------------------------------------------------------------
nn_conn_comps_object = c(get_nn_conn_comps(object = cuomo@reductions$pca@cell.embeddings,
                                           n_neigh_sequence = c(c(1,2,3,4), seq(from = 5, to = 30, by = 5)),
                                           n_repetitions = n_repetitions,
                                           graph_reduction_type = "UMAP",
                                           ncores = n_cores,
                                           min_dist = 0.3,
                                           n_neighbors = 30,
                                           metric = "cosine"),
                         get_nn_conn_comps(object = cuomo@assays[["RNA"]]@scale.data,
                                           n_neigh_sequence = c(c(1,2,3,4), seq(from = 5, to = 30, by = 5)),
                                           n_repetitions = n_repetitions,
                                           graph_reduction_type = "PCA",
                                           ncores = n_cores,
                                           nv = 30))

## ----comps_evo, fig.width=7, fig.height=5-------------------------------------
plot_connected_comps_evolution(nn_conn_comps_object)

## ---- include = F-------------------------------------------------------------
stop_time_nn_conn = Sys.time()

## ---- include = FALSE---------------------------------------------------------
start_time_nn_importance = Sys.time()

## ----nn_importance------------------------------------------------------------
nn_importance_object = mapply(c,
                              get_nn_importance(object = cuomo@assays[["RNA"]]@scale.data,
                                         n_neigh_sequence = seq(from = 5, to = 30, by = 5),
                                         n_repetitions = n_repetitions,
                                         graph_reduction_type = "PCA",
                                         ecs_thresh = 1,
                                         ncores = n_cores,
                                         algorithm = 1,
                                         nv = 30),
                              get_nn_importance(object = cuomo@reductions$pca@cell.embeddings,
                                         n_neigh_sequence = seq(from = 5, to = 30, by = 5),
                                         n_repetitions = n_repetitions,
                                         graph_reduction_type = "UMAP",
                                         ecs_thresh = 1,
                                         ncores = n_cores,
                                         algorithm = 1,
                                         min_dist = 0.3,
                                         n_neighbors = 30,
                                         metric = "cosine"),
                              SIMPLIFY = FALSE
)

## ----n_k_corr, fig.width=7, fig.height=5--------------------------------------
plot_n_neigh_k_correspondence(nn_importance_object)

## ----n_neigh_ecs, fig.width=7, fig.height=5-----------------------------------
plot_n_neigh_ecs(nn_importance_object)

## ---- include = F-------------------------------------------------------------
stop_time_nn_importance = Sys.time()

## ----clustering_diff----------------------------------------------------------
cuomo = RunPCA(cuomo,
               npcs = 30,
               approx = F,
               verbose = F)
cuomo = RunUMAP(cuomo,reduction = "pca",
                dims = 1:30,
                n.neighbors = 30,
                min.dist = 0.3,
                metric = "cosine",
                verbose = F)
adj_matrix = FindNeighbors(cuomo@reductions$umap@cell.embeddings, k.param = 25, nn.method = "rann", verbose = F)$snn

start_time_clustering_method = Sys.time()
clustering_diff_obj = get_clustering_difference(graph_adjacency_matrix = adj_matrix,
                                                resolution = seq(from = 0.5, to = 1, by = 0.1),
                                                n_repetitions = n_repetitions,
                                                ecs_thresh = 1,
                                                ncores = n_cores,
                                                algorithm = 1:4)

## ----diff_boxplot, fig.width=7, fig.height=5----------------------------------
plot_clustering_difference_boxplot(clustering_diff_obj)

## ----diff_facet, fig.width=7, fig.height=10-----------------------------------
plot_clustering_difference_facet(clustering_diff_obj, cuomo@reductions$umap@cell.embeddings)

## ---- include = F-------------------------------------------------------------
stop_time_clustering_method = Sys.time()

## ---- include = F-------------------------------------------------------------
start_time_resolution = Sys.time()

## ----resolution_importance----------------------------------------------------
resolution_gridsearch = get_resolution_importance(embedding = cuomo@reductions$umap@cell.embeddings,
                                                  resolution = seq(from = 0.5, to = 1, by = 0.1),
                                                  n_neigh = 30,
                                                  n_repetitions = n_repetitions,
                                                  clustering_method = c(1,2),
                                                  graph_type = 2,
                                                  ecs_thresh = 1,
                                                  ncores = n_cores)

## ----k_res_corr_1, fig.width=7, fig.height=5----------------------------------
plot_k_resolution_corresp(resolution_gridsearch) +
  ggtitle("resolution - k correspondence with ecs threshold = 1")

## ----k_n_part_1, fig.width=7, fig.height=5------------------------------------
plot_k_n_partitions(resolution_gridsearch) + ggtitle("k - # partitions correspondence with ecs threshold = 1")

## ----merge_partitions---------------------------------------------------------
resolution_gridsearch_thresh_99 = merge_partitions(resolution_gridsearch,
                                                   ecs_thresh = 0.99,
                                                   ncores = n_cores)

## ----k_res_corr_99, fig.width=7, fig.height=5---------------------------------
plot_k_resolution_corresp(resolution_gridsearch_thresh_99) +
  ggtitle("resolution - k correspondence with ecs threshold = 0.99")

## ----k_n_part_99, fig.width=7, fig.height=5-----------------------------------
plot_k_n_partitions(resolution_gridsearch_thresh_99) +
  ggtitle("k - # partitions correspondence with ecs threshold = 0.99")

## ---- include = F-------------------------------------------------------------
stop_time_resolution = Sys.time()

## -----------------------------------------------------------------------------
paste("Feature stability methods runtime:",
      format(as.numeric(stop_time_feature_stability - start_time_feature_stability,
                        units = "mins")), "minutes")
paste("NN - # connected components methods runtime:",
      format(as.numeric(stop_time_nn_conn - start_time_nn_conn,
                        units = "mins")), "minutes")
paste("NN importance methods runtime:",
      format(as.numeric(stop_time_nn_importance - start_time_nn_importance,
                        units = "mins")), "minutes")
paste("Clustering importance methods runtime:",
      format(as.numeric(stop_time_clustering_method - start_time_clustering_method,
                        units = "mins")), "minutes")
paste("Resolution gridsearch methods runtime:",
      format(as.numeric(stop_time_resolution - start_time_resolution,
                        units = "mins")), "minutes")

## -----------------------------------------------------------------------------
sessionInfo()

