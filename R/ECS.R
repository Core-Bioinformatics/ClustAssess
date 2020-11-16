#' Title
#'
#' @param clustering1
#' @param clustering2
#' @param alpha
#' @param r
#' @param r2
#' @param rescale_path_type
#' @param ppr_implementation
#'
#' @return
#' @export
#'
#' @examples
element.sim = function(clustering1, clustering2, alpha=0.9, r=1., r2=NULL,
                       rescale_path_type='max', ppr_implementation='prpack'){
  results = element.sim.elscore(clustering1, clustering2, alpha=alpha,
                                     r=r, r2=r2,
                                     rescale_path_type=rescale_path_type,
                                     ppr_implementation=ppr_implementation)
  elementScores = results$scores
  return(mean(elementScores))
}

#' Title
#'
#' @param clustering1
#' @param clustering2
#' @param alpha
#' @param r
#' @param r2
#' @param rescale_path_type
#' @param relabeled_elements
#' @param ppr_implementation
#'
#' @return
#' @export
#'
#' @examples
element.sim.elscore = function(clustering1, clustering2, alpha=0.9, r=1.,
                               r2=NULL, rescale_path_type='max',
                               relabeled_elements=NULL,
                               ppr_implementation='prpack'){
  # Error handling for comparisons
  if (length(clustering1) != length(clustering2)){
    stop('clustering1 and clustering2 do not have the same length')
  } else if (any(names(clustering1) != names(clustering2))){
    stop('Not all elements of clustering1 and clustering2 are the same')
  }
  # convert to character strings
  clustering1 = as.character(clustering1)
  clustering2 = as.character(clustering2)

  # the rows and columns of the affinity matrix correspond to relabeled elements
  #if (is.null(relabeled_elements)){
  #  relabeled_elements = relabel_objects(clustering1)
  #}

  if (is.null(r2)){
    r2 = r
  }

  # make the two affinity matrices
  clu_affinity_matrix1 = make.affinity.matrix(clustering1, alpha=alpha, r=r,
                                       rescale_path_type=rescale_path_type,
                                       relabeled_elements=relabeled_elements,
                                       ppr_implementation=ppr_implementation)
  clu_affinity_matrix2 = make.affinity.matrix(clustering2, alpha=alpha, r=r2,
                                       rescale_path_type=rescale_path_type,
                                       relabeled_elements=relabeled_elements,
                                       ppr_implementation=ppr_implementation)

  # use the corrected L1 similarity
  nodeScores = corrected.L1(clu_affinity_matrix1, clu_affinity_matrix2,
                            alpha=alpha)

  return(list(scores=nodeScores, elements=relabeled_elements))
}

clu2elm_dict = function(clustering){
  clu2elm_dict = list()
  for (clust in unique(clustering)){
    clu2elm_dict[[clust]] = which(clustering==clust)
  }
  return(clu2elm_dict)
}

# needs to be fixed for overlapping clustering functionality probably
elm2clu_dict = function(clustering){
  clu2elm_dict = list()
  for (clust in unique(clustering)){
    clustering.clu2elm_dict[[clust]] = which(clustering==clust)
  }
  return(clu2elm_dict)
}

corrected.L1 = function(x, y, alpha){
  res = 1 - 1/(2 * alpha) * Matrix::rowSums(abs(x - y))
  return(res)
}

make.affinity.matrix = function(clustering, alpha=0.9, r=1.,
                                rescale_path_type='max',
                                relabeled_elements=NULL,
                                ppr_implementation='prpack'){
  # check if the clustering is a partition
  clustering.is_disjoint = TRUE
  clustering.is_hierarchical = FALSE
  if (clustering.is_disjoint & !clustering.is_hierarchical){
    pprscore = ppr.partition(clustering=clustering, alpha=alpha,
                             relabeled_elements=relabeled_elements)
  } else {
    # otherwise we have to create the cielg and numberically solve for the
    # personalize page-rank distribution
    cielg = make_cielg(clustering=clustering, r=r,
                       rescale_path_type=rescale_path_type,
                       relabeled_elements=relabeled_elements)
    pprscore = numerical_ppr_scores(cielg, clustering, alpha=alpha,
                                    relabeled_elements=relabeled_elements,
                                    ppr_implementation=ppr_implementation)
  }
  return(pprscore)
}

# should be finished
ppr.partition = function(clustering, alpha=0.9, relabeled_elements=NULL){
  clustering.clu2elm_dict = clu2elm_dict(clustering)

  n.elements = length(clustering)
  ppr = matrix(0, nrow=n.elements, ncol=n.elements)

  for (i in 1:length(clustering.clu2elm_dict)){
    clustername=names(clustering.clu2elm_dict)[i]
    clusterlist = clustering.clu2elm_dict[[i]]
    Csize = length(clusterlist)
    ppr_result = alpha / Csize * matrix(1, nrow=Csize, ncol=Csize) +
      diag(Csize) * (1.0 - alpha)
    ppr[clusterlist, clusterlist] = ppr_result
  }

  return(ppr)
}

# Create the cluster-induced element graph for a Clustering
make.cielg = function(clustering, r=1.0, rescale_path_type='max',
                      relabeled_elements=NULL){
  # the rows and columns of the affinity matrix correspond to relabeled
  # elements
  if (is.null(relabeled_elements)){
    relabeled_elements = relabel_objects(clustering1)
  }

  # TODO: fix rest of this function
  # the hierarchical weight function
  if (clustering.is_hierarchical){
    cluster_height = clustering.hier_graph.rescale(
      rescale_path_type=rescale_path_type)

    weight_function = function(c){
      return(exp(r * (cluster_height.get(c, 0.0))))
    }
    clu2elm_dict = clustering.hier_clusdict()
  } else {
      weight_function = function(c){
        return(1)
      }
      clu2elm_dict = clustering.clu2elm_dict
  }

  relabeled_clusters = relabel_objects(clustering.clusters)

  edge_seq = list()
  edge_weight_seq = list()
  for (i in 1:length(clu2elm_dict.items)){
    c = names(clu2elm_dict.items)[i]
    element_list = clu2elm_dict.items[[i]]

    cstrength = weight_function(c)
    for (el in element_list){
      edge_seq.append(c(relabeled_elements[el], relabeled_clusters[c]))
      edge_weight_seq.append(cstrength)
    }
  }

  edge_seq = unlist(edge_seq)

  # use sparseMatrix
  bipartite_adj = spsparse.coo_matrix(c(edge_weight_seq,
                                       edge_seq[, 1], edge_seq[, 2]),
                                      dims=c(clustering.n_elements,
                                             clustering.n_clusters))


  proj1 = spsparse.coo_matrix(bipartite_adj / bipartite_adj.sum(axis=1))
  proj2 = spsparse.coo_matrix(bipartite_adj / bipartite_adj.sum(axis=0))
  projected_adj = proj1.dot(proj2.T)$tocoo()
  cielg = igraph.Graph(list(zip(projected_adj.row.tolist(), projected_adj.col.tolist())),
                       edge_attrs={'weight': projected_adj.data.tolist()}, directed=TRUE)
  return(cielg)
}

# TODO: finish below
find_groups_in_cluster = function(clustervs, elementgroupList){
  vertices = igraph::V(clustervs)

  if (length(intersect(vertices, x))>0){
    #?
  }
}

# TODO: finish below
numerical_ppr_scores = function(cielg, clustering, alpha=0.9,
                                relabeled_elements=NULL,
                                ppr_implementation='prpack'){
  # xx
}

# this function should be finished
get_sparse_transition_matrix = function(graph){
  transition_matrix = igraph::as_adjacency_matrix(graph, attr='weight',
                                                  sparse=TRUE)
  transition_matrix = transition_matrix / Matrix::rowSums(transition_matrix)
  return(transition_matrix)
}

# should be finished, tests look OK
calculate_ppr_with_power_iteration = function(W_matrix, index, alpha=0.9,
                                              repetition=1000, th=1e-4){
  total_length = nrow(W_matrix)
  e_s = Matrix::sparseMatrix(i=1, j=index, x=1, dims=c(1, total_length))
  p = Matrix::sparseMatrix(i=1, j=index, x=1, dims=c(1, total_length))
  for (i in 1:repetition){
    new_p =  ((1-alpha) * e_s) + ((alpha) * (p %*% W_matrix))
    if (max(abs(new_p - p)) < th){
      p = new_p
      break
    }
    p = new_p
  }

  return(as.vector(p))
}

# TODO: possibly change so sim_matrix is actual matrix, or dataframe
# otherwise finished I think
element_sim_matrix = function(clustering_list, alpha=0.9, r=1.,
                              rescale_path_type='max'){
  # relabeled_elements = relabel_objects(clustering_list[0].elements)

  affinity_matrix_list = lapply(clustering_list, function(x)
    make_affinity_matrix(clustering, alpha=alpha, r=r,
    rescale_path_type=rescale_path_type, relabeled_elements=relabeled_elements))

  Nclusterings = length(clustering_list)
  sim_matrix = rep(0, floor(Nclusterings * (Nclusterings - 1) / 2))
  icompare = 0
  for (i.clustering in 1:Nclusterings){
    for (j.clustering in (i+1):Nclusterings){
      if (j>i){
        break
      }
      sim_matrix[icompare] = mean(cL1(affinity_matrix_list[i.clustering],
                                      affinity_matrix_list[j.clustering],
                                      alpha))
      icompare = icompare + 1
    }
  }

  return(sim_matrix)
}


method.heatmaps = function(df.ecs){
  ggplot2::ggplot(data = subset(df.ecs, .data$n1==.data$n2),
                  ggplot2::aes(x=.data$method1, y=.data$method2,
                               fill=.data$ecs)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(.~n1) +
    ggplot2::labs(title='EC similarity across methods') +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#012345",
                                 limits=c(0,1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
}

n.genes.heatmaps = function(df.ecs){
  ggplot2::ggplot(data = subset(df.ecs, .data$method1==.data$method2),
                  ggplot2::aes(x=.data$n1, y=.data$n2, fill=.data$ecs)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(.~method1) +
    ggplot2::labs(title='EC similarity across n.genes') +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#012345",
                                 limits=c(0,1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1))
}
