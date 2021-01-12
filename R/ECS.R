#' The Element-Centric Clustering Similarity
#'
#' @description Calculates the average element-centric similarity between two
#' Clustering objects.
#'
#' @param clustering1 The first Clustering.
#' @param clustering2 The second Clustering.
#'
#' @return The average element-wise similarity between the two Clusterings.
#' @export
#'
#' @examples
#' km.res = kmeans(iris[,1:4], 3)$cluster
#' km.clustering = create_clustering(km.res)
#' hc.res = hclust(dist(iris[,1:4]))
#' hc.clustering = create_clustering(hc.res)
#' element_sim(km.clustering, hc.clustering)
element_sim = function(clustering1,
                       clustering2){
  element.scores = element_sim_elscore(clustering1,
                                       clustering2)
  return(mean(element.scores))
}

#' The Element-Centric Clustering Similarity for each Element
#'
#' @description Calculates the element-wise element-centric similarity between
#' two Clustering objects.
#'
#' @param clustering1 The first Clustering.
#' @param clustering2 The second Clustering.
#'
#' @return Vector of element-centric similarity between the two clusterings for
#'  each element.
#' @export
#'
#' @examples
#' km.res = kmeans(iris[,1:4], 3)$cluster
#' km.clustering = create_clustering(km.res)
#' hc.res = hclust(dist(iris[,1:4]))
#' hc.clustering = create_clustering(hc.res)
#' element_sim_elscore(km.clustering, hc.clustering)
element_sim_elscore = function(clustering1, clustering2){
  # Make sure clusterings are comparable
  if (clustering1@n_elements != clustering2@n_elements){
    stop('clustering1 and clustering2 do not have the same length.')
  } else if (any(names(clustering1) != names(clustering2))){
    stop('Not all elements of clustering1 and clustering2 are the same.')
  } else if (clustering1@alpha != clustering2@alpha){
    stop('clustering1 and clustering2 were created using different alpha.')
  }

  # use the corrected L1 similarity
  node.scores = corrected_L1(clustering1@affinity_matrix,
                             clustering2@affinity_matrix,
                             clustering1@alpha)

  names(node.scores) = names(clustering1)
  return(node.scores)
}

# create the cluster -> element dictionary (list) for flat, disjoint clustering
create_clu2elm_dict = function(clustering){
  clu2elm_dict = list()
  for (clust in unique(clustering)){
    clu2elm_dict[[clust]] = which(clustering==clust)
  }
  return(clu2elm_dict)
}

# create the cluster -> element dictionary (list) for hierarchical clustering
create_clu2elm_dict_hierarchical = function(clustering){
  n.clusters = nrow(clustering$merge)
  clu2elm_dict = list()

  # add NA to initialize list
  clu2elm_dict[[2*n.clusters+1]] = NA

  # add single member clusters for leaf nodes
  for (i in 1:(n.clusters+1)){
    clu2elm_dict[[i]] = i
  }

  # traverse hc tree and add nodes to clu2elm_dict
  for (i in 1:n.clusters){
    for (j in 1:2){
      index = clustering$merge[i,j]
      if (index<0){ # for leaf nodes, just append them
        clu2elm_dict[[i+n.clusters+1]] = c(clu2elm_dict[[i+n.clusters+1]],
                                           -index)
      } else { # for previous clusters, append all leaf nodes in cluster
        clu2elm_dict[[i+n.clusters+1]] = c(clu2elm_dict[[i+n.clusters+1]],
                                           clu2elm_dict[[index+n.clusters+1]])
      }
    }
  }
  # remove the NA that was added in the beginning
  clu2elm_dict[[2*n.clusters+1]] = clu2elm_dict[[2*n.clusters+1]][-1]

  return(clu2elm_dict)
}

# create element -> cluster dictionary (list) for flat, overlapping clustering
create_elm2clu_dict_overlapping = function(clustering){
  elm2clu_dict = list()
  for (i in 1:nrow(clustering)){
    elm2clu_dict[[i]] = which(clustering[i,]>0)
  }
  return(elm2clu_dict)
}

# create element -> cluster dictionary (list) for hierarchical clustering: only
# maps elements to leaf node singleton clusters
create_elm2clu_dict_hierarchical = function(clustering){
  elm2clu_dict = list()
  for (i in 1:length(clustering$order)){
    elm2clu_dict[[i]] = i
  }
  return(elm2clu_dict)
}

# corrected L1 distance
corrected_L1 = function(x, y, alpha){
  res = 1 - 1/(2 * alpha) * Matrix::rowSums(abs(x - y))
  return(res)
}

# PPR calculation for partition
ppr_partition = function(clustering, alpha=0.9){
  ppr = Matrix::Matrix(0,
                       nrow=clustering@n_elements,
                       ncol=clustering@n_elements,
                       sparse=TRUE)

  for (i in 1:length(clustering@clu2elm_dict)){
    clustername=names(clustering@clu2elm_dict)[i]
    clusterlist = clustering@clu2elm_dict[[i]]
    Csize = length(clusterlist)
    ppr_result = alpha / Matrix::Matrix(Csize, nrow=Csize, ncol=Csize) +
      Matrix::Diagonal(Csize) * (1.0 - alpha)
    ppr[clusterlist, clusterlist] = ppr_result
  }

  return(ppr)
}

# Create the cluster-induced element graph for overlapping clustering
make_cielg_overlapping = function(bipartite_adj){
  proj1 = Matrix::Diagonal(x = 1 / Matrix::rowSums(bipartite_adj)) %*%
    bipartite_adj
  proj2 = bipartite_adj %*%
    Matrix::Diagonal(x = 1 / Matrix::colSums(bipartite_adj))
  projected_adj = proj1 %*% Matrix::t(proj2)

  cielg = igraph::graph_from_adjacency_matrix(projected_adj,
                                              weighted=TRUE,
                                              mode='directed')
  return(cielg)
}

# create linkage distances vector
get_linkage_dist = function(hierarchical_clustering, dist_rescaled){
  n.elements = length(hierarchical_clustering$order)

  if (dist_rescaled){
    maxdist = max(hierarchical_clustering$height)
    # add 1s at beggining to account for leaf node singleton clusters
    linkage_dist = c(rep(1, n.elements),
                     1.0 - hierarchical_clustering$height/maxdist)
  } else {
    # add 0s at beginning to account for leaf node singleton clusters
    linkage_dist = c(rep(0, n.elements), hierarchical_clustering$height)
  }
  return(linkage_dist)
}

# helper function to calculate path distances in graph
dist_path = function(hierarchical_clustering,
                     rescale_path_type){
  n.elements = length(hierarchical_clustering$order)
  # dist to node from final cluster containing all elements
  path_to_node = rep(0, 2*n.elements-1)
  for (i in rev(1:(n.elements-1))){
    for (j in 1:2){
      index = hierarchical_clustering$merge[i,j]
      if (index < 0){
        index = -index
      } else {
        index = index + n.elements
      }
      if (rescale_path_type=='min'){
        if (path_to_node[index] == 0){
          path_to_node[index] = path_to_node[i + n.elements] + 1
        } else {
          path_to_node[index] = min(path_to_node[index],
                                    path_to_node[i + n.elements] + 1)
        }
      } else {
        path_to_node[index] = max(path_to_node[index],
                                  path_to_node[i + n.elements] + 1)
      }

    }
  }

  # dist to node from leaves
  path_from_node = rep(0, 2*n.elements-1)
  for (i in 1:(n.elements-1)){
    index = hierarchical_clustering$merge[i,]
    for (j in 1:2){
      if (index[j] < 0){
        index[j] = -index[j]
      } else {
        index[j] = index[j] + n.elements
      }
    }
    if (rescale_path_type=='min'){
      path_from_node[i + n.elements] = min(path_from_node[index[1]],
                                           path_from_node[index[2]]) + 1

    } else {
      path_from_node[i + n.elements] = max(path_from_node[index[1]],
                                           path_from_node[index[2]]) + 1
    }
  }

  return(list(to=path_to_node, from=path_from_node))
}

# rescale the hierarchical clustering path
rescale_path = function(hierarchical_clustering,
                        rescale_path_type){
  # get path distances
  path.dist = dist_path(hierarchical_clustering, rescale_path_type)
  path_from_node = path.dist$from
  path_to_node = path.dist$to

  n.elements = length(hierarchical_clustering$order)
  rescaled_level = rep(0, n.elements)
  for (i in 1:(2*n.elements-1)){
    total_path_length = path_from_node[i] + path_to_node[i]
    if (total_path_length != 0){
      rescaled_level[i] = path_to_node[i] / total_path_length
    }
  }
  return(rescaled_level)
}

# Create the cluster-induced element graph for hierarchical clustering
make_cielg_hierarchical = function(hierarchical_clustering,
                                   clu2elm_dict,
                                   r,
                                   rescale_path_type,
                                   dist_rescaled=FALSE){
  n.elements = length(hierarchical_clustering$order)

  # rescale paths
  if (rescale_path_type == 'linkage'){
    cluster_height = get_linkage_dist(hierarchical_clustering,
                                      dist_rescaled)
  } else if (rescale_path_type %in% c('min', 'max')) {
    cluster_height = rescale_path(hierarchical_clustering,
                                  rescale_path_type)
  } else {
    stop(paste0("rescale_path_type must be one of linkage, min or max."))
  }

  # weight function for different heights
  weight_function = function(clust) return(exp(r * (cluster_height[clust])))

  edge_i=c()
  edge_j=c()
  edge_weights=c()
  for (i in 1:length(clu2elm_dict)){
    element_list = clu2elm_dict[[i]]

    cstrength = weight_function(i)
    for (el in element_list){
      edge_i = c(edge_i, el)
      edge_j = c(edge_j, i)
      edge_weights = c(edge_weights, cstrength)
    }
  }

  bipartite_adj = Matrix::sparseMatrix(i=edge_i, j=edge_j, x=edge_weights,
                                       dims=c(n.elements,
                                              2*n.elements-1))

  proj1 = Matrix::Diagonal(x = 1 / Matrix::rowSums(bipartite_adj)) %*%
    bipartite_adj
  proj2 = bipartite_adj %*%
    Matrix::Diagonal(x = 1 / Matrix::colSums(bipartite_adj))
  projected_adj = proj1 %*% Matrix::t(proj2)

  cielg = igraph::graph_from_adjacency_matrix(projected_adj,
                                              weighted=TRUE,
                                              mode='directed')
  return(cielg)
}


# utility function to find vertices with the same cluster memberships.
find_groups_in_cluster = function(clustervs, elementgroupList){
  clustervertex = unique(clustervs)

  groupings = list()
  index = 1
  for (i in 1:length(elementgroupList)){
    current = elementgroupList[[i]]
    if (length(intersect(current, clustervertex))>0){
      groupings[[index]] = current
      index = index + 1
    }
  }

  return(groupings)
}

# numerical calculation of affinity matrix using PPR
numerical_ppr_scores = function(cielg, clustering, ppr_implementation='prpack'){
  if (!(ppr_implementation %in% c('prpack', 'power_iteration'))){
    stop('ppr_implementation argument must be one of prpack or power_iteration')
  }

  # keep track of all clusters an element is a member of
  elementgroupList = list()
  for (i in 1:length(clustering@elm2clu_dict)){
    element = names(clustering)[i]
    # collapse clusters into single string
    clusters = paste(clustering@elm2clu_dict[[i]], collapse=",")
    print(clusters)

    elementgroupList[[clusters]] = c(elementgroupList[[clusters]], element)
  }

  ppr_scores = Matrix::Matrix(0,
                              nrow=igraph::vcount(cielg),
                              ncol=igraph::vcount(cielg),
                              sparse=TRUE)

  # we have to calculate the ppr for each connected component
  comps = igraph::components(cielg)
  for (i.comp in 1:comps$no){
    members = which(comps$membership == i.comp)
    clustergraph = igraph::induced_subgraph(cielg, members, impl='auto')
    cc_ppr_scores = Matrix::Matrix(0,
                                   nrow=igraph::vcount(clustergraph),
                                   ncol=igraph::vcount(clustergraph),
                                   sparse=TRUE)

    if (ppr_implementation == 'power_iteration'){
      W_matrix = get_sparse_transition_matrix(clustergraph)
    }

    elementgroups = find_groups_in_cluster(members, elementgroupList)
    for (elementgroup in elementgroups){
      # we only have to solve for the ppr distribution once per group
      vertex.index = match(elementgroup[1], members)
      # vertex = V(clustergraph)[index]
      # vertex = clustergraph.vs[cluster.index(elementgroup[0])] #????
      if (ppr_implementation == 'prpack'){
        pers = rep(0, length=length(members))
        pers[vertex.index] = 1
        ppr_res = igraph::page_rank(clustergraph,
                                    directed=TRUE,
                                    weights=NULL,
                                    damping=clustering@alpha,
                                    personalized=pers,
                                    algo='prpack')
        cc_ppr_scores[vertex.index,] = ppr_res$vector
      } else if (ppr_implementation == 'power_iteration'){
        cc_ppr_scores[vertex.index,] = calculate_ppr_with_power_iteration(
          W_matrix,
          vertex.index,
          alpha=clustering@alpha,
          repetition=1000,
          th=0.0001)
      }

      # the other vertices in the group are permutations of that solution
      elgroup.size = length(elementgroup)
      if (elgroup.size >= 2){
        for (i2 in 2:elgroup.size){
          v2 = match(elementgroup[i2], members)
          cc_ppr_scores[v2,] = cc_ppr_scores[vertex.index,]
          cc_ppr_scores[v2, vertex.index] = cc_ppr_scores[vertex.index, v2]
          cc_ppr_scores[v2, v2] = cc_ppr_scores[vertex.index, vertex.index]
        }
      }
    }

    ppr_scores[members, members] = cc_ppr_scores
  }
  return(ppr_scores)
}

# utility function to get a row-normalized sparse transition matrix
get_sparse_transition_matrix = function(graph){
  transition_matrix = igraph::as_adjacency_matrix(graph, attr='weight',
                                                  sparse=TRUE)
  transition_matrix = Matrix::Diagonal(x = 1 /
                                         Matrix::rowSums(transition_matrix)) %*%
    transition_matrix
  return(transition_matrix)
}

# power iteration calculation of PPR
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

  return(p)
}


#' Pairwise Comparison of Clusterings
#' @description Compare a set of clusterings by calculating their pairwise
#' average element-centric clustering similarities.
#'
#' @param clustering_list A list of Clustering objects to be compared with
#' element-centric similarity.
#' @param output_type A string specifying whether the output should be a
#' matrix or a data.frame.
#'
#' @return A matrix or data.frame containing the pairwise ECS values.
#' @export
#'
#' @examples
#' clustering.list = list()
#' for (i in 1:20){
#'   km.res = kmeans(iris[,1:4], 3)$cluster
#'   clustering.list[[i]] = create_clustering(km.res)
#' }
#' element_sim_matrix(clustering.list, output_type='matrix')
element_sim_matrix = function(clustering_list, output_type='matrix'){
  if (!(output_type %in% c('data.frame', 'matrix'))){
    stop('output_type must be data.frame or matrix.')
  }
  # make sure all clusterings have same alpha
  alphas = sapply(clustering_list, function(x) x@alpha)
  if (length(unique(alphas)) != 1){
    stop('all clusterings in clustering_list must be created with same alpha')
  }

  n.clusterings = length(clustering_list)
  sim_matrix = matrix(NA, nrow=n.clusterings, ncol=n.clusterings)
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_list[[i]]@affinity_matrix
    for (j in (i+1):n.clusterings){
      j.aff = clustering_list[[j]]@affinity_matrix

      sim_matrix[i, j] = mean(corrected_L1(i.aff, j.aff, alphas[1]))
    }
  }
  diag(sim_matrix) = 1

  if (output_type=='data.frame'){
    sim_matrix = data.frame(i.clustering=rep(1:n.clusterings,
                                             n.clusterings),
                            j.clustering=rep(1:n.clusterings,
                                             each=n.clusterings),
                            element_sim=c(sim_matrix))
  }

  return(sim_matrix)
}


#' Element-Wise Frustration Between a Set of Clusterings
#' @description Inspect the consistency of a set of clusterings by calculating
#' their element-wise clustering frustration.
#'
#' @param clustering_list A list of Clustering objects used to calculate
#' the element-wise frustration.
#'
#' @return a vector containing the element-wise frustration.
#' @export
#'
#' @examples
#' clustering.list = list()
#' for (i in 1:20){
#'   km.res = kmeans(iris[,1:4], 3)$cluster
#'   clustering.list[[i]] = create_clustering(km.res)
#' }
#' element_frustration(clustering.list)
element_frustration = function(clustering_list){
  # make sure all clusterings have same alpha
  alphas = sapply(clustering_list, function(x) x@alpha)
  if (length(unique(alphas)) != 1){
    stop('all clusterings in clustering_list must be created with same alpha')
  }

  n.clusterings = length(clustering_list)
  frustration = rep(0, n.clusterings)
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_list[[i]]@affinity_matrix
    for (j in (i+1):n.clusterings){
      j.aff = clustering_list[[j]]@affinity_matrix

      frustration = frustration + corrected_L1(i.aff, j.aff, alphas[1])
    }
  }
  frustration = frustration / (n.clusterings * (n.clusterings-1) / 2)
  return(frustration)
}

#' Element-Wise Average Agreement Between a Set of Clusterings
#' @description Inspect how consistently of a set of clusterings agree with
#' a reference clustering by calculating their element-wise average agreement.
#'
#' @param reference_clustering A Clustering objects for the reference clustering
#' that each clustering in clustering_list is compared to.
#' @param clustering_list A list of Clustering objects used to calculate
#' the element-wise average agreement
#'
#' @return A vector containing the element-wise average agreement.
#' @export
#'
#' @examples
#' reference.clustering = create_clustering(iris$Species)
#' clustering.list = list()
#' for (i in 1:20){
#'   km.res = kmeans(iris[,1:4], 3)$cluster
#'   clustering.list[[i]] = create_clustering(km.res)
#' }
#' element_agreement(reference.clustering, clustering.list)
element_agreement = function(reference_clustering, clustering_list){
  # make sure all clusterings have same alpha
  alphas = sapply(clustering_list, function(x) x@alpha)
  if ((reference_clustering@alpha != alphas[1]) |
      (length(unique(alphas)) != 1)){
    stop('reference_clustering and all clusterings in clustering_list must be
         created with same alpha')
  }

  n.clusterings = length(clustering_list)
  avg_agreement = rep(0, n.clusterings)
  ref.aff = reference_clustering@affinity_matrix
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_list[[i]]@affinity_matrix
    avg_agreement = avg_agreement + corrected_L1(ref.aff, i.aff, alphas[1])
  }
  avg_agreement = avg_agreement / n.clusterings
  return(avg_agreement)
}


#' @importClassesFrom Matrix Matrix
setOldClass("Matrix::Matrix")
#' The Clustering Class
#'
#' @slot names A character vector of element names; will be 1:n_elements if no
#' names were available when creating the Clustering object.
#' @slot n_elements A numeric giving the number of elements.
#' @slot is_hierarchical A logical indicating whether the clustering is
#' hierarchical or flat.
#' @slot is_disjoint A logical indicating whether the clustering is disjoint or
#' overlapping.
#' @slot alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#' @slot r A numeric hierarchical scaling parameter.
#' @slot elm2clu_dict A list giving the clusters each element is a member of.
#' @slot clu2elm_dict A list giving the element members of each cluster.
#' @slot affinity_matrix A Matrix containing the personalized pagerank
#' equilibrium distribution.
#'
#' @return
#' @export
#'
#' @examples
setClass("Clustering", representation(names="character",
                                      n_elements="numeric",
                                      is_hierarchical="logical",
                                      is_disjoint="logical",
                                      alpha="numeric",
                                      r="numeric",
                                      elm2clu_dict="list",
                                      clu2elm_dict="list",
                                      affinity_matrix="Matrix"))

setGeneric("create_clustering", function(clustering_result, ...)
  standardGeneric("create_clustering") )

#' Title
#'
#' @param clustering_result A numeric vector of cluster labels for each element.
#'
#' @return A Clustering object.
#' @export
#'
#' @examples
setMethod("create_clustering", signature(clustering_result="numeric"),
          function(clustering_result, alpha=0.9) {
  # convert to character
  clustering_result = as.character(clustering_result)

  # create Clustering object in separate function
  return(create_flat_disjoint_clustering(clustering_result, alpha))
})

setMethod("create_clustering", signature(clustering_result="character"),
          function(clustering_result, alpha=0.9) {
  # create Clustering object in separate function
  return(create_flat_disjoint_clustering(clustering_result, alpha))
})

setMethod("create_clustering", signature(clustering_result="factor"),
          function(clustering_result, alpha=0.9) {
  # convert to character
  clustering_result = as.character(clustering_result)

  # create Clustering object in separate function
  return(create_flat_disjoint_clustering(clustering_result, alpha))
})

create_flat_disjoint_clustering = function(clustering_result, alpha){
  # disjoint partitions
  n_elements = length(clustering_result)

  # names of elements
  if (is.null(names(clustering_result))){
    element.names = as.character(1:n_elements)
  } else {
    element.names = names(clustering_result)
  }

  # create dictionaries, or lists
  clu2elm_dict = create_clu2elm_dict(clustering_result)
  elm2clu_dict = list()

  # create object
  clustering = methods::new('Clustering',
                            n_elements=n_elements,
                            is_hierarchical=FALSE,
                            is_disjoint=TRUE,
                            names=element.names,
                            alpha=alpha,
                            r=0,
                            elm2clu_dict=elm2clu_dict,
                            clu2elm_dict=clu2elm_dict,
                            affinity_matrix=Matrix::Matrix())

  # calculate affinity matrix
  clustering@affinity_matrix = ppr_partition(clustering, alpha)
  return(clustering)
}


setMethod("create_clustering", signature(clustering_result="matrix"),
          function(clustering_result, alpha=0.9, ppr_implementation='prpack',
                   row_normalize=TRUE) {
  clustering_result = Matrix::Matrix(clustering_result, sparse=TRUE)
  return(create_flat_overlapping_clustering(clustering_result,
                                            alpha,
                                            ppr_implementation,
                                            row_normalize))
})

setMethod("create_clustering", signature(clustering_result="Matrix::Matrix"),
          function(clustering_result, alpha=0.9, ppr_implementation='prpack',
                   row_normalize=TRUE) {
  # convert clustering_result to sparse if it is not already
  if (!methods::is(clustering_result, 'sparseMatrix')){
    clustering_result = Matrix::Matrix(clustering_result, sparse=TRUE)
  }
  return(create_flat_overlapping_clustering(clustering_result,
                                            alpha,
                                            ppr_implementation,
                                            row_normalize))
})

create_flat_overlapping_clustering = function(clustering_result,
                                              alpha,
                                              ppr_implementation,
                                              row_normalize){
  # soft clustering
  n_elements = nrow(clustering_result)

  # check that entries are non-negative
  if (any(clustering_result<0)){
    stop('clustering_result should only contain non-negative entries.')
  }
  # check that each element belongs to at least one cluster
  if (any(Matrix::rowSums(clustering_result)==0)){
    stop('Every element of clustering_result must be in at least one cluster.')
  }

  if (row_normalize){
    # normalize clustering_result so each row sums to one
    clustering_result = clustering_result / Matrix::rowSums(clustering_result)
  }

  # names of elements
  if (is.null(rownames(clustering_result))){
    element.names = as.character(1:n_elements)
  } else {
    element.names = rownames(clustering_result)
  }

  # create dictionaries, or lists
  clu2elm_dict = list()
  elm2clu_dict = create_elm2clu_dict_overlapping(clustering_result)

  # create object
  clustering = methods::new('Clustering',
                            n_elements=n_elements,
                            is_hierarchical=FALSE,
                            is_disjoint=FALSE,
                            names=element.names,
                            alpha=alpha,
                            r=0,
                            elm2clu_dict=elm2clu_dict,
                            clu2elm_dict=clu2elm_dict,
                            affinity_matrix=Matrix::Matrix())

  # create affinity matrix
  cielg = make_cielg_overlapping(clustering_result)
  clustering@affinity_matrix = numerical_ppr_scores(cielg,
                                                    clustering,
                                                    ppr_implementation=
                                                      ppr_implementation)
  return(clustering)
}

setOldClass("stats::hclust")
setOldClass("hclust")
setMethod("create_clustering", signature(clustering_result="hclust"),
          function(clustering_result,
                   alpha=0.9,
                   r=1,
                   rescale_path_type='max',
                   ppr_implementation='prpack',
                   dist_rescaled=FALSE) {
  # hierarchical partitions
  n_elements = length(clustering_result$order)

  # names of elements
  if (is.null(clustering_result$labels)){
    element.names = as.character(1:n_elements)
  } else {
    element.names = clustering_result$labels
  }

  # create dictionaries, or lists
  clu2elm_dict = create_clu2elm_dict_hierarchical(clustering_result)
  elm2clu_dict = create_elm2clu_dict_hierarchical(clustering_result)

  # create object
  clustering = methods::new('Clustering',
                            n_elements=n_elements,
                            is_hierarchical=TRUE,
                            is_disjoint=TRUE,
                            names=element.names,
                            alpha=alpha,
                            r=r,
                            elm2clu_dict=elm2clu_dict,
                            clu2elm_dict=clu2elm_dict,
                            affinity_matrix=Matrix::Matrix())

  # create affinity matrix
  cielg = make_cielg_hierarchical(clustering_result,
                                  clu2elm_dict,
                                  r,
                                  rescale_path_type,
                                  dist_rescaled)
  clustering@affinity_matrix = numerical_ppr_scores(cielg,
                                                    clustering,
                                                    ppr_implementation=
                                                      ppr_implementation)
  return(clustering)
})

setGeneric("length")
setMethod("length", signature(x="Clustering"),
          function(x) {
  return(x@n_elements)
})

# setGeneric("show")
setMethod("show", signature(object="Clustering"),
          function(object) {
  hierarchical = if (object@is_hierarchical) 'hierarchical' else 'flat'
  overlapping = if (object@is_disjoint) 'disjoint' else 'overlapping'
  cat(paste0('A ', hierarchical, ', ', overlapping, ' clustering of ',
             object@n_elements, ' elements.\n'))
})

setGeneric("print")
setMethod("print", signature(x="Clustering"),
          function(x) {
  hierarchical = if (x@is_hierarchical) 'hierarchical' else 'flat'
  overlapping = if (x@is_disjoint) 'disjoint' else 'overlapping'
  print(paste0('A ', hierarchical, ', ', overlapping, ' clustering of ',
               x@n_elements, ' elements.'))
})

