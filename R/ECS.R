#' @importFrom foreach %dopar%
NULL

#' The Element-Centric Clustering Similarity
#'
#' @description Calculates the average element-centric similarity between two
#' Clustering objects.
#'
#' @param clustering1 The first Clustering.
#' @param clustering2 The second Clustering.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#'
#' @return The average element-wise similarity between the two Clusterings.
#' @export
#'
#' @examples
#' km.res = kmeans(mtcars, 3)$cluster
#' km.clustering = create_clustering(km.res)
#' hc.res = hclust(dist(mtcars))
#' hc.clustering = create_clustering(hc.res)
#' element_sim(km.clustering, hc.clustering)
element_sim = function(clustering1,
                       clustering2,
                       alpha = 0.9){
  element.scores = element_sim_elscore(clustering1,
                                       clustering2,
                                       alpha)
  return(mean(element.scores))
}

#' The Element-Centric Clustering Similarity for each Element
#'
#' @description Calculates the element-wise element-centric similarity between
#' two Clustering objects.
#'
#' @param clustering1 The first Clustering.
#' @param clustering2 The second Clustering.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#'
#' @return Vector of element-centric similarity between the two clusterings for
#'  each element.
#' @export
#'
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
#'
#' @examples
#' km.res = kmeans(iris[,1:4], centers=8)$cluster
#' km.clustering = create_clustering(km.res)
#' hc.res = hclust(dist(iris[,1:4]))
#' hc.clustering = create_clustering(hc.res)
#' element_sim_elscore(km.clustering, hc.clustering)
element_sim_elscore = function(clustering1, clustering2, alpha = 0.9){

  # if both clusterings are membership vectors, calculate the ecs without
  # creating a ClustAssess object
  if(any(class(clustering1) %in% c("numeric", "integer", "factor", "character")) &&
     any(class(clustering2) %in% c("numeric", "integer", "factor", "character"))) {

    if (length(clustering1) != length(clustering2)){
      stop('clustering1 and clustering2 do not have the same length.')
    }
    if (any(names(clustering1) != names(clustering2))) {
      stop('Not all elements of clustering1 and clustering2 are the same.')
    }

    node.scores = corrected_l1_mb(clustering1,
                                  clustering2)

    names(node.scores) = names(clustering1)
    return(node.scores)
  }


  if(class(clustering1) != "Clustering") {
    clustering1 = create_clustering(clustering1)
  }

  if(class(clustering2) != "Clustering") {
    clustering2 = create_clustering(clustering2)
  }

  # Make sure clusterings are comparable
  if (clustering1@n_elements != clustering2@n_elements){
    stop('clustering1 and clustering2 do not have the same length.')
  } else if (any(names(clustering1) != names(clustering2))){
    stop('Not all elements of clustering1 and clustering2 are the same.')
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
  for (i in 1:length(clustering)){
    clust = clustering[i]
    clu2elm_dict[[clust]] = c(clu2elm_dict[[clust]], i)
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

corrected_l1_mb = function(mb1, mb2, alpha = 0.9) {
  n = length(mb1)

  if(class(mb1) != "character") {
    mb1 = as.character(mb1)
  }

  if(class(mb2) != "character") {
    mb2 = as.character(mb2)
  }

  clu2elm_dict_1 = create_clu2elm_dict(mb1)
  clu2elm_dict_2 = create_clu2elm_dict(mb2)

  # the possible number of different ecs values is n x m
  # where n is the number of clusters of the first partition
  # and n the number of clusters of the second partition
  unique_ecs_vals = matrix(NA, nrow = length(clu2elm_dict_1), ncol = length(clu2elm_dict_2))
  rownames(unique_ecs_vals) = names(clu2elm_dict_1)
  colnames(unique_ecs_vals) = names(clu2elm_dict_2)

  ecs = rep(0, n)

  ppr1 = rep(0, n)
  ppr2 = rep(0, n)


  # iterate through each point of the membership vector
  for(i in 1:n) {
    # check if the similarity between of this pair of clusters was already calculated
    if(is.na(unique_ecs_vals[mb1[i], mb2[i]])) {
      clusterlist1 = clu2elm_dict_1[[mb1[i]]]
      Csize1 = length(clusterlist1)

      # get the values from the affinity matrix of the first partition
      ppr1[clusterlist1] = alpha / Csize1
      ppr1[i] = 1.0 - alpha + alpha / Csize1

      clusterlist2 = clu2elm_dict_2[[mb2[i]]]
      Csize2 = length(clusterlist2)

      # get the values from the affinity matrix of the second partition
      ppr2[clusterlist2] = alpha / Csize2
      ppr2[i] = 1.0 - alpha + alpha / Csize2


      # calculate the sum of the difference between the affinity matrices
      ecs[i] = sum(abs(ppr1-ppr2)) # could be optimized?

      ppr1[clusterlist1] = 0
      ppr2[clusterlist2] = 0

      # store the similarity between the pair of clusters
      unique_ecs_vals[mb1[i], mb2[i]] = ecs[i]
    } else {
      # if yes, just copy the value
      ecs[i] = unique_ecs_vals[mb1[i], mb2[i]]
    }
  }

  # perform the last calculations to obtain the ECS at each point
  return(1 - 1 / (2 * alpha) * ecs)
}

# PPR calculation for partition
ppr_partition = function(clustering, alpha=0.9){
  ppr = matrix(0, nrow=clustering@n_elements, ncol=clustering@n_elements)

  for (i in 1:length(clustering@clu2elm_dict)){
    clusterlist = clustering@clu2elm_dict[[i]]
    Csize = length(clusterlist)
    ppr_res = matrix(alpha/Csize, nrow=Csize, ncol=Csize)
    diag(ppr_res) = 1.0 - alpha + alpha/Csize
    ppr[clusterlist, clusterlist] = ppr_res
  }

  ppr = Matrix::Matrix(ppr,
                       sparse=TRUE)
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
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
#'
#' @examples
#' clustering.list = list()
#' for (i in 1:20){
#'   km.res = kmeans(mtcars, 3)$cluster
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

element_sim_matrix_new = function(clustering_list,
                                  output_type='matrix',
                                  alpha = 0.9,
                                  r = 1,
                                  rescale_path_type = "max",
                                  ppr_implementation = "prpack",
                                  dist_rescaled = FALSE,
                                  row_normalize = TRUE,
                                  ncores = 1) {
  if (!(output_type %in% c('data.frame', 'matrix'))){
    stop('output_type must be data.frame or matrix.')
  }

  # check if all objects are flat disjoint membership vectors
  are_all_flat_disjoint = sapply(clustering_list, function(x) {
    any(class(x) %in% c("numeric", "integer", "factor", "character"))
  })

  # if the condition is met, perform element consistency using only the membership vector
  if(all(are_all_flat_disjoint == TRUE)) {
    element_sim_matrix_flat_disjoint(mb_list = clustering_list,
                                     ncores = ncores,
                                     alpha = alpha,
                                     output_type = output_type)
  }

  # create clustassess objects

  my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
  )

  doParallel::registerDoParallel(cl = my_cluster)

  clustering_object_list = foreach::foreach(obj = clustering_list,
                                   .noexport = all_vars[!(all_vars %in% needed_vars)],
                                   .packages = "ClustAssess") %dopar% {
                                     if(class(obj) == "Clustering") {
                                       return(obj)
                                     }

                                     if(class(obj) == "hclust") {
                                       return(create_clustering(clustering_result = obj,
                                                                alpha = alpha,
                                                                r = r,
                                                                rescale_path_type = rescale_path_type,
                                                                ppr_implementation = ppr_implementation,
                                                                dist_rescaled = dist_rescaled))
                                     }

                                     if(any(class(obj) %in% c("matrix", "Matrix"))) {
                                       return(create_clustering(clustering_result = obj,
                                                                alpha = alpha,
                                                                ppr_implementation = ppr_implementation,
                                                                row_normalize = row_normalize))
                                     }

                                     create_clustering(clustering_result = obj,
                                                       alpha = alpha)
                                   }

  parallel::stopCluster(cl = my_cluster)

  n.clusterings = length(clustering_object_list)
  sim_matrix = matrix(NA, nrow=n.clusterings, ncol=n.clusterings)
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_object_list[[i]]@affinity_matrix
    for (j in (i+1):n.clusterings){
      j.aff = clustering_object_list[[j]]@affinity_matrix

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


element_sim_matrix_flat_disjoint = function(mb_list, ncores = 2, alpha = 0.9, output_type='matrix') {
  if (!(output_type %in% c('data.frame', 'matrix'))){
    stop('output_type must be data.frame or matrix.')
  }


  n_clusterings = length(mb_list)
  first_index = unlist(sapply(1:(n_clusterings-1), function(i) { rep(i, n_clusterings-i)}))
  second_index = unlist(sapply(1:(n_clusterings-1), function(i) { (i+1):n_clusterings}))
  n_combinations = n_clusterings * (n_clusterings-1) / 2
  ncores = min(ncores, n_combinations)

  my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my_cluster)

  ecs_values = foreach::foreach(i = 1:n_combinations, .export = c("corrected_l1_mb", "create_clu2elm_dict"), .noexport = c("my_cluster"), .combine = "c") %dopar% {
    mean(corrected_l1_mb(mb_list[[first_index[i]]],
                         mb_list[[second_index[i]]],
                         alpha))
  }

  parallel::stopCluster(cl = my_cluster)

  sim_matrix = matrix(NA, nrow=n_clusterings, ncol=n_clusterings)
  sim_matrix[lower.tri(sim_matrix, diag=FALSE)]= ecs_values
  sim_matrix = t(sim_matrix)
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




get_membership = function(clust_obj) {
  membership = rep(0,  clust_obj@n_elements)
  no_clusters = length(clust_obj@clu2elm_dict)

  for(i in 1:no_clusters) {
    membership[clust_obj@clu2elm_dict[[i]] ] = i
  }

  membership
}



# determine whether two flat membership vectors are identical
are_identical_memberships = function(mb1, mb2) {
  contingency_table = (table(mb1, mb2) != 0)

  if(nrow(contingency_table) != ncol(contingency_table)) {
    return(FALSE)
  }

  no_different_elements = colSums(contingency_table)

  if(any(no_different_elements != 1)) {
    return(FALSE)
  }

  no_different_elements = rowSums(contingency_table)

  return(all(no_different_elements == 1))
}



# merge partitions that are identical (meaning the ecs threshold is 1)
# the order parameter is indicating whether to sort the merged objects
# based on their frequency
merge_identical_partitions = function(clustering_list, order = T) {
  # merge the same memberships into the same object
  merged_partitions = list(clustering_list[[1]])
  no_partitions = length(clustering_list)

  if(no_partitions == 1) {
    return(clustering_list)
  }

  for(i in 2:no_partitions) {
    no_partitions = length(merged_partitions)

    assigned_partition = -1

    for(j in 1:no_partitions) {
      are_identical = are_identical_memberships(merged_partitions[[j]]$mb, clustering_list[[i]]$mb)
      if(are_identical) {
        assigned_partition = j
        break
      }
    }

    if(assigned_partition != -1) {
      merged_partitions[[assigned_partition]]$freq = merged_partitions[[assigned_partition]]$freq + clustering_list[[i]]$freq
    } else {
      merged_partitions[[no_partitions+1]] = clustering_list[[i]]
    }
  }

  if(order) {
    ordered_indices = order(sapply(merged_partitions, function(x) { x$freq }), decreasing = T)
    merged_partitions = order_list(merged_partitions, ordered_indices)
  }


  merged_partitions
}


# merge the partitions when the ecs threshold is not 1
merge_partitions_ecs = function(partition_list, ecs_thresh = 0.99, ncores = 2, order = T) {

  partition_groups = list()
  nparts = length(partition_list)

  if(nparts == 1) {
    return(partition_list)
  }

  sim_matrix = element_sim_matrix_flat_disjoint(lapply(partition_list, function(x) { x$mb }), ncores = ncores) # use `create_clustering` for the original version


  for(i in 1:nparts) {
    sim_matrix[i,i] = NA
    partition_groups[[as.character(i)]] = i
  }

  for(no_changes in 1:(nparts-1)) {
    if(max(sim_matrix, na.rm = T) < ecs_thresh) {
      break
    }



    index = which.max(sim_matrix)[1]
    first_cluster = index %% nparts

    second_cluster = index %/% nparts + 1
    if(first_cluster == 0) {
      first_cluster = nparts
      second_cluster = second_cluster - 1
    }


    partition_groups[[as.character(second_cluster)]] = NULL
    partition_groups[[as.character(first_cluster)]] = c(partition_groups[[as.character(first_cluster)]], second_cluster)

    for(i in 1:nparts) {
      if(first_cluster < i) {
        if(!is.na(sim_matrix[first_cluster, i])) {
          sim_matrix[first_cluster, i] = min(sim_matrix[first_cluster, i], sim_matrix[second_cluster, i], sim_matrix[i, second_cluster], na.rm = T)
        }
      } else {
        if(!is.na(sim_matrix[i, first_cluster])) {
          sim_matrix[i, first_cluster] = min(sim_matrix[i, first_cluster], sim_matrix[second_cluster, i], sim_matrix[i, second_cluster], na.rm = T)
        }
      }

    }

    sim_matrix[second_cluster, ] = NA
    sim_matrix[, second_cluster] = NA
  }

  merged_index = 1
  merged_partitions = list()

  for(kept_partition in names(partition_groups)) {
    merged_partitions[[merged_index]] = partition_list[[as.numeric(kept_partition)]]
    merged_partitions[[merged_index]]$freq = sum(sapply(partition_groups[[kept_partition]], function(x) { partition_list[[as.numeric(x)]]$freq }))

    merged_index = merged_index +1
  }

  if(order) {
    ordered.indices = order(sapply(merged_partitions, function(x) { x$freq }), decreasing = T)
    merged_partitions = order_list(merged_partitions, ordered.indices)
  }

  merged_partitions
}

#' Merge partitions
#' @description Merge partitions whose ECS score is above a given threshold.
#' The merging is done using a complete linkeage approach.
#'
#' @param partition_list A list of flat Clustering objects
#' @param ecs_thresh the ecs threshold object
#'
#' @return a list of the merged partitions
#' @export
#' @examples
#' initial_list = list(c(1,1,2), c(2,2,2))
#' merge_partitions(initial_list, 0.99)
#'

merge_partitions = function(partition_list, ecs_thresh = 0.99, ncores = 2) {
  # check the type of object that is provided in the list
  if(class(partition_list[[1]]) != "list") {
    partition_list = lapply(partition_list, function(x) {
      list(mb = x,
           freq = 1)
    })
  } else {
    if(!all(c("mb", "freq") %in% names(partition_list[[1]]))) {
      return(lapply(partition_list, function(sublist) {
        merge_partitions(sublist, ecs_thresh)
      }))
    }
  }

  if(ecs_thresh == 1) {
    return(merge_identical_partitions(partition_list))
  }

  merge_partitions_ecs(partition_list, ecs_thresh, ncores = ncores)
}


#' Element-Wise Consistency Between a Set of Clusterings
#' @description Inspect the consistency of a set of clusterings by calculating
#' their element-wise clustering consistency (also known as element-wise frustration).
#'
#' @param clustering_list A list of Clustering objects used to calculate
#' the element-wise consistency.
#'
#' @return a vector containing the element-wise consistency.
#' @export
#'
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
#'
#' @examples
#' clustering.list = list()
#' for (i in 1:20){
#'   km.res = kmeans(mtcars, 3)$cluster
#'   clustering.list[[i]] = create_clustering(km.res)
#' }
#' element_consistency(clustering.list)
element_consistency = function(clustering_list){
  # make sure all clusterings have same alpha
  alphas = sapply(clustering_list, function(x) x@alpha)
  if (length(unique(alphas)) != 1){
    stop('all clusterings in clustering_list must be created with same alpha')
  }

  n.clusterings = length(clustering_list)
  consistency = rep(0, length(clustering_list[[1]]))
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_list[[i]]@affinity_matrix
    for (j in (i+1):n.clusterings){
      j.aff = clustering_list[[j]]@affinity_matrix

      consistency = consistency + corrected_L1(i.aff, j.aff, alphas[1])
    }
  }
  consistency = consistency / (n.clusterings * (n.clusterings-1) / 2)
  return(consistency)
}

element_consistency_new = function(clustering_list,
                                   alpha = 0.9,
                                   r = 1,
                                   rescale_path_type = "max",
                                   ppr_implementation = "prpack",
                                   dist_rescaled = FALSE,
                                   row_normalize = TRUE,
                                   ecs_thresh = NULL,
                                   ncores = 1) {

  # check if all objects are flat disjoint membership vectors
  are_all_flat_disjoint = sapply(clustering_list, function(x) {
    any(class(x) %in% c("numeric", "integer", "factor", "character"))
  })

  # if the condition is met, perform element consistency using only the membership vector
  if(all(are_all_flat_disjoint == TRUE)) {
    if(is.null(ecs_thresh)) {
      return(weighted_element_consistency_new(clustering_list = clustering_list,
                                              ncores = ncores))
    } else {
      final_clustering_list = merge_partitions(partition_list = clustering_list,
                                               ecs_thresh = ecs_thresh,
                                               ncores = ncores)
      return(weighted_element_consistency_new(clustering_list = lapply(final_clustering_list, function(x) { x$mb }),
                                              weights = sapply(final_clustering_list, function(x) { x$freq }),
                                              ncores = ncores))
    }

  }

  # create clustassess objects

  my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
  )

  doParallel::registerDoParallel(cl = my_cluster)


  # all_vars = ls()
  # needed_vars = c("r", "rescale_path_type", "ppr_implementation", "dist_rescaled", "alpha", "row_normalize")
  # required_methods = c("create_clustering", "create_flat_overlapping_clustering",
  #                      "create_flat_disjoint_clustering", "create_clu2elm_dict",
  #                      "ppr_partition", "create_elm2clu_dict_overlapping",
  #                      "make_cielg_overlapping", "numerical_ppr_scores")

  clustering_object_list = foreach::foreach(obj = clustering_list,
                                   .noexport = all_vars[!(all_vars %in% needed_vars)],
                                   .packages = "ClustAssess") %dopar% {
                                     # if the object is already a ClustAssess one, we will just return it
                                     if(class(obj) == "Clustering") {
                                       return(obj)
                                     }

                                     if(class(obj) == "hclust") {
                                       return(create_clustering(clustering_result = obj,
                                                                alpha = alpha,
                                                                r = r,
                                                                rescale_path_type = rescale_path_type,
                                                                ppr_implementation = ppr_implementation,
                                                                dist_rescaled = dist_rescaled))
                                     }

                                     if(any(class(obj) %in% c("matrix", "Matrix"))) {
                                       return(create_clustering(clustering_result = obj,
                                                                alpha = alpha,
                                                                ppr_implementation = ppr_implementation,
                                                                row_normalize = row_normalize))
                                     }

                                     create_clustering(clustering_result = obj,
                                                       alpha = alpha)
                                   }

  parallel::stopCluster(cl = my_cluster)

  # calculate the consistency between the ClustAssess objects

  n.clusterings = length(clustering_object_list)
  consistency = rep(0, length(clustering_object_list[[1]]))
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_object_list[[i]]@affinity_matrix
    for (j in (i+1):n.clusterings){
      j.aff = clustering_object_list[[j]]@affinity_matrix

      consistency = consistency + corrected_L1(i.aff, j.aff, alpha)
    }
  }
  consistency = consistency / (n.clusterings * (n.clusterings-1) / 2)
  return(consistency)
}

#' Element-Wise Consistency Between a Weighted Set of Disjoint Clusterings
#' @description Inspect the consistency of a set of clusterings by calculating
#' their element-wise clustering consistency (also known as element-wise frustration).
#'
#' @param clustering_list A list of Clustering objects used to calculate
#' the element-wise consistency.
#'
#' @return a vector containing the element-wise consistency.
#' @export
#'
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y

element_consistency_disjoint = function(clustering_list) {

  n.clusterings = length(clustering_list)

  # if(class(clustering_list[[1]]) != "list" && all(c("mb", "freq") %in% names(clustering_list[[1]]))) {
  #   clustering_list = lapply(1:n.clusterings, function(index) {
  #     list(mb = clustering_list[[index]],
  #          freq = 1)
  #   })
  # }

  # merge the same memberships into the same object
  final_clustering_list = merge_identical_partitions(clustering_list)


  # order the objects decreasing based on their frequency
  # ordered_indices = order(sapply(final_clustering_list, function(x) { x$freq }), decreasing = T)
  # final_clustering_list = order_list(final_clustering_list, ordered_indices)


  weighted_element_consistency(lapply(final_clustering_list, function(x) { create_clustering(x$mb) }),
                               sapply(final_clustering_list, function(x) { x$freq }))

}




#' Element-Wise Consistency Between a Weighted Set of Clusterings
#' @description Inspect the consistency of a set of clusterings by calculating
#' their element-wise clustering consistency (also known as element-wise frustration).
#'
#' @param clustering_list A list of Clustering objects used to calculate
#' the element-wise consistency.
#' @param weights A list of weights representing the frequency of each clustering
#' object
#'
#' @return a vector containing the element-wise consistency.
#' @export
#'
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
#'
#' @examples
#' clustering.list = list()
#' for (i in 1:20){
#'   km.res = kmeans(mtcars, 3)$cluster
#'   clustering.list[[i]] = create_clustering(km.res)
#' }
#' weighted_element_consistency(clustering.list, rep(1,20))
#'

weighted_element_consistency = function(clustering_list, weights = NULL){
  # make sure all clusterings have same alpha
  alphas = sapply(clustering_list, function(x) x@alpha)
  if (length(unique(alphas)) != 1){
    stop('all clusterings in clustering_list must be created with same alpha')
  }

  n.clusterings = length(clustering_list)

  if(n.clusterings == 1) {
    return(rep(1, length(clustering_list[[1]])))
  }

  consistency = rep(0, length(clustering_list[[1]]))

  if(is.null(weights)) {
    weights = rep(1, n.clusterings)
  }

  no_identical_comparisons = 0
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_list[[i]]@affinity_matrix
    for (j in (i+1):n.clusterings){
      j.aff = clustering_list[[j]]@affinity_matrix

      # compute the consistency between the two different partitions i and j
      # the consistency is multiplied by the weights (or frequencies) of the two partitions
      # in order to cover all pairwise combinations
      consistency = consistency + corrected_L1(i.aff, j.aff, alphas[1]) * weights[i] * weights[j]
    }

    # if a partition has a weight greater than 1, it means, in the unweighted case,
    # having to calculate the ECS between `weight` identical partitions
    # weights[i] * (weights[i]-1) / 2 denotes the number of all possible combinations
    # of pairing those identical memberships
    no_identical_comparisons = no_identical_comparisons + weights[i] * (weights[i]-1) / 2
  }

  no_identical_comparisons = no_identical_comparisons + weights[n.clusterings] * (weights[n.clusterings]-1) / 2

  # at the end, add the results of comparing a partition to itself (resulting in a vector
  # of ones)
  consistency = consistency + rep(1, length(consistency)) * no_identical_comparisons


  # divide over the number of all possible combinations to normalize the results
  # between 0 and 1
  consistency = consistency / (sum(weights) * (sum(weights)-1) / 2)
  return(consistency)
}

weighted_element_consistency_new = function(clustering_list,
                                            weights = NULL,
                                            ncores = 2) {
  n_clusterings = length(clustering_list)

  if(n_clusterings == 1) {
    return(rep(1, length(clustering_list[[1]])))
  }



  if(is.null(weights)) {
    weights = rep(1, n_clusterings)
  }

  first_index = unlist(sapply(1:(n_clusterings-1), function(i) { rep(i, n_clusterings-i)}))
  second_index = unlist(sapply(1:(n_clusterings-1), function(i) { (i+1):n_clusterings}))
  n_combinations = n_clusterings * (n_clusterings-1) / 2
  needed_vars = c("first_index", "second_index", "clustering_list", "weights")

  ncores = min(ncores, n_combinations)

  my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
  )

  doParallel::registerDoParallel(cl = my_cluster)


  all_vars = ls()

  consistency = foreach::foreach(i = 1:n_combinations, .noexport = all_vars[!(all_vars %in% needed_vars)], .export = c("corrected_l1_mb", "create_clu2elm_dict"), .combine = "+") %dopar% {
    corrected_l1_mb(clustering_list[[first_index[i]]],
                    clustering_list[[second_index[i]]]) * weights[first_index[i]] * weights[second_index[i]]
  }


  parallel::stopCluster(cl = my_cluster)

  no_identical_comparisons = sum(sapply(weights, function(w) { w*(w-1) / 2}))

  consistency = consistency + rep(1, length(consistency)) * no_identical_comparisons

  consistency = consistency / (sum(weights) * (sum(weights)-1) / 2)

  return(consistency)
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
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
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
  avg_agreement = rep(0, length(clustering_list[[1]]))
  ref_aff = reference_clustering@affinity_matrix
  for (i in 1:(n.clusterings-1)){
    i.aff = clustering_list[[i]]@affinity_matrix
    avg_agreement = avg_agreement + corrected_L1(ref_aff, i.aff, alphas[1])
  }
  avg_agreement = avg_agreement / n.clusterings
  return(avg_agreement)
}


#' @importClassesFrom Matrix Matrix
setOldClass("Matrix::Matrix")
#' The Clustering Class
#' @description A class containing relevant data for comparing clusterings,
#' including the affinity matrix for the Clustering.
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
#' @export
#'
#' @examples
#' km.res = kmeans(mtcars, 3)$cluster
#' km.clustering = create_clustering(km.res)
#' hc.res = hclust(dist(mtcars))
#' hc.clustering = create_clustering(hc.res)
#' element_sim(km.clustering, hc.clustering)
Clustering <- setClass("Clustering",
                       representation(names="character",
                                      n_elements="numeric",
                                      is_hierarchical="logical",
                                      is_disjoint="logical",
                                      alpha="numeric",
                                      r="numeric",
                                      elm2clu_dict="list",
                                      clu2elm_dict="list",
                                      affinity_matrix="Matrix"))

#' Create Clustering Object
#' @description Creates a Clustering object from the output of a clustering
#' method.
#'
#' @param clustering_result The clustering result, either:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#' @param ppr_implementation Choose a implementation for personalized
#' page-rank calcuation:
#' * 'prpack': use PPR alogrithms in igraph.
#' * 'power_iteration': use power_iteration method.
#' @param row_normalize Whether to normalize all rows in clustering_result
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r A numeric hierarchical scaling parameter.
#' @param rescale_path_type A string; rescale the hierarchical height by:
#' * 'max' : the maximum path from the root.
#' * 'min' : the minimum path form the root.
#' * 'linkage' : use the linkage distances in the clustering.
#' @param dist_rescaled A logical: if TRUE, the linkage distances are linearly
#' rescaled to be in-between 0 and 1.
#' @param ... This argument is not used.
#'
#' @return A Clustering object.
#' @export
#'
#' @md
#' @examples
#' km.res = kmeans(mtcars, 3)$cluster
#' km.clustering = create_clustering(km.res)
#' hc.res = hclust(dist(mtcars))
#' hc.clustering = create_clustering(hc.res)
#' element_sim(km.clustering, hc.clustering)
setGeneric("create_clustering",
           function(clustering_result, ...)
             standardGeneric("create_clustering") )

#' @describeIn create_clustering Create Clustering Object from Numeric Vector
setMethod("create_clustering",
          signature(clustering_result="numeric"),
          function(clustering_result,
                   alpha=0.9) {
            # convert to character
            element.names = names(clustering_result)
            clustering_result = as.character(clustering_result)
            names(clustering_result) = element.names

            # create Clustering object in separate function
            return(create_flat_disjoint_clustering(clustering_result,
                                                   alpha))
          })

#' @describeIn create_clustering Create Clustering Object from Character Vector
setMethod("create_clustering",
          signature(clustering_result="character"),
          function(clustering_result,
                   alpha=0.9) {
            # create Clustering object in separate function
            return(create_flat_disjoint_clustering(clustering_result,
                                                   alpha))
          })

#' @describeIn create_clustering Create Clustering Object from Factor Vector
setMethod("create_clustering",
          signature(clustering_result="factor"),
          function(clustering_result,
                   alpha=0.9) {
            # convert to character
            element.names = names(clustering_result)
            clustering_result = as.character(clustering_result)
            names(clustering_result) = element.names

            # create Clustering object in separate function
            return(create_flat_disjoint_clustering(clustering_result,
                                                   alpha))
          })

#' @importFrom methods new
create_flat_disjoint_clustering = function(clustering_result,
                                           alpha){
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
  clustering@affinity_matrix = ppr_partition(clustering,
                                             alpha)
  return(clustering)
}

#' @describeIn create_clustering Create Clustering Object from base matrix
setMethod("create_clustering",
          signature(clustering_result="matrix"),
          function(clustering_result,
                   alpha=0.9,
                   ppr_implementation='prpack',
                   row_normalize=TRUE) {
            clustering_result = Matrix::Matrix(clustering_result,
                                               sparse=TRUE)
            return(create_flat_overlapping_clustering(clustering_result,
                                                      alpha,
                                                      ppr_implementation,
                                                      row_normalize))
          })

#' @describeIn create_clustering Create Clustering Object from Matrix::Matrix
setMethod("create_clustering",
          signature(clustering_result="Matrix"),
          function(clustering_result,
                   alpha=0.9,
                   ppr_implementation='prpack',
                   row_normalize=TRUE) {
            # convert clustering_result to sparse if it is not already
            if (!methods::is(clustering_result,
                             'sparseMatrix')){
              clustering_result = Matrix::Matrix(clustering_result,
                                                 sparse=TRUE)
            }
            return(create_flat_overlapping_clustering(clustering_result,
                                                      alpha,
                                                      ppr_implementation,
                                                      row_normalize))
          })

#' @importFrom methods new
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

setOldClass("hclust")
#' @describeIn create_clustering Create Clustering Object from hclust
#' @importFrom methods new
setMethod("create_clustering",
          signature(clustering_result="hclust"),
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
#' Length of an Object
#' @description Get the number of elements in the Clustering.
#'
#' @param x The Clustering object.
#'
#' @return The number of elements.
#' @export
#'
#' @examples
#' km.res = kmeans(mtcars, 3)$cluster
#' km.clustering = create_clustering(km.res)
#' length(km.clustering)
setMethod("length",
          signature(x="Clustering"),
          function(x) {
            return(x@n_elements)
          })

setMethod("show",
          signature(object="Clustering"),
          function(object) {
            hierarchical = if (object@is_hierarchical) 'hierarchical' else 'flat'
            overlapping = if (object@is_disjoint) 'disjoint' else 'overlapping'
            cat(paste0('A ', hierarchical, ', ', overlapping, ' clustering of ',
                       object@n_elements, ' elements.\n'))
          })

setGeneric("print")
#' Print an Object
#' @description Prints out information about the Clustering, including
#' number of elements.
#'
#' @param x The Clustering object.
#'
#' @return The printed character string.
#' @export
#'
#' @examples
#' km.res = kmeans(mtcars, 3)$cluster
#' km.clustering = create_clustering(km.res)
#' print(km.clustering)
setMethod("print",
          signature(x="Clustering"),
          function(x) {
            hierarchical = if (x@is_hierarchical) 'hierarchical' else 'flat'
            overlapping = if (x@is_disjoint) 'disjoint' else 'overlapping'
            print(paste0('A ', hierarchical, ', ', overlapping, ' clustering of ',
                         x@n_elements, ' elements.'))
          })

