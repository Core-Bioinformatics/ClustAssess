#' @importFrom foreach %dopar% %do%
#' @importFrom dplyr %>%
#' @importFrom methods new
NULL

#' The Element-Centric Clustering Similarity
#'
#' @description Calculates the average element-centric similarity between two
#' clustering results
#'
#' @param clustering1 The first clustering result, which can be one of:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param clustering2 The second clustering result, which can be one of:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#' @param ppr_implementation_cl1 Choose a implementation for personalized
#' page-rank calculation for the first clustering:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize_cl1 Whether to normalize all rows in the first clustering
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r_cl1 A numeric hierarchical scaling parameter for the first clustering.
#' @param rescale_path_type_cl1 A string; rescale the hierarchical height of
#' the first clustering by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled_cl1 A logical: if TRUE, the linkage distances of the first
#' clustering are linearly rescaled to be in-between 0 and 1.
#' @param ppr_implementation_cl2 Choose a implementation for personalized
#' page-rank calculation for the second clustering:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize_cl2 Whether to normalize all rows in the second clustering
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r_cl2 A numeric hierarchical scaling parameter for the second clustering.
#' @param rescale_path_type_cl2 A string; rescale the hierarchical height of
#' the second clustering by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled_cl2 A logical: if TRUE, the linkage distances of the second
#' clustering are linearly rescaled to be in-between 0 and 1.
#'
#' @return The average element-wise similarity between the two Clusterings.
#' @export
#'
#' @examples
#' km.res <- kmeans(mtcars, centers = 3)$cluster
#' hc.res <- hclust(dist(mtcars))
#' element_sim(km.res, hc.res)
element_sim <- function(clustering1,
                        clustering2,
                        alpha = 0.9,
                        r_cl1 = 1,
                        rescale_path_type_cl1 = "max",
                        ppr_implementation_cl1 = "prpack",
                        dist_rescaled_cl1 = FALSE,
                        row_normalize_cl1 = TRUE,
                        r_cl2 = 1,
                        rescale_path_type_cl2 = "max",
                        ppr_implementation_cl2 = "prpack",
                        dist_rescaled_cl2 = FALSE,
                        row_normalize_cl2 = TRUE) {
    element_scores <- element_sim_elscore(
        clustering1,
        clustering2,
        alpha,
        r_cl1,
        rescale_path_type_cl1,
        ppr_implementation_cl1,
        dist_rescaled_cl1,
        row_normalize_cl1,
        r_cl2,
        rescale_path_type_cl2,
        ppr_implementation_cl2,
        dist_rescaled_cl2,
        row_normalize_cl2
    )
    return(mean(element_scores))
}

#' The Element-Centric Clustering Similarity for each Element
#'
#' @description Calculates the element-wise element-centric similarity between
#' two clustering results.
#'
#' @param clustering1 The first clustering result, which can be one of:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param clustering2 The second clustering result, which can be one of:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#' @param ppr_implementation_cl1 Choose a implementation for personalized
#' page-rank calculation for the first clustering:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize_cl1 Whether to normalize all rows in the first clustering
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r_cl1 A numeric hierarchical scaling parameter for the first clustering.
#' @param rescale_path_type_cl1 A string; rescale the hierarchical height of
#' the first clustering by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled_cl1 A logical: if TRUE, the linkage distances of the first
#' clustering are linearly rescaled to be in-between 0 and 1.
#' @param ppr_implementation_cl2 Choose a implementation for personalized
#' page-rank calculation for the second clustering:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize_cl2 Whether to normalize all rows in the second clustering
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r_cl2 A numeric hierarchical scaling parameter for the second clustering.
#' @param rescale_path_type_cl2 A string; rescale the hierarchical height of
#' the second clustering by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled_cl2 A logical: if TRUE, the linkage distances of the second
#' clustering are linearly rescaled to be in-between 0 and 1.
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
#' km.res <- kmeans(iris[, 1:4], centers = 8)$cluster
#' hc.res <- hclust(dist(iris[, 1:4]))
#' element_sim_elscore(km.res, hc.res)
element_sim_elscore <- function(clustering1,
                                clustering2,
                                alpha = 0.9,
                                r_cl1 = 1,
                                rescale_path_type_cl1 = "max",
                                ppr_implementation_cl1 = "prpack",
                                dist_rescaled_cl1 = FALSE,
                                row_normalize_cl1 = TRUE,
                                r_cl2 = 1,
                                rescale_path_type_cl2 = "max",
                                ppr_implementation_cl2 = "prpack",
                                dist_rescaled_cl2 = FALSE,
                                row_normalize_cl2 = TRUE) {
    # if both clusterings are membership vectors, calculate the ecs without
    # creating a Clustering object
    if (inherits(clustering1, c("numeric", "integer", "factor", "character")) &&
        inherits(clustering2, c("numeric", "integer", "factor", "character"))) {
        if (length(clustering1) != length(clustering2)) {
            stop("clustering1 and clustering2 do not have the same length.")
        }
        if (any(names(clustering1) != names(clustering2))) {
            stop("Not all elements of clustering1 and clustering2 are the same.")
        }

        if (inherits(clustering1, "character") || inherits(clustering2, "character")) {
            node.scores <- corrected_l1_mb(clustering1, clustering2, alpha)
        } else {
            node.scores <- disjointECS(
                clustering1,
                clustering2
            )
        }

        names(node.scores) <- names(clustering1)
        return(node.scores)
    }

    if (inherits(clustering1, "hclust")) {
        clustering1 <- create_clustering(
            clustering_result = clustering1,
            alpha = alpha,
            r = r_cl1,
            rescale_path_type = rescale_path_type_cl1,
            ppr_implementation = ppr_implementation_cl1,
            dist_rescaled = dist_rescaled_cl1
        )
    } else if (is.matrix(clustering1) | inherits(clustering1, "Matrix")) {
        clustering1 <- create_clustering(
            clustering_result = clustering1,
            alpha = alpha,
            ppr_implementation = ppr_implementation_cl1,
            row_normalize = row_normalize_cl1
        )
    } else {
        clustering1 <- create_clustering(clustering1,
            alpha = alpha
        )
    }

    if (inherits(clustering2, "hclust")) {
        clustering2 <- create_clustering(
            clustering_result = clustering2,
            alpha = alpha,
            r = r_cl2,
            rescale_path_type = rescale_path_type_cl2,
            ppr_implementation = ppr_implementation_cl2,
            dist_rescaled = dist_rescaled_cl2
        )
    } else if (is.matrix(clustering2) | inherits(clustering2, "Matrix")) {
        clustering2 <- create_clustering(
            clustering_result = clustering2,
            alpha = alpha,
            ppr_implementation = ppr_implementation_cl2,
            row_normalize = row_normalize_cl2
        )
    } else {
        clustering2 <- create_clustering(clustering2,
            alpha = alpha
        )
    }

    # Make sure clusterings are comparable
    if (clustering1@n_elements != clustering2@n_elements) {
        stop("clustering1 and clustering2 do not have the same length.")
    } else if (any(names(clustering1) != names(clustering2))) {
        stop("Not all elements of clustering1 and clustering2 are the same.")
    }

    # use the corrected L1 similarity
    node.scores <- corrected_L1(
        clustering1@affinity_matrix,
        clustering2@affinity_matrix,
        clustering1@alpha
    )

    names(node.scores) <- names(clustering1)
    return(node.scores)
}

# create the cluster -> element dictionary (list) for flat, disjoint clustering
create_clu2elm_dict <- function(clustering) {
    clu2elm_dict <- list()
    for (i in unique(clustering)) {
        clu2elm_dict[[i]] <- which(clustering == i)
    }

    return(clu2elm_dict)
}

# create the cluster -> element dictionary (list) for hierarchical clustering
create_clu2elm_dict_hierarchical <- function(clustering) {
    n.clusters <- nrow(clustering$merge)
    clu2elm_dict <- list()

    # add NA to initialize list
    clu2elm_dict[[2 * n.clusters + 1]] <- NA

    # add single member clusters for leaf nodes
    for (i in 1:(n.clusters + 1)) {
        clu2elm_dict[[i]] <- i
    }

    # traverse hc tree and add nodes to clu2elm_dict
    for (i in 1:n.clusters) {
        for (j in 1:2) {
            index <- clustering$merge[i, j]
            if (index < 0) { # for leaf nodes, just append them
                clu2elm_dict[[i + n.clusters + 1]] <- c(
                    clu2elm_dict[[i + n.clusters + 1]],
                    -index
                )
            } else { # for previous clusters, append all leaf nodes in cluster
                clu2elm_dict[[i + n.clusters + 1]] <- c(
                    clu2elm_dict[[i + n.clusters + 1]],
                    clu2elm_dict[[index + n.clusters + 1]]
                )
            }
        }
    }
    # remove the NA that was added in the beginning
    clu2elm_dict[[2 * n.clusters + 1]] <- clu2elm_dict[[2 * n.clusters + 1]][-1]

    return(clu2elm_dict)
}

# create element -> cluster dictionary (list) for flat, overlapping clustering
create_elm2clu_dict_overlapping <- function(clustering) {
    elm2clu_dict <- list()
    for (i in seq_len(nrow(clustering))) {
        elm2clu_dict[[i]] <- which(clustering[i, ] > 0)
    }
    return(elm2clu_dict)
}

# create element -> cluster dictionary (list) for hierarchical clustering: only
# maps elements to leaf node singleton clusters
create_elm2clu_dict_hierarchical <- function(clustering) {
    elm2clu_dict <- list()
    for (i in seq_along(clustering$order)) {
        elm2clu_dict[[i]] <- i
    }
    return(elm2clu_dict)
}

# corrected L1 distance
corrected_L1 <- function(x, y, alpha) {
    res <- 1 - 1 / (2 * alpha) * rowSums(abs(x - y))
    return(res)
}

# corrected L1 distance for two membership vectors
corrected_l1_mb_very_old <- function(mb1, mb2, alpha = 0.9) {
    n <- length(mb1)

    if (!is.character(mb1)) {
        mb1 <- as.character(mb1)
    }

    if (!is.character(mb2)) {
        mb2 <- as.character(mb2)
    }

    clu2elm_dict_1 <- create_clu2elm_dict(mb1)
    clu2elm_dict_2 <- create_clu2elm_dict(mb2)

    # the possible number of different ecs values is n x m
    # where n is the number of clusters of the first partition
    # and n the number of clusters of the second partition
    unique_ecs_vals <- matrix(NA, nrow = length(clu2elm_dict_1), ncol = length(clu2elm_dict_2))
    rownames(unique_ecs_vals) <- names(clu2elm_dict_1)
    colnames(unique_ecs_vals) <- names(clu2elm_dict_2)

    ecs <- rep(0, n)

    ppr1 <- rep(0, n)
    ppr2 <- rep(0, n)

    # iterate through each point of the membership vector
    for (i in 1:n) {
        # check if the similarity between of this pair of clusters was already calculated
        if (is.na(unique_ecs_vals[mb1[i], mb2[i]])) {
            clusterlist1 <- clu2elm_dict_1[[mb1[i]]]
            Csize1 <- length(clusterlist1)

            # get the values from the affinity matrix of the first partition
            ppr1[clusterlist1] <- alpha / Csize1
            ppr1[i] <- 1.0 - alpha + alpha / Csize1

            clusterlist2 <- clu2elm_dict_2[[mb2[i]]]
            Csize2 <- length(clusterlist2)

            # get the values from the affinity matrix of the second partition
            ppr2[clusterlist2] <- alpha / Csize2
            ppr2[i] <- 1.0 - alpha + alpha / Csize2

            # calculate the sum of the difference between the affinity matrices
            ecs[i] <- sum(abs(ppr1 - ppr2)) # could be optimized?

            ppr1[clusterlist1] <- 0
            ppr2[clusterlist2] <- 0

            # store the similarity between the pair of clusters
            unique_ecs_vals[mb1[i], mb2[i]] <- ecs[i]
        } else {
            # if yes, just copy the value
            ecs[i] <- unique_ecs_vals[mb1[i], mb2[i]]
        }
    }

    # perform the last calculations to obtain the ECS at each point
    return(1 - 1 / (2 * alpha) * ecs)
}

corrected_l1_mb_ <- function(mb1, mb2, alpha = 0.9) {
    n <- length(mb1)

    if (!is.character(mb1)) {
        mb1 <- as.character(mb1)
    }

    if (!is.character(mb2)) {
        mb2 <- as.character(mb2)
    }

    clu2elm_dict_1 <- create_clu2elm_dict(mb1)
    clu2elm_dict_2 <- create_clu2elm_dict(mb2)

    # the possible number of different ecs values is n x m
    # where n is the number of clusters of the first partition
    # and n the number of clusters of the second partition
    nclust1 <- length(clu2elm_dict_1)
    nclust2 <- length(clu2elm_dict_2)
    # unique_ecs_vals = matrix(NA, nrow = nclust1, ncol = nclust2)
    # rownames(unique_ecs_vals) = names(clu2elm_dict_1)
    # colnames(unique_ecs_vals) = names(clu2elm_dict_2)

    # cont_table = table(mb1, mb2)
    # order_clust1 = rownames(cont_table)
    # order_clust2 = colnames(cont_table)

    # ecs = vector(mode = "double", length = n)
    ecs <- rep(0, n)

    for (i in 1:nclust1) {
        cluster1 <- clu2elm_dict_1[[i]]
        size_cluster1 <- length(cluster1)
        for (j in 1:nclust2) {
            cluster2 <- clu2elm_dict_2[[j]]
            size_cluster2 <- length(cluster2)
            cluster_intersection <- match(cluster1, cluster2, nomatch = 0) > 0
            intersection_points <- cluster1[cluster_intersection]
            # cluster_intersection2 = match(cluster2, intersection_points, nomatch = 0) > 0
            # clu2elm_dict_2[[j]] = cluster2[!cluster_intersection2]

            # cluster_intersection = intersect(cluster1, cluster2)

            # cluster_points = cluster1[cluster_intersection]
            cluster1 <- cluster1[!cluster_intersection]
            # print(length(cluster1))

            intersection_size <- length(intersection_points)

            ecs_val <- intersection_size * (1 / size_cluster1 - 1 / size_cluster2)

            if (ecs_val < 0) {
                ecs_val <- -ecs_val
            }

            ecs_val <- ecs_val +
                (size_cluster1 - intersection_size) * 1 / size_cluster1 +
                (size_cluster2 - intersection_size) * 1 / size_cluster2

            ecs_val <- 1 - 1 / 2 * ecs_val

            ecs[intersection_points] <- ecs_val
        }
    }

    # iterate through each point of the membership vector
    # for(i in 1:n) {
    # ecs[i] = 1#unique_ecs_vals[mb1[i], mb2[i]]
    # }

    # ecs = vapply(1:n, function(i) { unique_ecs_vals[mb1[i], mb2[i]]}, double(1L))

    return(ecs)
}

corrected_l1_mb <- function(mb1, mb2, alpha = 0.9) {
    n <- length(mb1)
    ecs <- rep(0, n)

    cont_table <- table(mb1, mb2)
    rownames(cont_table) <- as.character(rownames(cont_table))
    colnames(cont_table) <- as.character(colnames(cont_table))
    c1_sizes <- rowSums(cont_table)
    c2_sizes <- colSums(cont_table)
    unique_ecs_vals <- matrix(NA, nrow = length(c1_sizes), ncol = length(c2_sizes))
    rownames(unique_ecs_vals) <- rownames(cont_table)
    colnames(unique_ecs_vals) <- colnames(cont_table)

    # iterate through each point of the membership vector
    for (i in 1:n) {
        # check if the similarity between of this pair of clusters was already calculated
        mb1_elem <- as.character(mb1[i])
        mb2_elem <- as.character(mb2[i])
        if (is.na(unique_ecs_vals[mb1_elem, mb2_elem])) {
            Csize1 <- c1_sizes[mb1_elem]
            Csize2 <- c2_sizes[mb2_elem]

            intersect_size <- cont_table[mb1_elem, mb2_elem]
            temp_ecs <- (1 / Csize1 - 1 / Csize2) * intersect_size
            if (temp_ecs < 0) {
                temp_ecs <- -temp_ecs
            }
            # calculate the sum of the difference between the affinity matrices
            ecs[i] <- 1 - 0.5 * (temp_ecs + (Csize1 - intersect_size) / Csize1 + (Csize2 - intersect_size) / Csize2)

            # store the similarity between the pair of clusters
            unique_ecs_vals[mb1_elem, mb2_elem] <- ecs[i]
        } else {
            # if yes, just copy the value
            ecs[i] <- unique_ecs_vals[mb1_elem, mb2_elem]
        }
    }

    return(ecs)
}

# PPR calculation for partition
ppr_partition <- function(clustering, alpha = 0.9) {
    ppr <- matrix(0, nrow = clustering@n_elements, ncol = clustering@n_elements)

    for (i in seq_along(clustering@clu2elm_dict)) {
        clusterlist <- clustering@clu2elm_dict[[i]]
        Csize <- length(clusterlist)
        ppr_res <- matrix(alpha / Csize, nrow = Csize, ncol = Csize)
        diag(ppr_res) <- 1.0 - alpha + alpha / Csize
        ppr[clusterlist, clusterlist] <- ppr_res
    }

    return(ppr)
}

# Create the cluster-induced element graph for overlapping clustering
make_cielg_overlapping <- function(bipartite_adj) {
    proj1 <- diag(x = 1 / rowSums(bipartite_adj), nrow = nrow(bipartite_adj)) %*%
        bipartite_adj

    proj2 <- bipartite_adj %*%
        diag(x = 1 / colSums(bipartite_adj), nrow = ncol(bipartite_adj))
    projected_adj <- proj1 %*% t(proj2)

    cielg <- igraph::graph_from_adjacency_matrix(projected_adj,
        weighted = TRUE,
        mode = "directed"
    )
    return(cielg)
}

# create linkage distances vector
get_linkage_dist <- function(hierarchical_clustering, dist_rescaled) {
    n.elements <- length(hierarchical_clustering$order)

    if (dist_rescaled) {
        maxdist <- max(hierarchical_clustering$height)
        # add 1s at beggining to account for leaf node singleton clusters
        linkage_dist <- c(
            rep(1, n.elements),
            1.0 - hierarchical_clustering$height / maxdist
        )
    } else {
        # add 0s at beginning to account for leaf node singleton clusters
        linkage_dist <- c(rep(0, n.elements), hierarchical_clustering$height)
    }
    return(linkage_dist)
}

# helper function to calculate path distances in graph
dist_path <- function(hierarchical_clustering,
                      rescale_path_type) {
    n.elements <- length(hierarchical_clustering$order)
    # dist to node from final cluster containing all elements
    path_to_node <- rep(0, 2 * n.elements - 1)
    for (i in rev(1:(n.elements - 1))) {
        for (j in 1:2) {
            index <- hierarchical_clustering$merge[i, j]
            if (index < 0) {
                index <- -index
            } else {
                index <- index + n.elements
            }
            if (rescale_path_type == "min") {
                if (path_to_node[index] == 0) {
                    path_to_node[index] <- path_to_node[i + n.elements] + 1
                } else {
                    path_to_node[index] <- min(
                        path_to_node[index],
                        path_to_node[i + n.elements] + 1
                    )
                }
            } else {
                path_to_node[index] <- max(
                    path_to_node[index],
                    path_to_node[i + n.elements] + 1
                )
            }
        }
    }

    # dist to node from leaves
    path_from_node <- rep(0, 2 * n.elements - 1)
    for (i in 1:(n.elements - 1)) {
        index <- hierarchical_clustering$merge[i, ]
        for (j in 1:2) {
            if (index[j] < 0) {
                index[j] <- -index[j]
            } else {
                index[j] <- index[j] + n.elements
            }
        }
        if (rescale_path_type == "min") {
            path_from_node[i + n.elements] <- min(
                path_from_node[index[1]],
                path_from_node[index[2]]
            ) + 1
        } else {
            path_from_node[i + n.elements] <- max(
                path_from_node[index[1]],
                path_from_node[index[2]]
            ) + 1
        }
    }

    return(list(to = path_to_node, from = path_from_node))
}

# rescale the hierarchical clustering path
rescale_path <- function(hierarchical_clustering,
                         rescale_path_type) {
    # get path distances
    path.dist <- dist_path(hierarchical_clustering, rescale_path_type)
    path_from_node <- path.dist$from
    path_to_node <- path.dist$to

    n.elements <- length(hierarchical_clustering$order)
    rescaled_level <- rep(0, n.elements)
    for (i in 1:(2 * n.elements - 1)) {
        total_path_length <- path_from_node[i] + path_to_node[i]
        if (total_path_length != 0) {
            rescaled_level[i] <- path_to_node[i] / total_path_length
        }
    }
    return(rescaled_level)
}

# Create the cluster-induced element graph for hierarchical clustering
make_cielg_hierarchical <- function(hierarchical_clustering,
                                    clu2elm_dict,
                                    r,
                                    rescale_path_type,
                                    dist_rescaled = FALSE) {
    n.elements <- length(hierarchical_clustering$order)

    # rescale paths
    if (rescale_path_type == "linkage") {
        cluster_height <- get_linkage_dist(
            hierarchical_clustering,
            dist_rescaled
        )
    } else if (rescale_path_type %in% c("min", "max")) {
        cluster_height <- rescale_path(
            hierarchical_clustering,
            rescale_path_type
        )
    } else {
        stop(paste0("rescale_path_type must be one of linkage, min or max."))
    }

    # weight function for different heights
    weight_function <- function(clust) {
        return(exp(r * (cluster_height[clust])))
    }

    bipartite_adj <- matrix(0, nrow = n.elements, ncol = 2 * n.elements - 1)
    for (i in seq_along(clu2elm_dict)) {
        element_list <- clu2elm_dict[[i]]

        cstrength <- weight_function(i)
        for (el in element_list) {
            bipartite_adj[el, i] <- cstrength
        }
    }

    proj1 <- diag(x = 1 / rowSums(bipartite_adj), nrow = nrow(bipartite_adj)) %*%
        bipartite_adj
    proj2 <- bipartite_adj %*%
        diag(x = 1 / colSums(bipartite_adj), nrow = ncol(bipartite_adj))
    projected_adj <- proj1 %*% t(proj2)

    cielg <- igraph::graph_from_adjacency_matrix(projected_adj,
        weighted = TRUE,
        mode = "directed"
    )
    return(cielg)
}

# utility function to find vertices with the same cluster memberships.
find_groups_in_cluster <- function(clustervs, elementgroupList) {
    clustervertex <- unique(clustervs)

    groupings <- list()
    index <- 1
    for (i in seq_along(elementgroupList)) {
        current <- elementgroupList[[i]]
        if (length(intersect(current, clustervertex)) > 0) {
            groupings[[index]] <- current
            index <- index + 1
        }
    }

    return(groupings)
}

# numerical calculation of affinity matrix using PPR
numerical_ppr_scores <- function(cielg, clustering, ppr_implementation = "prpack") {
    if (!(ppr_implementation %in% c("prpack", "power_iteration"))) {
        stop("ppr_implementation argument must be one of prpack or power_iteration")
    }

    # keep track of all clusters an element is a member of
    elementgroupList <- list()
    for (i in seq_along(clustering@elm2clu_dict)) {
        element <- names(clustering)[i]
        # collapse clusters into single string
        clusters <- paste(clustering@elm2clu_dict[[i]], collapse = ",")

        elementgroupList[[clusters]] <- c(elementgroupList[[clusters]], element)
    }

    ppr_scores <- matrix(0,
        nrow = igraph::vcount(cielg),
        ncol = igraph::vcount(cielg)
    )

    # we have to calculate the ppr for each connected component
    comps <- igraph::components(cielg)
    for (i.comp in 1:comps$no) {
        members <- which(comps$membership == i.comp)
        clustergraph <- igraph::induced_subgraph(cielg, members, impl = "auto")
        cc_ppr_scores <- matrix(0,
            nrow = igraph::vcount(clustergraph),
            ncol = igraph::vcount(clustergraph)
        )

        if (ppr_implementation == "power_iteration") {
            W_matrix <- get_sparse_transition_matrix(clustergraph)
        }

        elementgroups <- find_groups_in_cluster(members, elementgroupList)
        for (elementgroup in elementgroups) {
            # we only have to solve for the ppr distribution once per group
            vertex.index <- match(elementgroup[1], members)
            if (ppr_implementation == "prpack") {
                pers <- rep(0, length = length(members))
                pers[vertex.index] <- 1
                ppr_res <- igraph::page_rank(clustergraph,
                    directed = TRUE,
                    weights = NULL,
                    damping = clustering@alpha,
                    personalized = pers,
                    algo = "prpack"
                )
                cc_ppr_scores[vertex.index, ] <- ppr_res$vector
            } else if (ppr_implementation == "power_iteration") {
                cc_ppr_scores[vertex.index, ] <- calculate_ppr_with_power_iteration(
                    W_matrix,
                    vertex.index,
                    alpha = clustering@alpha,
                    repetition = 1000,
                    th = 0.0001
                )
            }

            # the other vertices in the group are permutations of that solution
            elgroup.size <- length(elementgroup)
            if (elgroup.size >= 2) {
                for (i2 in 2:elgroup.size) {
                    v2 <- match(elementgroup[i2], members)
                    cc_ppr_scores[v2, ] <- cc_ppr_scores[vertex.index, ]
                    cc_ppr_scores[v2, vertex.index] <- cc_ppr_scores[vertex.index, v2]
                    cc_ppr_scores[v2, v2] <- cc_ppr_scores[vertex.index, vertex.index]
                }
            }
        }

        ppr_scores[members, members] <- cc_ppr_scores
    }
    return(ppr_scores)
}

# utility function to get a row-normalized sparse transition matrix
get_sparse_transition_matrix <- function(graph) {
    transition_matrix <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = FALSE)
    transition_matrix <- diag(x = 1 / rowSums(transition_matrix), nrow = nrow(transition_matrix)) %*%
        transition_matrix
    return(transition_matrix)
}

# power iteration calculation of PPR
calculate_ppr_with_power_iteration <- function(W_matrix, index, alpha = 0.9,
                                               repetition = 1000, th = 1e-4) {
    total_length <- nrow(W_matrix)
    e_s <- matrix(0, nrow = 1, ncol = total_length)
    e_s[1, index] <- 1
    p <- matrix(0, nrow = 1, ncol = total_length)
    p[1, index] <- 1
    for (i in 1:repetition) {
        new_p <- ((1 - alpha) * e_s) + ((alpha) * (p %*% W_matrix))
        if (max(abs(new_p - p)) < th) {
            p <- new_p
            break
        }
        p <- new_p
    }

    return(p)
}

# TODO idea: create a similarity graph given a list of clusterings; the purpose: identifying the partitions that relate the most with the others; caveats: the same as in the ecs threshold: if A and B are similar, and B and C are similar, that doesn't mean that A and C are similar

#' Pairwise Comparison of Clusterings
#' @description Compare a set of clusterings by calculating their pairwise
#' average element-centric clustering similarities.
#'
#' @param clustering_list The list of clustering results, each of which is either:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param output_type A string specifying whether the output should be a
#' matrix or a data.frame.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#' @param ppr_implementation Choose a implementation for personalized
#' page-rank calculation:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize Whether to normalize all rows in clustering_result
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r A numeric hierarchical scaling parameter.
#' @param rescale_path_type A string; rescale the hierarchical height by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled A logical: if TRUE, the linkage distances are linearly
#' rescaled to be in-between 0 and 1.
#'
#' @return A matrix or data.frame containing the pairwise ECS values.
#' @export
#'
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
#'
#' @examples
#' # cluster across 20 random seeds
#' clustering.list <- lapply(1:20, function(x) kmeans(mtcars, centers = 3)$cluster)
#' element_sim_matrix(clustering.list, output_type = "matrix")
element_sim_matrix <- function(clustering_list,
                               output_type = "matrix",
                               alpha = 0.9,
                               r = 1,
                               rescale_path_type = "max",
                               ppr_implementation = "prpack",
                               dist_rescaled = FALSE,
                               row_normalize = TRUE) {
    if (!(output_type %in% c("data.frame", "matrix"))) {
        stop("output_type must be data.frame or matrix.")
    }

    # check if all objects are flat disjoint membership vectors
    are_all_flat_disjoint <- sapply(clustering_list, function(x) {
        inherits(x, c("numeric", "integer", "factor", "character"))
    })

    are_all_viable_objects <- sapply(clustering_list, function(x) {
        is.matrix(x) || inherits(x, "hclust") || inherits(x, "Matrix")
    })

    if (any(are_all_flat_disjoint == FALSE) &&
        any(are_all_viable_objects == FALSE)) {
        stop("Please ensure each entry in clustering_list is from these following classes: numeric, factor, character, matrix, Matrix, hclust.")
    }

    # if the condition is met, perform element frustration using only the membership vector
    if (all(are_all_flat_disjoint == TRUE)) {
        return(element_sim_matrix_flat_disjoint(
            mb_list = clustering_list,
            alpha = alpha,
            output_type = output_type
        ))
    }

    clustering_list <- create_clustering_list(
        object_list = clustering_list,
        alpha = alpha,
        r = r,
        rescale_path_type = rescale_path_type,
        ppr_implementation = ppr_implementation,
        dist_rescaled = dist_rescaled,
        row_normalize = row_normalize
    )

    n.clusterings <- length(clustering_list)
    sim_matrix <- matrix(NA, nrow = n.clusterings, ncol = n.clusterings)
    for (i in 1:(n.clusterings - 1)) {
        i.aff <- clustering_list[[i]]@affinity_matrix
        for (j in (i + 1):n.clusterings) {
            j.aff <- clustering_list[[j]]@affinity_matrix

            sim_matrix[i, j] <- mean(corrected_L1(i.aff, j.aff, alpha))
        }
    }
    diag(sim_matrix) <- 1

    if (output_type == "data.frame") {
        sim_matrix <- data.frame(
            i.clustering = rep(
                1:n.clusterings,
                n.clusterings
            ),
            j.clustering = rep(1:n.clusterings,
                each = n.clusterings
            ),
            element_sim = c(sim_matrix)
        )
    }

    return(sim_matrix)
}

# calculate the similarity matrix when all the objects are flat disjoint memberships
element_sim_matrix_flat_disjoint <- function(mb_list, alpha = 0.9, output_type = "matrix") {
    if (!(output_type %in% c("data.frame", "matrix"))) {
        stop("output_type must be data.frame or matrix.")
    }

    n_clusterings <- length(mb_list)

    if (n_clusterings == 1) {
        return(matrix(1, nrow = 1, ncol = 1))
    }

    first_index <- unlist(sapply(1:(n_clusterings - 1), function(i) {
        rep(i, n_clusterings - i)
    }))
    second_index <- unlist(sapply(1:(n_clusterings - 1), function(i) {
        (i + 1):n_clusterings
    }))
    n_combinations <- n_clusterings * (n_clusterings - 1) / 2

    i <- NA

    # calculate the ecs between every pair of partitions
    ecs_values <- foreach::foreach(i = 1:n_combinations, .export = c("corrected_l1_mb", "create_clu2elm_dict"), .noexport = c("my_cluster"), .combine = "c") %dopar% {
        mean(element_sim_elscore(
            mb_list[[first_index[i]]],
            mb_list[[second_index[i]]],
            alpha
        ))
    }

    sim_matrix <- matrix(NA, nrow = n_clusterings, ncol = n_clusterings)
    sim_matrix[lower.tri(sim_matrix, diag = FALSE)] <- ecs_values
    sim_matrix <- t(sim_matrix)
    diag(sim_matrix) <- 1

    if (output_type == "data.frame") {
        sim_matrix <- data.frame(
            i.clustering = rep(
                1:n_clusterings,
                n_clusterings
            ),
            j.clustering = rep(1:n_clusterings,
                each = n_clusterings
            ),
            element_sim = c(sim_matrix)
        )
    }

    return(sim_matrix)
}

# determine whether two flat disjoint membership vectors are identical
are_identical_memberships <- function(mb1, mb2) {
    # calculate the contingency table between the two partitions
    if (inherits(mb1, c("character", "factor")) || inherits(mb2, c("character", "factor"))) {
        contingency_table <- table(mb1, mb2)
    } else {
        contingency_table <- myContTable(mb1, mb2, min(mb1), min(mb2))
    }

    contingency_table <- (contingency_table != 0)

    if (nrow(contingency_table) != ncol(contingency_table)) {
        return(FALSE)
    }

    # there are only singleton clusters
    if (nrow(contingency_table) == length(mb1)) {
        return(TRUE)
    }

    no_different_elements <- colSums(contingency_table)

    # if a column has two or more nonzero entries, it means the cluster is split
    # thus the two partitions are not identical
    if (any(no_different_elements != 1)) {
        return(FALSE)
    }

    no_different_elements <- rowSums(contingency_table)

    # if a row has two or more nonzero entries, it means the cluster is split
    # thus the two partitions are not identical
    return(all(no_different_elements == 1))
}

# merge partitions that are identical (meaning the ecs threshold is 1)
# the order parameter is indicating whether to sort the merged objects
# based on their frequency
merge_identical_partitions <- function(clustering_list) {
    # merge the same memberships into the same object
    merged_partitions <- list(clustering_list[[1]])
    n_partitions <- length(clustering_list)

    if (n_partitions == 1) {
        return(clustering_list)
    }

    for (i in 2:n_partitions) {
        n_merged_partitions <- length(merged_partitions)

        assigned_partition <- -1

        # for each partition, look into the merged list and see if there is a perfect match
        for (j in seq_len(n_merged_partitions)) {
            are_identical <- are_identical_memberships(merged_partitions[[j]]$mb, clustering_list[[i]]$mb)
            if (are_identical) {
                assigned_partition <- j
                break
            }
        }

        if (assigned_partition != -1) {
            # if a match was found, update the frequency of the partition
            merged_partitions[[assigned_partition]]$freq <- merged_partitions[[assigned_partition]]$freq + clustering_list[[i]]$freq
        } else {
            # if not, add the partition to the list of merged partitions
            merged_partitions[[n_merged_partitions + 1]] <- clustering_list[[i]]
        }
    }

    merged_partitions
}

# merge the partitions when the ecs threshold is not 1
merge_partitions_ecs <- function(partition_list, ecs_thresh = 0.99) {
    partition_groups <- list()
    nparts <- length(partition_list)

    if (nparts == 1) {
        return(partition_list)
    }

    # calculate the pairwise ecs between the partitions
    sim_matrix <- element_sim_matrix_flat_disjoint(
        lapply(partition_list, function(x) {
            x$mb
        })
    )

    for (i in 1:nparts) {
        sim_matrix[i, i] <- NA
        partition_groups[[as.character(i)]] <- i
    }

    while (length(partition_groups) > 1) {
        # if the similarity between any pair of partitions is below the threshold
        # we do not perform any merging
        if (max(sim_matrix, na.rm = TRUE) < ecs_thresh) {
            break
        }

        # find the pair of partitions that are the most similar wrt ECS
        index <- which.max(sim_matrix)[1]
        first_cluster <- index %% nparts

        second_cluster <- index %/% nparts + 1
        if (first_cluster == 0) {
            first_cluster <- nparts
            second_cluster <- second_cluster - 1
        }

        # update the partition group of the first partition, where we add the second one
        partition_groups[[as.character(first_cluster)]] <- c(
            partition_groups[[as.character(first_cluster)]],
            partition_groups[[as.character(second_cluster)]]
        )
        # remove the partition group of the second partition that will be merged
        partition_groups[[as.character(second_cluster)]] <- NULL

        # update the similarities in a single-linkeage fashion
        # iterate only through the names of the groups of partitions, as the others
        # are already filled with NAs
        for (i in as.numeric(names(partition_groups))) {
            if (first_cluster < i) {
                if (!is.na(sim_matrix[first_cluster, i])) {
                    sim_matrix[first_cluster, i] <- min(sim_matrix[first_cluster, i], sim_matrix[second_cluster, i], sim_matrix[i, second_cluster], na.rm = T)
                }
            } else {
                if (!is.na(sim_matrix[i, first_cluster])) {
                    sim_matrix[i, first_cluster] <- min(sim_matrix[i, first_cluster], sim_matrix[second_cluster, i], sim_matrix[i, second_cluster], na.rm = T)
                }
            }
        }

        sim_matrix[second_cluster, ] <- NA
        sim_matrix[, second_cluster] <- NA
    }

    merged_index <- 1
    merged_partitions <- list()

    # update the frequencies of the partition groups
    for (kept_partition in names(partition_groups)) {
        merged_partitions[[merged_index]] <- partition_list[[as.numeric(kept_partition)]]
        merged_partitions[[merged_index]]$freq <- sum(sapply(partition_groups[[kept_partition]], function(x) {
            partition_list[[as.numeric(x)]]$freq
        }))

        merged_index <- merged_index + 1
    }

    merged_partitions
}

#' Merge Partitions
#' @description Merge flat disjoint clusterings whose pairwise ECS score is
#' above a given threshold. The merging is done using a complete linkage approach.
#'
#' @param partition_list A list of flat disjoint membership vectors.
#' @param ecs_thresh A numeric: the ecs threshold.
#' @param order_logic Variable indicating the method of ordering the partitions.
#' It can take these three values:
#' - "freq": order the partitions based on their frequencies. The partition with
#' the highest frequency will be the first on the list (default).
#' - "avg_agreement": order the partitions based on their average agreement index.
#' The average agreement index of a partition is calculated as the mean of
#' the ECS scores between that partition and the other partitions from the list.
#' The partition with the highest agreement will be  the first on the list.
#' - "none": do not perform any ordering (not recommended). If selected, the
#' average agreement scores will not be calculated.
#' @param return_ecs_matrix A logical: if TRUE, the function will add the
#' ECS matrix to the return list. Defaults to FALSE.
#' @param check_ties A logical value that indicates whether to check for ties
#' in the highest frequency partitions or not. If TRUE, the function will put
#' at the first position the partition that has the highest similarity
#' with the other partitions. Defaults to `FALSE`.
#' @return a list of the merged partitions, together with their associated
#' ECC score. If `return_ecs_matrix` is set to TRUE, the function will also
#' return the ECS matrix.
#' @export
#' @examples
#' initial_list <- list(c(1, 1, 2), c(2, 2, 2), c("B", "B", "A"))
#' merge_partitions(initial_list, 1)
merge_partitions <- function(partition_list,
                             ecs_thresh = 1,
                             order_logic = c("freq", "avg_agreement", "none"),
                             return_ecs_matrix = FALSE,
                             check_ties = TRUE) {
    # check the parameters
    if (!is.numeric(ecs_thresh) || length(ecs_thresh) > 1) {
        stop("ecs_thresh parameter should be numeric")
    }

    order_logic <- order_logic[1]
    if (!(order_logic %in% c("freq", "avg_agreement", "none"))) {
        stop("order parameter should be either `freq`, `avg_agreement` or `none`")
    }

    # check the type of object that is provided in the list
    if (!inherits(partition_list[[1]], "list")) { # the elements should be membership vectors
        partition_list <- lapply(partition_list, function(x) {
            list(
                mb = x,
                freq = 1
            )
        })
    } else {
        # the list contains the partition field, created by `get_resolution_importance` method
        if ("partitions" %in% names(partition_list)) {
            return(merge_partitions(partition_list$partitions, ecs_thresh, order_logic))
        }
        # the case with list of lists of partitions
        if (!all(c("mb", "freq") %in% names(partition_list[[1]]))) {
            return(lapply(partition_list, function(sublist) {
                merge_partitions(sublist, ecs_thresh, order_logic)
            }))
        }
    }

    # if the threshold is 1, apply the identical merging
    if (ecs_thresh == 1) {
        part_list <- merge_identical_partitions(partition_list)
    } else {
        # otherwise merge the partitions using the `merge_partitions_ecs` function
        part_list <- merge_partitions_ecs(partition_list, ecs_thresh)
    }

    if (length(part_list) == 1) {
        return(list(partitions = part_list, ecc = rep(1, length(part_list[[1]]$mb))))
    }

    consistency <- weighted_element_consistency(
        clustering_list = lapply(part_list, function(x) {
            x$mb
        }),
        weights = sapply(part_list, function(x) {
            x$freq
        }),
        calculate_sim_matrix = TRUE
    )

    if (!(order_logic %in% c("freq", "avg_agreement")) && !return_ecs_matrix) {
        return(list(partitions = part_list, ecc = consistency$ecc))
    }

    if (!(order_logic %in% c("freq", "avg_agreement")) && return_ecs_matrix) {
        return(list(partitions = part_list, ecc = consistency$ecc, ecs_matrix = consistency$ecs_matrix))
    }

    # add average aggreement index for each membership - how similar it is
    # with the other clusterings
    n_unique_parts <- nrow(consistency$ecs_matrix)
    n_total_parts <- sum(sapply(part_list, function(x) x$freq))
    average_agreement <- sapply(
        seq_len(n_unique_parts),
        function(i) {

            total_sim <- (partition_list[[i]]$freq - 1) / (n_total_parts - 1)
            for (j in seq_len(ncol(consistency$ecs_matrix))) {
                if (i == j) {
                    next
                }

                total_sim <- total_sim + consistency$ecs_matrix[i, j] * partition_list[[j]]$freq / (n_total_parts - 1)
            }

            return(total_sim)
        }
    )

    for (i in seq_len(nrow(consistency$ecs_matrix))) {
        part_list[[i]]$avg_agreement <- average_agreement[i]
    }

    # order the newly merged partitions based on their frequencies / agreement
    ordered_indices <- order(sapply(part_list, function(x) {
        x[[order_logic]]
    }), decreasing = TRUE)
    part_list <- part_list[ordered_indices]

    if (return_ecs_matrix) {
        consistency$ecs_matrix <- consistency$ecs_matrix[ordered_indices, ordered_indices]
    }

    check_ties <- check_ties && length(part_list) > 1

    if (!check_ties && !return_ecs_matrix) {
        return(list(partitions = part_list, ecc = consistency$ecc))
    }

    if (!check_ties && return_ecs_matrix) {
        return(list(partitions = part_list, ecc = consistency$ecc, ecs_matrix = consistency$ecs_matrix))
    }

    second_key <- setdiff(c("freq", "avg_agreement"), order_logic)
    reference_value <- part_list[[1]][[order_logic]]

    for (i in seq(from = 2, to = length(part_list), by = 1)) {
        if (part_list[[i]][[order_logic]] != reference_value) {
            break
        }
    }

    i <- i - 1
    if (i > 1) {
        ordered_indices <- order(sapply(part_list[seq_len(i)], function(x) {
            x[[second_key]]
        }), decreasing = TRUE)

        part_list[seq_len(i)] <- part_list[ordered_indices]
    }

    if (!return_ecs_matrix) {
        return(list(partitions = part_list, ecc = consistency$ecc))
    }

    if (i > 1) {
        consistency$ecs_matrix[seq_len(i), seq_len(i)] <- consistency$ecs_matrix[ordered_indices, ordered_indices]
    }

    return(list(partitions = part_list, ecc = consistency$ecc, ecs_matrix = consistency$ecs_matrix))
}

#' Element-Wise Consistency Between a Set of Clusterings
#' @description Inspect the consistency of a set of clusterings by calculating
#' their element-wise clustering consistency (also known as element-wise frustration).
#'
#' @param clustering_list The list of clustering results, each of which is either:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#' @param ppr_implementation Choose a implementation for personalized
#' page-rank calculation:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize Whether to normalize all rows in clustering_result
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r A numeric hierarchical scaling parameter.
#' @param rescale_path_type A string; rescale the hierarchical height by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled A logical: if TRUE, the linkage distances are linearly
#' rescaled to be in-between 0 and 1.
#'
#' @return A vector containing the element-wise consistency. If
#' `calculate_sim_matrix` is set to `TRUE`, the element similarity matrix
#' will be returned as well.
#' @export
#'
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
#' @examples
#' # cluster across 20 random seeds
#' clustering.list <- lapply(1:20, function(x) kmeans(mtcars, centers = 3)$cluster)
#' element_consistency(clustering.list)
element_consistency <- function(clustering_list,
                                alpha = 0.9,
                                r = 1,
                                rescale_path_type = "max",
                                ppr_implementation = "prpack",
                                dist_rescaled = FALSE,
                                row_normalize = TRUE) {
    # check if all objects are flat disjoint membership vectors
    are_all_flat_disjoint <- sapply(clustering_list, function(x) {
        inherits(x, c("numeric", "integer", "factor", "character"))
    })

    are_all_viable_objects <- sapply(clustering_list, function(x) {
        is.matrix(x) | inherits(x, "Matrix") | inherits(x, "hclust")
    })

    if (any(are_all_flat_disjoint == FALSE) &&
        any(are_all_viable_objects == FALSE)) {
        stop("Please ensure each entry in clustering_list is from these following classes: numeric, factor, character, matrix, Matrix, hclust.")
    }

    # if the condition is met, perform element consistency using only the membership vector
    if (all(are_all_flat_disjoint == TRUE)) {
        # merge the identical partitions into the same object
        final_clustering_list <- merge_partitions(
            partition_list = clustering_list,
            ecs_thresh = 1,
            order_logic = "none",
            return_ecs_matrix = FALSE
        )
        
        return(final_clustering_list$ecc)
    }

    clustering_list <- create_clustering_list(
        object_list = clustering_list,
        alpha = alpha,
        r = r,
        rescale_path_type = rescale_path_type,
        ppr_implementation = ppr_implementation,
        dist_rescaled = dist_rescaled,
        row_normalize = row_normalize
    )

    # calculate the consistency between the Clustering objects
    n.clusterings <- length(clustering_list)
    consistency <- rep(0, length(clustering_list[[1]]))
    for (i in 1:(n.clusterings - 1)) {
        i.aff <- clustering_list[[i]]@affinity_matrix
        for (j in (i + 1):n.clusterings) {
            j.aff <- clustering_list[[j]]@affinity_matrix

            consistency <- consistency + corrected_L1(i.aff, j.aff, alpha)
        }
    }
    consistency <- consistency / (n.clusterings * (n.clusterings - 1) / 2)
    return(consistency)
}

#' Weighted Element-Centric Consistency
#' @description Calculate the weighted element-centric consistency of a set of
#' clusterings. The weights are used to give more importance to some clusterings
#' over others.
#'
#' @param clustering_list The list of clustering results, each of which is either:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param weights A numeric vector of weights for each clustering in
#' `clustering_list`. If `NULL`, then all weights will be equal to 1. Defaults
#' to `NULL`.
#' @param calculate_sim_matrix A logical value that indicates whether to
#' calculate the similarity matrix or not along with the consistency score.
#' Defaults to `FALSE`.
#'
#' @return A vector containing the weighted element-wise consistency. If
#' `calculate_sim_matrix` is set to `TRUE`, the element similarity matrix
#' will be returned as well.
#'
#' @note The weighted ECC will be calculated as \eqn{\displaystyle \frac{\sum_{i} \sum_{j} w_i w_j ECS(i, j)}{\sum_{i} w_i}}
#'
#' @export
#' @examples
#' # cluster across 20 random seeds
#' clustering_list <- lapply(1:20, function(x) kmeans(mtcars, centers = 3)$cluster)
#' weights <- sample(1:10, 20, replace = TRUE)
#' weighted_element_consistency(clustering_list, weights = weights)
weighted_element_consistency <- function(clustering_list,
                                         weights = NULL,
                                         calculate_sim_matrix = FALSE) {
    n_clusterings <- length(clustering_list)

    if (calculate_sim_matrix) {
        ecs_sim_matrix <- matrix(0, nrow = n_clusterings, ncol = n_clusterings)
    }

    if (n_clusterings == 1) {
        if (calculate_sim_matrix) {
            ecs_matrix[1, 1] <- 1
            return(list(ecc = rep(1, length(clustering_list[[1]])), ecs_matrix = ecs_sim_matrix))
        }

        return(rep(1, length(clustering_list[[1]])))
    }

    if (is.null(weights)) {
        weights <- rep(1, n_clusterings)
    }

    consistency <- 0
    total_weights <- sum(weights)
    norm_factor <- total_weights * (total_weights - 1) / 2
    for (i in seq_len(n_clusterings - 1)) {
        # add the ecc of 1 for identical partitions
        consistency <- consistency + weights[i] / norm_factor * (weights[i] - 1) / 2

        for (j in seq(from = i + 1, to = n_clusterings)) {
            current_ecs <- element_sim_elscore(clustering_list[[i]], clustering_list[[j]])
            if (calculate_sim_matrix) {
                ecs_sim_matrix[i, j] <- mean(current_ecs)
                ecs_sim_matrix[j, i] <- ecs_sim_matrix[i, j]
            }

            consistency <- consistency + current_ecs * weights[i] / norm_factor * weights[j]
        }
        
        if (calculate_sim_matrix) {
            ecs_sim_matrix[i, i] <- 1
        }
    }
    consistency <- consistency + weights[n_clusterings] / norm_factor * (weights[n_clusterings] - 1) / 2

    if (!calculate_sim_matrix) {
        return(consistency)
    }

    ecs_sim_matrix[n_clusterings, n_clusterings] <- 1

    return(list(ecc = consistency, ecs_matrix = ecs_sim_matrix))
}

#' Element-Wise Average Agreement Between a Set of Clusterings
#' @description Inspect how consistently of a set of clusterings agree with
#' a reference clustering by calculating their element-wise average agreement.
#'
#' @param reference_clustering The reference clustering, that each clustering in
#' clustering_list is compared to. It can be either:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param clustering_list The list of clustering results, each of which is either:
#' * A numeric/character/factor vector of cluster labels for each element.
#' * A samples x clusters matrix/Matrix::Matrix of nonzero membership values.
#' * An hclust object.
#' @param alpha A numeric giving the personalized PageRank damping factor;
#' 1 - alpha is the restart probability for the PPR random walk.
#' @param ppr_implementation Choose a implementation for personalized
#' page-rank calculation:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize Whether to normalize all rows in clustering_result
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r A numeric hierarchical scaling parameter.
#' @param rescale_path_type A string; rescale the hierarchical height by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled A logical: if TRUE, the linkage distances are linearly
#' rescaled to be in-between 0 and 1.
#'
#' @return A vector containing the element-wise average agreement.
#' @export
#'
#' @references Gates, A. J., Wood, I. B., Hetrick, W. P., & Ahn, Y. Y. (2019).
#' Element-centric clustering comparison unifies overlaps and hierarchy.
#' Scientific reports, 9(1), 1-13. https://doi.org/10.1038/s41598-019-44892-y
#'
#' @examples
#' # perform k-means clustering across 20 random seeds
#' reference.clustering <- iris$Species
#' clustering.list <- lapply(1:20, function(x) kmeans(iris[, 1:4], centers = 3)$cluster)
#' element_agreement(reference.clustering, clustering.list)
element_agreement <- function(reference_clustering,
                              clustering_list,
                              alpha = 0.9,
                              r = 1,
                              rescale_path_type = "max",
                              ppr_implementation = "prpack",
                              dist_rescaled = FALSE,
                              row_normalize = TRUE) {
    # check if all objects are flat disjoint membership vectors
    are_all_flat_disjoint <- sapply(clustering_list, function(x) {
        inherits(x, c("numeric", "integer", "factor", "character"))
    })

    are_all_viable_objects <- sapply(clustering_list, function(x) {
        is.matrix(x) | inherits(x, "Matrix") | inherits(x, "hclust")
    })

    # if the condition is met, perform element frustration using only the membership vector
    if (all(are_all_flat_disjoint == TRUE) &&
        inherits(reference_clustering, c("numeric", "integer", "factor", "character"))) {
        return(element_agreement_flat_disjoint(
            reference_clustering = reference_clustering,
            clustering_list = clustering_list,
            alpha = alpha
        ))
    }

    if (any(are_all_flat_disjoint == FALSE &&
        are_all_viable_objects == FALSE) ||
        !(inherits(reference_clustering, c("numeric", "integer", "factor", "character", "matrix", "Matrix", "hclust")))) {
        stop("Please ensure reference_clustering and each entry in clustering_list is from these following classes: numeric, factor, character, matrix, Matrix, hclust.")
    }

    clustering_list <- create_clustering_list(
        object_list = clustering_list,
        alpha = alpha,
        r = r,
        rescale_path_type = rescale_path_type,
        ppr_implementation = ppr_implementation,
        dist_rescaled = dist_rescaled,
        row_normalize = row_normalize
    )

    # create Clustering object for the reference clustering
    if (inherits(reference_clustering, "hclust")) {
        reference_clustering <- create_clustering(
            clustering_result = reference_clustering,
            alpha = alpha,
            r = r,
            rescale_path_type = rescale_path_type,
            ppr_implementation = ppr_implementation,
            dist_rescaled = dist_rescaled
        )
    } else if (is.matrix(reference_clustering) | inherits(reference_clustering, "Matrix")) {
        reference_clustering <- create_clustering(
            clustering_result = reference_clustering,
            alpha = alpha,
            ppr_implementation = ppr_implementation,
            row_normalize = row_normalize
        )
    } else {
        reference_clustering <- create_clustering(
            clustering_result = reference_clustering,
            alpha = alpha
        )
    }

    # calculate the average element agreement between the clustering list and the reference
    n.clusterings <- length(clustering_list)
    avg_agreement <- rep(0, length(clustering_list[[1]]))
    ref_aff <- reference_clustering@affinity_matrix
    for (i in 1:n.clusterings) {
        i.aff <- clustering_list[[i]]@affinity_matrix
        avg_agreement <- avg_agreement + corrected_L1(ref_aff, i.aff, alpha)
    }
    avg_agreement <- avg_agreement / n.clusterings
    return(avg_agreement)
}

# calculate the element agreement when the clustering list and the reference
# are flat disjoint membership vectors
element_agreement_flat_disjoint <- function(reference_clustering,
                                            clustering_list,
                                            alpha = 0.9) {
    n_points <- length(reference_clustering)
    for (mb_vector in clustering_list) {
        if (length(mb_vector) != n_points) {
            stop("The partitions do not have the same length")
        }
    }

    obj <- NA

    # calculate the average of the ecs between one member of the clustering list
    # and the reference clustering
    avg_agreement <- foreach::foreach(
        obj = clustering_list,
        .export = c("corrected_l1_mb", "create_clu2elm_dict"),
        .combine = "+"
    ) %do% {
        corrected_l1_mb(reference_clustering, obj, alpha)
    }

    return(avg_agreement / length(clustering_list))
}

# create a list of Clustering objects
create_clustering_list <- function(object_list,
                                   alpha = 0.9,
                                   r = 1,
                                   rescale_path_type = "max",
                                   ppr_implementation = "prpack",
                                   dist_rescaled = FALSE,
                                   row_normalize = TRUE) {
    obj <- NA
    clustering_list <- foreach::foreach(
        obj = object_list,
        .packages = "ClustAssess",
        .export = "create_clustering"
    ) %dopar% {
        if (inherits(obj, "hclust")) {
            return(create_clustering(
                clustering_result = obj,
                alpha = alpha,
                r = r,
                rescale_path_type = rescale_path_type,
                ppr_implementation = ppr_implementation,
                dist_rescaled = dist_rescaled
            ))
        }

        if (is.matrix(obj) | inherits(obj, "Matrix")) {
            return(create_clustering(
                clustering_result = obj,
                alpha = alpha,
                ppr_implementation = ppr_implementation,
                row_normalize = row_normalize
            ))
        }

        create_clustering(
            clustering_result = obj,
            alpha = alpha
        )
    }

    return(clustering_list)
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
#' @keywords internal
Clustering <- setClass(
    "Clustering",
    representation(
        names = "character",
        n_elements = "numeric",
        is_hierarchical = "logical",
        is_disjoint = "logical",
        alpha = "numeric",
        r = "numeric",
        elm2clu_dict = "list",
        clu2elm_dict = "list",
        affinity_matrix = "matrix"
    )
)

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
#' page-rank calculation:
#' * "prpack": use PPR algorithms in igraph.
#' * "power_iteration": use power_iteration method.
#' @param row_normalize Whether to normalize all rows in clustering_result
#' so they sum to one before calculating ECS. It is recommended to set this to
#' TRUE, which will lead to slightly different ECS values compared to clusim.
#' @param r A numeric hierarchical scaling parameter.
#' @param rescale_path_type A string; rescale the hierarchical height by:
#' * "max" : the maximum path from the root.
#' * "min" : the minimum path form the root.
#' * "linkage" : use the linkage distances in the clustering.
#' @param dist_rescaled A logical: if TRUE, the linkage distances are linearly
#' rescaled to be in-between 0 and 1.
#' @param ... This argument is not used.
#'
#' @return A Clustering object.
#'
#' @keywords internal
setGeneric(
    "create_clustering",
    function(clustering_result, ...) {
        standardGeneric("create_clustering")
    }
)

#' @describeIn create_clustering Create Clustering Object from Numeric Vector
setMethod(
    "create_clustering",
    signature(clustering_result = "numeric"),
    function(clustering_result,
             alpha = 0.9) {
        # convert to character
        element.names <- names(clustering_result)
        clustering_result <- as.character(clustering_result)
        names(clustering_result) <- element.names

        # create Clustering object in separate function
        return(create_flat_disjoint_clustering(
            clustering_result,
            alpha
        ))
    }
)

#' @describeIn create_clustering Create Clustering Object from Integer Vector
setMethod(
    "create_clustering",
    signature(clustering_result = "integer"),
    function(clustering_result,
             alpha = 0.9) {
        # convert to character
        element.names <- names(clustering_result)
        clustering_result <- as.character(clustering_result)
        names(clustering_result) <- element.names

        # create Clustering object in separate function
        return(create_flat_disjoint_clustering(
            clustering_result,
            alpha
        ))
    }
)

#' @describeIn create_clustering Create Clustering Object from Character Vector
setMethod(
    "create_clustering",
    signature(clustering_result = "character"),
    function(clustering_result,
             alpha = 0.9) {
        # create Clustering object in separate function
        return(create_flat_disjoint_clustering(
            clustering_result,
            alpha
        ))
    }
)

#' @describeIn create_clustering Create Clustering Object from Factor Vector
setMethod(
    "create_clustering",
    signature(clustering_result = "factor"),
    function(clustering_result,
             alpha = 0.9) {
        # convert to character
        element.names <- names(clustering_result)
        clustering_result <- as.character(clustering_result)
        names(clustering_result) <- element.names

        # create Clustering object in separate function
        return(create_flat_disjoint_clustering(
            clustering_result,
            alpha
        ))
    }
)

create_flat_disjoint_clustering <- function(clustering_result,
                                            alpha) {
    # disjoint partitions
    n_elements <- length(clustering_result)

    # names of elements
    if (is.null(names(clustering_result))) {
        element.names <- as.character(1:n_elements)
    } else {
        element.names <- names(clustering_result)
    }

    # create dictionaries, or lists
    clu2elm_dict <- create_clu2elm_dict(clustering_result)
    elm2clu_dict <- list()

    # create object
    clustering <- methods::new("Clustering",
        n_elements = n_elements,
        is_hierarchical = FALSE,
        is_disjoint = TRUE,
        names = element.names,
        alpha = alpha,
        r = 0,
        elm2clu_dict = elm2clu_dict,
        clu2elm_dict = clu2elm_dict,
        affinity_matrix = matrix()
    )

    # calculate affinity matrix
    clustering@affinity_matrix <- ppr_partition(
        clustering,
        alpha
    )
    return(clustering)
}

#' @describeIn create_clustering Create Clustering Object from base matrix
setMethod(
    "create_clustering",
    signature(clustering_result = "matrix"),
    function(clustering_result,
             alpha = 0.9,
             ppr_implementation = "prpack",
             row_normalize = TRUE) {
        return(create_flat_overlapping_clustering(
            clustering_result,
            alpha,
            ppr_implementation,
            row_normalize
        ))
    }
)

#' @describeIn create_clustering Create Clustering Object from Matrix::Matrix
setMethod(
    "create_clustering",
    signature(clustering_result = "Matrix"),
    function(clustering_result,
             alpha = 0.9,
             ppr_implementation = "prpack",
             row_normalize = TRUE) {
        # convert clustering_result to base matrix
        clustering_result <- as.matrix(clustering_result)
        return(create_flat_overlapping_clustering(
            clustering_result,
            alpha,
            ppr_implementation,
            row_normalize
        ))
    }
)

create_flat_overlapping_clustering <- function(clustering_result,
                                               alpha,
                                               ppr_implementation,
                                               row_normalize) {
    # print(clustering_result)
    # soft clustering
    n_elements <- nrow(clustering_result)

    # check that entries are non-negative
    if (any(clustering_result < 0)) {
        stop("clustering_result should only contain non-negative entries.")
    }
    # check that each element belongs to at least one cluster
    if (any(rowSums(clustering_result) == 0)) {
        stop("Every element of clustering_result must be in at least one cluster.")
    }

    if (row_normalize) {
        # normalize clustering_result so each row sums to one
        clustering_result <- clustering_result / rowSums(clustering_result)
    }

    # names of elements
    if (is.null(rownames(clustering_result))) {
        element.names <- as.character(1:n_elements)
    } else {
        element.names <- rownames(clustering_result)
    }

    # create dictionaries, or lists
    clu2elm_dict <- list()
    elm2clu_dict <- create_elm2clu_dict_overlapping(clustering_result)

    # create object
    clustering <- methods::new("Clustering",
        n_elements = n_elements,
        is_hierarchical = FALSE,
        is_disjoint = FALSE,
        names = element.names,
        alpha = alpha,
        r = 0,
        elm2clu_dict = elm2clu_dict,
        clu2elm_dict = clu2elm_dict,
        affinity_matrix = matrix()
    )

    # create affinity matrix
    cielg <- make_cielg_overlapping(clustering_result)
    clustering@affinity_matrix <- numerical_ppr_scores(cielg,
        clustering,
        ppr_implementation =
            ppr_implementation
    )
    return(clustering)
}

setOldClass("hclust")
#' @describeIn create_clustering Create Clustering Object from hclust
setMethod(
    "create_clustering",
    signature(clustering_result = "hclust"),
    function(clustering_result,
             alpha = 0.9,
             r = 1,
             rescale_path_type = "max",
             ppr_implementation = "prpack",
             dist_rescaled = FALSE) {
        # hierarchical partitions
        n_elements <- length(clustering_result$order)

        # names of elements
        if (is.null(clustering_result$labels)) {
            element.names <- as.character(1:n_elements)
        } else {
            element.names <- clustering_result$labels
        }

        # create dictionaries, or lists
        clu2elm_dict <- create_clu2elm_dict_hierarchical(clustering_result)
        elm2clu_dict <- create_elm2clu_dict_hierarchical(clustering_result)

        # create object
        clustering <- methods::new("Clustering",
            n_elements = n_elements,
            is_hierarchical = TRUE,
            is_disjoint = TRUE,
            names = element.names,
            alpha = alpha,
            r = r,
            elm2clu_dict = elm2clu_dict,
            clu2elm_dict = clu2elm_dict,
            affinity_matrix = matrix()
        )

        # create affinity matrix
        cielg <- make_cielg_hierarchical(
            clustering_result,
            clu2elm_dict,
            r,
            rescale_path_type,
            dist_rescaled
        )
        clustering@affinity_matrix <- numerical_ppr_scores(cielg,
            clustering,
            ppr_implementation =
                ppr_implementation
        )
        return(clustering)
    }
)

setGeneric("length")
#' Length of an Object
#' @description Get the number of elements in the Clustering.
#'
#' @param x The Clustering object.
#'
#' @return The number of elements.
#'
#' @keywords internal
setMethod(
    "length",
    signature(x = "Clustering"),
    function(x) {
        return(x@n_elements)
    }
)

setMethod(
    "show",
    signature(object = "Clustering"),
    function(object) {
        hierarchical <- if (object@is_hierarchical) "hierarchical" else "flat"
        overlapping <- if (object@is_disjoint) "disjoint" else "overlapping"
        cat(paste0(
            "A ", hierarchical, ", ", overlapping, " clustering of ",
            object@n_elements, " elements.\n"
        ))
    }
)

setGeneric("print")
#' Print an Object
#' @description Prints out information about the Clustering, including
#' number of elements.
#'
#' @param x The Clustering object.
#'
#' @return The printed character string.
#'
#' @keywords internal
setMethod(
    "print",
    signature(x = "Clustering"),
    function(x) {
        hierarchical <- if (x@is_hierarchical) "hierarchical" else "flat"
        overlapping <- if (x@is_disjoint) "disjoint" else "overlapping"
        print(paste0(
            "A ", hierarchical, ", ", overlapping, " clustering of ",
            x@n_elements, " elements."
        ))
    }
)
