#' @useDynLib ClustAssess
#' @importFrom Rcpp sourceCpp
NULL

calculate_dist = function(x, dist_method, p=2){
  if (dist_method %in% c("euclidean", "maximum", "manhattan", "canberra",
                         "binary", "minkowski")){
    d <- stats::dist(x, method=dist_method, p=p)
  } else if (dist_method=='pearson'){
    d <- 1 - stats::cor(x, method='pearson')
  } else {
    stop('Distance measure not recognized. Please choose one of euclidean,
         maximum, canberra, binary, minkowski and pearson')
  }

  return(d)
}

#' Consensus Clustering and Proportion of Ambiguously Clustered Pairs
#' @description Calculate consensus clustering and proportion of ambiguously
#' clustered pairs (PAC) with hierarchical clustering.
#'
#' @param x A samples x features normalized data matrix.
#' @param k_min The minimum number of clusters calculated.
#' @param k_max The maximum number of clusters calculated.
#' @param n_reps The total number of subsamplings and reclusterings of the data;
#' this value needs to be high enough to ensure PAC converges; convergence can
#' be assessed with pac_convergence.
#' @param p_sample The proportion of samples included in each subsample.
#' @param p_feature The proportion of features included in each subsample.
#' @param dist_method The distance measure for the distance matrix used in
#' hclust; must be one of "euclidean", "maximum", "manhattan", "canberra",
#' "binary" or "minkowski".
#' @param linkage The linkage method used in hclust; must be one of "ward.D",
#' "ward.D2", "single", "complete", "average", "mcquitty", "median" or
#' "centroid"
#' @param upper_lim The upper limit for determining whether a pair is
#' clustered ambiguously; the higher this value, the higher the PAC.
#' @param lower_lim The lower limit for determining whether a pair is
#' clustered ambiguously; the lower this value, the higher the PAC.
#' @param p_minkowski The power of the Minkowski distance.
#' @param verbose Logical value used for choosing to display a progress bar or not.
#'
#' @return A data.frame with PAC values across iterations, as well as parameter
#' values used when calling the method.
#'
#' @references Monti, S., Tamayo, P., Mesirov, J., & Golub, T. (2003).
#' Consensus clustering: a resampling-based method for class discovery and
#' visualization of gene expression microarray data. Machine learning, 52(1),
#' 91-118. https://doi.org/10.1023/A:1023949509487
#' @references Senbabaoglu, Y., Michailidis, G., & Li, J. Z. (2014).
#' Critical limitations of consensus clustering in class discovery.
#' Scientific reports, 4(1), 1-13. https://doi.org/10.1038/srep06207
#'
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' pac.res = consensus_cluster(iris[,1:4], k_max=20)
#' pac_convergence(pac.res, k_plot=c(3,5,7,9))
consensus_cluster <- function(x, k_min=3, k_max=100, n_reps=100, p_sample=0.8,
                              p_feature=1.0, p_minkowski=2,
                              dist_method='euclidean', linkage='complete',
                              lower_lim=0.1, upper_lim=0.9,
                              verbose=TRUE){
  n.samples = floor(p_sample * nrow(x))
  n.features = floor(p_feature * ncol(x))
  # threshold for k_max so we don't accidentally look for too many clusters
  if (k_max > n.samples){
    k_max = n.samples
  }

  indicator = matrix(0, nrow = nrow(x), ncol = nrow(x))

  # initialize connectivity matrices
  connectivity = list()
  for (i in k_min:k_max){
    connectivity[[i]] = matrix(0, nrow = nrow(x), ncol = nrow(x))
  }

  pac.matrix = matrix(0, nrow=(k_max-k_min+1), ncol=n_reps)

  if(verbose) {
    message("Calculating consensus clustering")
    pb <- progress::progress_bar$new(
      format = "[:bar] :current/:total eta: :eta  total elapsed:  :elapsed",
      total = n_reps,
      clear = FALSE,
      width = 80
    )

    pb$tick(0)
  }

  for (i in 1:n_reps){
    if(verbose)
      pb$tick()

    sample.indices = sample(x=nrow(x), size=n.samples, replace=FALSE)
    indicator[sample.indices, sample.indices] =
      indicator[sample.indices, sample.indices] + 1

    feature.indices = sample(x=ncol(x), size=n.features, replace=FALSE)
    d <- calculate_dist(x[sample.indices, feature.indices], dist_method,
                        p_minkowski)
    tree <- fastcluster::hclust(d, method=linkage)

    for (k in k_min:k_max){
      res <- stats::cutree(tree, k=k)

      connectivity[[k]] = update_connectivity_cpp(connectivity[[k]],
                                                   sample.indices,
                                                   res)

      pac.value = calculate_pac_cpp(indicator, connectivity[[k]], lower_lim,
                                    upper_lim)
      pac.matrix[(k-k_min+1), i] = pac.value
    }
  }

  convergence = data.frame(dist_method = dist_method,
                           linkage = linkage,
                           iteration = rep(1:n_reps, each=(k_max-k_min+1)),
                           pac = as.vector(pac.matrix),
                           lower_lim = lower_lim,
                           upper_lim = upper_lim,
                           n.clusters = rep(k_min:k_max, times=n_reps),
                           p_sample = p_sample,
                           p_feature = p_feature)

  return(convergence)
}


#' PAC Convergence Plot
#' @description Plot PAC across iterations for a set of k to assess convergence.
#'
#' @param pac_res The data.frame output by consensus_cluster.
#' @param k_plot A vector with values of k to plot.
#'
#' @return A ggplot2 object with the convergence plot.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' pac.res = consensus_cluster(iris[,1:4], k_max=20)
#' pac_convergence(pac.res, k_plot=c(3,5,7,9))
pac_convergence = function(pac_res, k_plot){
  conv = pac_res %>% dplyr::filter(.data$n.clusters %in% k_plot)
  conv$n.clusters = as.factor(conv$n.clusters)

  ggplot2::ggplot(data=conv,
                  ggplot2::aes(x=.data$iteration, y=.data$pac,
                               color=.data$n.clusters,
                               group=.data$n.clusters)) +
          ggplot2::geom_line() +
          ggplot2::labs(title='PAC convergence')
}


#' PAC Landscape Plot
#' @description Plot final PAC values across range of k to find optimal number
#' of clusters.
#'
#' @param pac_res The data.frame output by consensus_cluster.
#' @param n_shade The number of iterations to shade to show the
#' variability of PAC across the last n_shade iterations.
#'
#' @return A ggplot2 object with the final PAC vs k plot.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' pac.res = consensus_cluster(iris[,1:4], k_max=20)
#' pac_landscape(pac.res)
pac_landscape = function(pac_res, n_shade = max(pac_res$iteration)/5){
  # threshold n_shade
  if (n_shade > max(pac_res$iteration)){
    n_shade = max(pac_res$iteration)
  }

  # create data frame with shaded regions
  pac.df = pac_res %>% dplyr::group_by(.data$n.clusters) %>%
    dplyr::top_n(n_shade, wt=.data$iteration) %>%
    dplyr::mutate(pmin=min(.data$pac), pmax=max(.data$pac)) %>%
    dplyr::top_n(1, wt=.data$iteration)

  ggplot2::ggplot(data=pac.df,
                  ggplot2::aes(x=.data$n.clusters, y=.data$pac)) +
         ggplot2::geom_line(color='blue') +
         ggplot2::geom_ribbon(ggplot2::aes(ymin=pmin, ymax=pmax), alpha=0.2) +
         ggplot2::scale_color_gradient() +
         ggplot2::labs(title='Final PAC values')
}
