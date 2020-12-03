#' @useDynLib ClustAssess
#' @importFrom Rcpp sourceCpp
NULL

calculate.dist = function(x, dist.method){
  if (dist.method %in% c("euclidean", "maximum", "manhattan", "canberra",
                         "binary", "minkowski")){
    d <- stats::dist(x, method=dist.method)
  } else if (dist.method == 'pearson'){
    d <- 1 - stats::cor(x, method='pearson')
  } else {
    stop('Distance measure not recognized. Please choose one of euclidean,
         maximum, canberra, binary, minkowski and pearson')
  }

  return(d)
}

#' Title
#'
#' @param x
#' @param k.min
#' @param k.max
#' @param n.reps
#' @param p.item
#' @param p.feature
#' @param dist.method
#' @param linkage
#' @param upper.lim
#' @param lower.lim
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
consensus.cluster <- function(x, k.min=3, k.max=100, n.reps=100, p.item=0.8,
                              p.feature=1.0,
                              dist.method='euclidean', linkage='complete',
                              upper.lim=0.9, lower.lim=0.1,
                              calculate.dist.upfront=FALSE){
  n.items = floor(p.item * nrow(x))
  # threshold for k.max so we don't accidentally look for too many clusters
  if (k.max > n.items){
    k.max = n.items
  }
  indices = 1:nrow(x)

  if (calculate.dist.upfront){
    all.dist <- calculate.dist(x, dist.method)
    all.dist <- as.matrix(all.dist)
  }

  indicator = matrix(0, nrow = nrow(x), ncol = nrow(x))
  connectivity = list()

  # initialize connectivity matrices
  for (i in k.min:k.max){
    connectivity[[i]] = matrix(0, nrow = nrow(x), ncol = nrow(x))
  }

  convergence = NULL
  pac.matrix = matrix(0, nrow=(k.max-k.min+1), ncol=n.reps)

  for (i in 1:n.reps){
    item.indices = sample(x=indices, size=n.items, replace=FALSE)
    indicator[item.indices, item.indices] =
      indicator[item.indices, item.indices] + 1

    # for hclust
    if (calculate.dist.upfront){
      d <- all.dist[item.indices, item.indices]
      d <- stats::as.dist(d)
    } else {
      d <- calculate.dist(x[item.indices,], dist.method)
    }

    tree <- fastcluster::hclust(d, method=linkage)

    for (k in k.min:k.max){
      res <- stats::cutree(tree, k=k)
      #connectivity[[k]] = update.connectivity.sort(connectivity[[k]],
      #                                              item.indices,
      #                                              res)

      connectivity[[k]] = update_connectivity_cpp(connectivity[[k]],
                                                   item.indices,
                                                   res)

      pac.value = calculate_pac_cpp(indicator, connectivity[[k]], lower.lim, upper.lim)
      pac.matrix[(k-k.min+1), i] = pac.value
    }
  }

  convergence = data.frame(dist.method = dist.method,
                           linkage = linkage,
                           iteration = rep(1:n.reps, each=(k.max-k.min+1)),
                           pac = as.vector(pac.matrix),
                           lower.lim = lower.lim,
                           upper.lim = upper.lim,
                           n.clusters = rep(k.min:k.max, times=n.reps),
                           p.item = p.item)

  pac.res = list(convergence=convergence)
  return(pac.res)
}


#' Title
#'
#' @param pac.res
#' @param k.plot
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
pac.convergence = function(pac.res, k.plot){
  conv = pac.res$convergence %>% dplyr::filter(.data$n.clusters %in% k.plot)

  ggplot2::ggplot(data=conv,
                  ggplot2::aes(x = .data$iteration, y = .data$pac,
                               color = .data$n.clusters,
                               group=.data$n.clusters)) +
          ggplot2::geom_line() +
          ggplot2::scale_color_gradient(low = "#cde6fe", high = "#012345") +
          ggplot2::labs(title='PAC convergence')
}


#' Title
#'
#' @param pac.res
#' @param n.shade
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
pac.landscape = function(pac.res, n.shade = max(pac.res$convergence$iteration)/5){
  # threshold n.shade
  if (n.shade > max(pac.res$convergence$iteration)){
    n.shade = max(pac.res$convergence$iteration)
  }

  # create data frame with shaded regions
  pac.df = pac.res$convergence %>% dplyr::group_by(.data$n.clusters) %>%
    dplyr::top_n(n.shade, wt=.data$iteration) %>%
    dplyr::mutate(pmin = min(.data$pac), pmax=max(.data$pac)) %>%
    dplyr::top_n(1, wt=.data$iteration)

  ggplot2::ggplot(data=pac.df,
                  ggplot2::aes(x=.data$n.clusters, y=.data$pac)) +
         ggplot2::geom_line(color='blue') +
         ggplot2::geom_ribbon(ggplot2::aes(ymin=pmin, ymax=pmax), alpha=0.2) +
         ggplot2::scale_color_gradient() +
         ggplot2::labs(title='Final PAC values')
}
