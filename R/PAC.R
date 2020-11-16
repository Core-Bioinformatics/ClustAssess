#' @useDynLib ClustAssess
#' @importFrom Rcpp sourceCpp
NULL


consensus.cdf <- function(consensus, cdf.length = 1e3){
  consensus.order = consensus[order(consensus)]
  n.elements = length(consensus)
  consensus.index = seq(0, 1, length = min(n.elements, cdf.length))
  cdf = stats::ecdf(consensus.order)(consensus.index)
  cdf[1] = 0
  return(list(cdf = cdf, consensus.index = consensus.index))
}

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
                              calculate.cdf=FALSE, cdf.length = 1e3,
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
  # consensus = list()

  # initialize connectivity and consensus matrices
  for (i in k.min:k.max){
    connectivity[[i]] = matrix(0, nrow = nrow(x), ncol = nrow(x))
    # consensus[[i]] = matrix(0, nrow = nrow(x), ncol = nrow(x))
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

      #consensus[[k]] = connectivity[[k]] / indicator
      #consensus[[k]][indicator == 0] = 0
      #consensus.upper = consensus[[k]][!diag(nrow(x))]
      #n.elements = length(consensus.upper)

      # kinda slow
      # m = mean(lower.lim < consensus.upper & consensus.upper < upper.lim)

      # faster
      #lower.mask = (lower.lim < consensus.upper)
      #upper.and.lower.mask = (consensus.upper[lower.mask] < upper.lim)
      #pac.value = sum(upper.and.lower.mask) / n.elements

      pac.value = calculate_pac_cpp(indicator, connectivity[[k]], lower.lim, upper.lim)
      pac.matrix[(k-k.min+1), i] = pac.value

      # convergence = rbind(convergence, data.frame(dist.method = dist.method,
      #                                             linkage = linkage,
      #                                             iteration = i,
      #                                             pac = pac.value,
      #                                             lower.lim = lower.lim,
      #                                             upper.lim = upper.lim,
      #                                             n.clusters = k,
      #                                             p.item = p.item))

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

  # calculate CDF
  if (calculate.cdf){
    cdf.length = min(nrow(x)**2, cdf.length)
    index.vec = rep(0, (k.max-k.min+1)*cdf.length)
    cdf.vec = rep(0, (k.max-k.min+1)*cdf.length)
    for (k in k.min:k.max){
      res.cdf = consensus.cdf(consensus[[k]], cdf.length)
      index.vec[((k-k.min)*cdf.length+1):((k-k.min+1)*cdf.length)] = res.cdf$consensus.index
      cdf.vec[((k-k.min)*cdf.length+1):((k-k.min+1)*cdf.length)] = res.cdf$cdf
      # cdf = rbind(cdf, data.frame(consensus.index = res.cdf$consensus.index,
      #                             CDF = res.cdf$cdf,
      #                             n.clusters = k))
    }
    cdf = data.frame(consensus.index = cdf.vec,
                     CDF = cdf.vec,
                     n.clusters = rep(k.min:k.max, each=cdf.length))
  } else{
    cdf = NULL
  }

  pac.res = list(convergence=convergence, cdf=cdf)
  return(pac.res)
}

update.connectivity.naive = function(connectivity,
                                     sampling.indices,
                                     cluster.assignments){
  for (nk in unique(cluster.assignments)){
    same.cluster.indices = (cluster.assignments == nk)
    connectivity[sampling.indices, sampling.indices][same.cluster.indices, same.cluster.indices] =
      connectivity[sampling.indices, sampling.indices][same.cluster.indices, same.cluster.indices] + 1
  }
  return(connectivity)
}

update.connectivity.loop = function(connectivity,
                                     sampling.indices,
                                     cluster.assignments){
  n.samples = length(sampling.indices)
  for (i.cell in 1:n.samples){
    cell.1 = sampling.indices[i.cell]
    for (j.cell in (i.cell+1):n.samples){
      if (j.cell>n.samples){
        break
      }
      cell.2 = sampling.indices[j.cell]
      if (cluster.assignments[i.cell] == cluster.assignments[j.cell]){
        connectivity[cell.1, cell.2] = connectivity[cell.1, cell.2] + 1
        connectivity[cell.2, cell.1] = connectivity[cell.1, cell.2]
      }
    }
  }
  return(connectivity)
}

update.connectivity.ccp = function(connectivity,
                                     sampling.indices,
                                     cluster.assignments){
  names( cluster.assignments ) <- sampling.indices
  cls <- lapply( unique( cluster.assignments ), function(i) as.numeric( names( cluster.assignments[ cluster.assignments %in% i ] ) ) )  #list samples by clusterId
  nelts <- 1:ncol( connectivity )
  for ( i in 1:length( cls ) ) {
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    connectivity <- connectivity + updt
  }
  return(connectivity)
}

update.connectivity.sort = function(connectivity,
                                   sampling.indices,
                                   cluster.assignments){
  names(cluster.assignments) <- sampling.indices
  n.samples = length(cluster.assignments)
  n.items = nrow(connectivity)
  sorted.assignments = sort(cluster.assignments)
  for (i in 1:sorted.assignments[n.samples]){
    update.items = rep(0, n.items)

    indices = findInterval(c(i-.5, i+.5), sorted.assignments)
    update.indices = as.integer(names(sorted.assignments)[(indices[1]+1):indices[2]])
    update.items[update.indices]=1

    updt <- outer(update.items, update.items) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    connectivity <- connectivity + updt
  }
  return(connectivity)
}

# .onUnload <- function (libpath) {
#   library.dynam.unload("ClustAssess", libpath)
# }


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
#' @param k.plot
#'
#' @return
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
pac.cdf = function(pac.res, k.plot){
  cdf = pac.res$cdf %>% dplyr::filter(.data$n.clusters %in% k.plot)

  ggplot2::ggplot(cdf,
                  ggplot2::aes(x=.data$consensus.index, y=.data$CDF,
                               color=.data$n.clusters,
                               group=.data$n.clusters)) +
          ggplot2::geom_step() +
          ggplot2::scale_color_gradient(low = "#cde6fe", high = "#012345") +
          ggplot2::labs(title='PAC CDF')
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
