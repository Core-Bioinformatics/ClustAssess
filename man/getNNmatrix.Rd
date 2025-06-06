% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getNNmatrix}
\alias{getNNmatrix}
\title{Computes the NN adjacency matrix given the neighbours}
\usage{
getNNmatrix(nnRanked, k = -1L, start = 0L, prune = 0)
}
\arguments{
\item{nnRanked}{A matrix with the lists of the nearest neighbours for each point}

\item{k}{The number of neighbours to consider. Defaults to \code{-1}, which means all neighbours.}

\item{start}{The index of the first neighbour to consider. Defaults to \code{0}.}

\item{prune}{The threshold to prune the SNN matrix. If -1, the function will only return the NN matrix. Defaults to \code{0}.}
}
\value{
A list with the NN and SNN adjacency matrices.
}
\description{
Computes the NN adjacency matrix given the neighbours
}
