# generate values from an interval that simulates the logarithmic scale
generate_breaks <- function(min.range, max.range) {
    start.point <- min.range + 5 - min.range %% 5
    breaks.list <- list()
    index <- 1
    step <- 5
    while (start.point < max.range) {
        breaks.list[[index]] <- start.point
        index <- index + 1

        start.point <- start.point + step

        # the step is doubled after every two iterations
        if (index %% 2) {
            step <- step * 2
        }
    }

    c(min.range, unlist(breaks.list), max.range)
}

# source: https://github.com/satijalab/seurat/blob/763259d05991d40721dee99c9919ec6d4491d15e/R/objects.R#L2415
get_count_features <- function(counts_matrix) {
    return(list(
        nCount = Matrix::colSums(counts_matrix),
        nFeature = Matrix::colSums(counts_matrix > 0)
    ))
}
