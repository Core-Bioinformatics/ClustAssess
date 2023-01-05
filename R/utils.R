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
