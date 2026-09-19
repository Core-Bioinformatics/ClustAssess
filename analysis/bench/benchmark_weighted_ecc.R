#!/usr/bin/env Rscript

env_int <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) default else as.integer(value)
}

env_num <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) default else as.numeric(value)
}

env_flag <- function(name, default = FALSE) {
    value <- tolower(Sys.getenv(name, unset = if (default) "true" else "false"))
    value %in% c("1", "true", "yes", "on")
}

lib_path <- Sys.getenv("CA_ECC_LIB", unset = "")
if (nzchar(lib_path)) {
    .libPaths(c(lib_path, .libPaths()))
}

suppressPackageStartupMessages(library(ClustAssess))

weighted_ecc_reference <- function(clustering_list,
                                   weights = NULL,
                                   calculate_sim_matrix = FALSE) {
    n_clusterings <- length(clustering_list)
    if (n_clusterings <= 1) {
        stop("At least two clusterings are required to calculate the consistency.")
    }
    if (is.null(weights)) {
        weights <- rep(1, n_clusterings)
    }

    consistency <- 0
    total_weights <- sum(weights)
    norm_factor <- total_weights * (total_weights - 1) / 2
    if (calculate_sim_matrix) {
        ecs_sim_matrix <- matrix(0, nrow = n_clusterings, ncol = n_clusterings)
    }

    for (i in seq_len(n_clusterings - 1)) {
        consistency <- consistency +
            weights[i] / norm_factor * (weights[i] - 1) / 2
        for (j in seq(from = i + 1, to = n_clusterings)) {
            current_ecs <- element_sim_elscore(
                clustering_list[[i]],
                clustering_list[[j]]
            )
            consistency <- consistency +
                current_ecs * weights[i] / norm_factor * weights[j]
            if (calculate_sim_matrix) {
                ecs_sim_matrix[i, j] <- mean(current_ecs)
                ecs_sim_matrix[j, i] <- ecs_sim_matrix[i, j]
            }
        }
        if (calculate_sim_matrix) {
            ecs_sim_matrix[i, i] <- 1
        }
    }

    consistency <- consistency +
        weights[n_clusterings] / norm_factor *
        (weights[n_clusterings] - 1) / 2

    if (!calculate_sim_matrix) {
        return(consistency)
    }
    ecs_sim_matrix[n_clusterings, n_clusterings] <- 1
    list(ecc = consistency, ecs_matrix = ecs_sim_matrix)
}

make_synthetic_input <- function() {
    n_elements <- env_int("CA_ECC_ELEMENTS", 2000L)
    n_clusterings <- env_int("CA_ECC_CLUSTERINGS", 20L)
    n_clusters <- env_int("CA_ECC_CLUSTERS", 20L)
    perturb_fraction <- env_num("CA_ECC_PERTURB", 0.10)
    seed <- env_int("CA_ECC_SEED", 1L)

    set.seed(seed)
    base <- rep(seq_len(n_clusters), length.out = n_elements)
    element_names <- paste0("cell_", seq_len(n_elements))
    n_perturb <- min(n_elements, max(0L, round(n_elements * perturb_fraction)))
    clustering_list <- lapply(seq_len(n_clusterings), function(i) {
        membership <- base
        if (n_perturb > 0L) {
            changed <- sample.int(n_elements, n_perturb)
            membership[changed] <- sample.int(n_clusters, n_perturb, replace = TRUE)
        }
        setNames(membership, element_names)
    })

    max_weight <- env_int("CA_ECC_MAX_WEIGHT", 3L)
    weights <- if (max_weight <= 1L) {
        rep(1, n_clusterings)
    } else {
        sample.int(max_weight, n_clusterings, replace = TRUE)
    }

    list(
        clustering_list = clustering_list,
        weights = weights,
        calculate_sim_matrix = env_flag("CA_ECC_SIM_MATRIX", FALSE),
        source = "synthetic"
    )
}

load_input <- function() {
    input_file <- Sys.getenv("CA_ECC_INPUT", unset = "")
    if (!nzchar(input_file)) {
        return(make_synthetic_input())
    }
    input <- readRDS(input_file)
    required <- c("clustering_list", "weights")
    if (!is.list(input) || !all(required %in% names(input))) {
        stop("CA_ECC_INPUT must contain clustering_list and weights.")
    }
    input$calculate_sim_matrix <- input$calculate_sim_matrix %||% FALSE
    input$source <- normalizePath(input_file)
    input
}

`%||%` <- function(x, y) if (is.null(x)) y else x

summarize_allocations <- function(fun, threshold) {
    profile_file <- tempfile("weighted-ecc-", fileext = ".Rprofmem")
    on.exit(unlink(profile_file), add = TRUE)
    utils::Rprofmem(profile_file, threshold = threshold)
    on.exit(utils::Rprofmem(NULL), add = TRUE)
    invisible(fun())
    utils::Rprofmem(NULL)

    lines <- readLines(profile_file, warn = FALSE)
    sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    sizes <- sizes[!is.na(sizes)]
    c(
        allocation_count = length(sizes),
        allocated_bytes = sum(sizes),
        largest_allocation_bytes = if (length(sizes)) max(sizes) else 0
    )
}

measure <- function(name,
                    fun,
                    iterations,
                    calls_per_iteration,
                    warmups,
                    allocation_threshold) {
    run_calls <- function() {
        for (call_index in seq_len(calls_per_iteration)) {
            invisible(fun())
        }
    }

    if (warmups > 0L) {
        for (i in seq_len(warmups)) {
            run_calls()
        }
    }

    elapsed <- numeric(iterations)
    for (i in seq_len(iterations)) {
        gc(FALSE)
        timing <- system.time(run_calls())
        elapsed[i] <- unname(timing[["elapsed"]]) / calls_per_iteration
    }
    allocations <- summarize_allocations(fun, allocation_threshold)

    data.frame(
        method = name,
        iterations = iterations,
        calls_per_iteration = calls_per_iteration,
        min_seconds = min(elapsed),
        median_seconds = stats::median(elapsed),
        max_seconds = max(elapsed),
        allocation_count = unname(allocations[["allocation_count"]]),
        allocated_mib = unname(allocations[["allocated_bytes"]]) / 1024^2,
        largest_allocation_mib =
            unname(allocations[["largest_allocation_bytes"]]) / 1024^2,
        stringsAsFactors = FALSE
    )
}

input <- load_input()
clustering_list <- input$clustering_list
weights <- input$weights
calculate_sim_matrix <- isTRUE(input$calculate_sim_matrix)
installed_weighted_ecc <- getFromNamespace(
    "weighted_element_consistency",
    "ClustAssess"
)

reference_result <- weighted_ecc_reference(
    clustering_list,
    weights,
    calculate_sim_matrix
)
installed_result <- installed_weighted_ecc(
    clustering_list,
    weights,
    calculate_sim_matrix
)
parity <- isTRUE(all.equal(
    installed_result,
    reference_result,
    tolerance = 1e-12,
    check.attributes = TRUE
))
if (!parity) {
    stop("Installed and reference weighted ECC results differ.")
}

iterations <- env_int("CA_ECC_ITERATIONS", 5L)
calls_per_iteration <- env_int("CA_ECC_CALLS_PER_ITERATION", 1L)
warmups <- env_int("CA_ECC_WARMUPS", 1L)
allocation_threshold <- env_num("CA_ECC_RPROFMEM_THRESHOLD", 0)
if (iterations < 1L || calls_per_iteration < 1L || warmups < 0L) {
    stop("Iterations and calls must be positive; warmups cannot be negative.")
}

methods <- list(
    installed = function() installed_weighted_ecc(
        clustering_list,
        weights,
        calculate_sim_matrix
    ),
    pairwise_reference = function() weighted_ecc_reference(
        clustering_list,
        weights,
        calculate_sim_matrix
    )
)

namespace <- asNamespace("ClustAssess")
if (exists("encode_flat_disjoint_memberships", envir = namespace, inherits = FALSE) &&
    exists("weightedDisjointECS", envir = namespace, inherits = FALSE)) {
    encode <- get("encode_flat_disjoint_memberships", envir = namespace)
    weighted_cpp <- get("weightedDisjointECS", envir = namespace)
    encoded <- encode(clustering_list)
    methods$encode_only <- function() encode(clustering_list)
    methods$cpp_preencoded <- function() weighted_cpp(
        encoded,
        weights,
        calculate_sim_matrix
    )
}

results <- do.call(rbind, lapply(names(methods), function(name) {
    measure(
        name,
        methods[[name]],
        iterations,
        calls_per_iteration,
        warmups,
        allocation_threshold
    )
}))

cat("source=", input$source, "\n", sep = "")
cat(
    "package_version=", as.character(utils::packageVersion("ClustAssess")),
    " elements=", length(clustering_list[[1]]),
    " clusterings=", length(clustering_list),
    " sim_matrix=", calculate_sim_matrix,
    " parity=", parity,
    "\n",
    sep = ""
)
print(results, row.names = FALSE, digits = 5)
