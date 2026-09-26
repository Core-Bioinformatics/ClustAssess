#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grepl(file_arg, args)][1])
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package `pkgload` is required to run this benchmark.")
}

pkgload::load_all(repo_root, quiet = TRUE)

profile_enabled <- tolower(Sys.getenv("CA_BENCH_PROFILE", "false")) %in%
    c("1", "true", "yes")
profile_dir <- Sys.getenv(
    "CA_BENCH_PROFILE_DIR",
    file.path(repo_root, "analysis", "bench", "profiles")
)
profile_interval <- as.numeric(Sys.getenv("CA_BENCH_PROFILE_INTERVAL", "0.01"))
profile_memory_threshold <- as.numeric(Sys.getenv("CA_BENCH_PROFILE_MEMORY_THRESHOLD", "0"))

env_flag <- function(name, default = FALSE) {
    value <- Sys.getenv(name, if (default) "true" else "false")
    tolower(value) %in% c("1", "true", "yes")
}

env_int <- function(name, default) {
    as.integer(Sys.getenv(name, as.character(default)))
}

safe_filename <- function(label) {
    gsub("[^A-Za-z0-9_.-]+", "_", label)
}

summarize_rprofmem <- function(memory_file) {
    lines <- readLines(memory_file, warn = FALSE)
    allocation_lines <- grepl("^[0-9]+[[:space:]]*:", lines)
    allocation_bytes <- suppressWarnings(as.numeric(
        sub("^([0-9]+)[[:space:]]*:.*", "\\1", lines[allocation_lines])
    ))

    data.frame(
        allocation_count = length(allocation_bytes),
        allocated_bytes = sum(allocation_bytes, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}

make_exact_pca_processor <- function(default_npcs) {
    make_pca_matrix_processing(actual_npcs = default_npcs, method = "exact")
}

make_irlba_pca_processor <- function(default_npcs) {
    make_pca_matrix_processing(actual_npcs = default_npcs, method = "irlba")
}

make_low_rank_expression <- function(n_genes, n_cells, n_modules = 6) {
    n_modules <- min(n_modules, n_genes, n_cells)
    gene_module <- rep(seq_len(n_modules), length.out = n_genes)
    cell_module <- rep(seq_len(n_modules), length.out = n_cells)
    signal <- outer(gene_module, cell_module, FUN = "==") * 4
    expression_matrix <- signal + matrix(stats::rnorm(n_genes * n_cells, sd = 0.35), nrow = n_genes)
    rownames(expression_matrix) <- paste0("gene_", seq_len(n_genes))
    colnames(expression_matrix) <- paste0("cell_", seq_len(n_cells))
    expression_matrix
}

run_case <- function(label, expr) {
    gc()
    profile_file <- NA_character_
    profile_total_file <- NA_character_
    profile_self_file <- NA_character_
    profile_memory_file <- NA_character_
    allocation_count <- NA_integer_
    allocated_bytes <- NA_real_

    if (profile_enabled) {
        dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
        profile_file <- file.path(profile_dir, paste0(safe_filename(label), ".Rprof"))
        profile_memory_file <- file.path(profile_dir, paste0(safe_filename(label), ".Rprofmem"))
        utils::Rprof(profile_file, interval = profile_interval, memory.profiling = TRUE)
        utils::Rprofmem(profile_memory_file, threshold = profile_memory_threshold)
        on.exit(utils::Rprofmem(NULL), add = TRUE)
        on.exit(utils::Rprof(NULL), add = TRUE)
    }

    elapsed <- system.time(force(expr), gcFirst = FALSE)[["elapsed"]]

    if (profile_enabled) {
        utils::Rprofmem(NULL)
        utils::Rprof(NULL)
        profile_summary <- tryCatch(
            utils::summaryRprof(profile_file, memory = "both"),
            error = function(e) NULL
        )

        if (!is.null(profile_summary)) {
            if (!is.null(profile_summary$by.total) && nrow(profile_summary$by.total) > 0) {
                profile_total_file <- file.path(profile_dir, paste0(safe_filename(label), "_by_total.csv"))
                utils::write.csv(profile_summary$by.total, profile_total_file)
            }
            if (!is.null(profile_summary$by.self) && nrow(profile_summary$by.self) > 0) {
                profile_self_file <- file.path(profile_dir, paste0(safe_filename(label), "_by_self.csv"))
                utils::write.csv(profile_summary$by.self, profile_self_file)
            }
        }

        memory_summary <- tryCatch(
            summarize_rprofmem(profile_memory_file),
            error = function(e) NULL
        )
        if (!is.null(memory_summary)) {
            allocation_count <- memory_summary$allocation_count
            allocated_bytes <- memory_summary$allocated_bytes
        }
    }

    data.frame(
        case = label,
        elapsed_seconds = elapsed,
        profile_file = profile_file,
        profile_total_file = profile_total_file,
        profile_self_file = profile_self_file,
        profile_memory_file = profile_memory_file,
        allocation_count = allocation_count,
        allocated_bytes = allocated_bytes,
        stringsAsFactors = FALSE
    )
}

set.seed(20260626)
n_cells <- env_int("CA_BENCH_N_CELLS", 200)
n_partitions <- env_int("CA_BENCH_N_PARTITIONS", 40)
n_clusters <- env_int("CA_BENCH_N_CLUSTERS", 8)
include_pipeline <- env_flag("CA_BENCH_INCLUDE_PIPELINE", TRUE)
include_representative_pipeline <- env_flag("CA_BENCH_INCLUDE_REPRESENTATIVE_PIPELINE", FALSE)
app_dir <- Sys.getenv("CA_BENCH_APP_DIR", "")

partitions <- lapply(seq_len(n_partitions), function(i) {
    sample.int(n_clusters, size = n_cells, replace = TRUE)
})

duplicated_partitions <- c(
    partitions,
    lapply(partitions[seq_len(min(10, length(partitions)))], function(x) {
        match(x, unique(x))
    })
)
unique_winner_partitions <- c(
    partitions,
    rep(list(match(partitions[[1]], unique(partitions[[1]]))), 10)
)

results <- list(
    run_case("getNNmatrix_computeSNN", {
        embedding <- matrix(stats::rnorm(n_cells * 10), ncol = 10)
        nn_idx <- RANN::nn2(embedding, k = 20)$nn.idx
        nn <- getNNmatrix(nn_idx, k = 20, start = 0, prune = -1)$nn
        computeSNN(nn, k = 20, prune = 0)
    }),
    run_case("merge_partitions_freq", {
        merge_partitions(duplicated_partitions, order_logic = "freq")
    }),
    run_case("merge_partitions_freq_unique_winner", {
        merge_partitions(
            unique_winner_partitions,
            order_logic = "freq",
            check_ties = TRUE
        )
    }),
    run_case("merge_partitions_avg_agreement", {
        merge_partitions(duplicated_partitions, order_logic = "avg_agreement")
    }),
    run_case("weighted_element_consistency", {
        weighted_element_consistency(partitions)
    }),
    run_case("element_sim_matrix", {
        element_sim_matrix(partitions)
    }),
    run_case("consensus_cluster", {
        consensus_cluster(
            matrix(stats::rnorm(120 * 12), nrow = 120),
            k_min = 2,
            k_max = 6,
            n_reps = 10,
            p_sample = 0.8,
            p_feature = 0.8
        )
    })
)

if (requireNamespace("Matrix", quietly = TRUE) &&
    requireNamespace("Seurat", quietly = TRUE) &&
    requireNamespace("RANN", quietly = TRUE)) {
    adjacency <- matrix(0, nrow = 20, ncol = 20)
    adjacency[1:10, 1:10] <- 1
    adjacency[11:20, 11:20] <- 1
    diag(adjacency) <- 0
    adjacency <- methods::as(Matrix::Matrix(adjacency, sparse = TRUE), "generalMatrix")
    cell_names <- paste0("cell_", seq_len(nrow(adjacency)))
    rownames(adjacency) <- cell_names
    colnames(adjacency) <- cell_names

    results[[length(results) + 1]] <- run_case("assess_clustering_stability_small", {
        assess_clustering_stability(
            graph_adjacency_matrix = adjacency,
            resolution = c(0.5, 1),
            seed_sequence = c(1, 101, 201),
            clustering_algorithm = 1,
            verbose = FALSE
        )
    })

    pca_embedding <- matrix(stats::rnorm(80 * 8), ncol = 8)
    rownames(pca_embedding) <- paste0("cell_", seq_len(nrow(pca_embedding)))
    colnames(pca_embedding) <- paste0("PC_", seq_len(ncol(pca_embedding)))

    results[[length(results) + 1]] <- run_case("assess_nn_stability_pca_small", {
        suppressWarnings(
            assess_nn_stability(
                embedding = pca_embedding,
                n_neigh_sequence = c(5, 10),
                seed_sequence = c(1, 101, 201),
                graph_reduction_type = "PCA",
                graph_type = 1,
                clustering_algorithm = 1
            )
        )
    })

    if (include_pipeline) {
        expression_matrix <- matrix(
            c(
                5, 5, 0, 0, 0, 0, 0, 0,
                4, 4, 0, 0, 0, 0, 0, 0,
                0, 0, 5, 5, 0, 0, 0, 0,
                0, 0, 4, 4, 0, 0, 0, 0,
                0, 0, 0, 0, 5, 5, 0, 0,
                0, 0, 0, 0, 4, 4, 0, 0,
                0, 0, 0, 0, 0, 0, 5, 5,
                0, 0, 0, 0, 0, 0, 4, 4
            ),
            nrow = 8,
            byrow = TRUE
        )
        rownames(expression_matrix) <- paste0("gene_", seq_len(nrow(expression_matrix)))
        colnames(expression_matrix) <- paste0("cell_", seq_len(ncol(expression_matrix)))

        matrix_processing <- make_exact_pca_processor(default_npcs = 2)

        results[[length(results) + 1]] <- run_case("automatic_stability_assessment_tiny", {
            suppressWarnings(
                automatic_stability_assessment(
                    expression_matrix = expression_matrix,
                    n_repetitions = 2,
                    n_neigh_sequence = c(2),
                    resolution_sequence = c(0.5),
                    features_sets = list(all = rownames(expression_matrix)),
                    steps = list(all = c(4, 8)),
                    seed_sequence = c(1, 101),
                    n_top_configs = 1,
                    matrix_processing = matrix_processing,
                    verbose = FALSE,
                    save_temp = FALSE
                )
            )
        })

        results[[length(results) + 1]] <- run_case("automatic_stability_assessment_tiny_no_umap_embeddings", {
            suppressWarnings(
                automatic_stability_assessment(
                    expression_matrix = expression_matrix,
                    n_repetitions = 2,
                    n_neigh_sequence = c(2),
                    resolution_sequence = c(0.5),
                    features_sets = list(all = rownames(expression_matrix)),
                    steps = list(all = c(4, 8)),
                    seed_sequence = c(1, 101),
                    n_top_configs = 1,
                    matrix_processing = matrix_processing,
                    compute_umap_embeddings = FALSE,
                    verbose = FALSE,
                    save_temp = FALSE
                )
            )
        })
    }

    if (include_representative_pipeline) {
        rep_genes <- env_int("CA_BENCH_REP_GENES", 400)
        rep_cells <- env_int("CA_BENCH_REP_CELLS", 200)
        rep_npcs <- env_int("CA_BENCH_REP_NPCS", 20)
        rep_step_small <- min(env_int("CA_BENCH_REP_STEP_SMALL", 100), rep_genes)
        rep_step_large <- min(env_int("CA_BENCH_REP_STEP_LARGE", 300), rep_genes)
        rep_steps <- sort(unique(c(rep_step_small, rep_step_large)))
        rep_expression <- make_low_rank_expression(rep_genes, rep_cells)
        rep_common <- list(
            expression_matrix = rep_expression,
            n_repetitions = env_int("CA_BENCH_REP_REPETITIONS", 3),
            n_neigh_sequence = c(5, 10),
            resolution_sequence = c(0.5, 1),
            features_sets = list(all = rownames(rep_expression)),
            steps = list(all = rep_steps),
            seed_sequence = seq(from = 1, by = 100, length.out = env_int("CA_BENCH_REP_REPETITIONS", 3)),
            n_top_configs = 1,
            compute_umap_embeddings = FALSE,
            return_nn_idx = FALSE,
            verbose = FALSE,
            save_temp = FALSE
        )

        results[[length(results) + 1]] <- run_case("automatic_stability_assessment_representative_exact_no_umap", {
            suppressWarnings(do.call(
                automatic_stability_assessment,
                c(rep_common, list(matrix_processing = make_exact_pca_processor(rep_npcs)))
            ))
        })

        if (requireNamespace("irlba", quietly = TRUE)) {
            results[[length(results) + 1]] <- run_case("automatic_stability_assessment_representative_irlba_no_umap", {
                suppressWarnings(do.call(
                    automatic_stability_assessment,
                    c(rep_common, list(matrix_processing = make_irlba_pca_processor(rep_npcs)))
                ))
            })
        }
    }
}

if (nzchar(app_dir) &&
    requireNamespace("rhdf5", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE) &&
    requireNamespace("RANN", quietly = TRUE)) {
    expression_h5 <- file.path(app_dir, "expression.h5")
    stability_h5 <- file.path(app_dir, "stability.h5")

    if (file.exists(expression_h5)) {
        app_expression <- rhdf5::h5read(expression_h5, "expression_matrix")
        max_cells <- min(ncol(app_expression), as.integer(Sys.getenv("CA_BENCH_APP_CELLS", "500")))
        max_genes <- min(nrow(app_expression), as.integer(Sys.getenv("CA_BENCH_APP_GENES", "1000")))
        app_expression <- app_expression[seq_len(max_genes), seq_len(max_cells), drop = FALSE]

        results[[length(results) + 1]] <- run_case("app_fixture_pca_nn_snn", {
            app_pca_rank <- min(20, ncol(app_expression), ncol(t(app_expression)) - 1)
            app_k <- min(20, ncol(app_expression) - 1)
            pca <- stats::prcomp(t(app_expression), rank. = app_pca_rank)$x
            nn_idx <- RANN::nn2(pca, k = app_k)$nn.idx
            nn <- getNNmatrix(nn_idx, k = app_k, start = 0, prune = -1)$nn
            computeSNN(nn, k = app_k, prune = 0)
        })

        if (requireNamespace("irlba", quietly = TRUE)) {
            results[[length(results) + 1]] <- run_case("app_fixture_pca_nn_snn_irlba", {
                app_pca_rank <- min(20, ncol(app_expression) - 1, nrow(app_expression) - 1)
                app_k <- min(20, ncol(app_expression) - 1)
                pca <- irlba::prcomp_irlba(t(app_expression), n = app_pca_rank)$x
                nn_idx <- RANN::nn2(pca, k = app_k)$nn.idx
                nn <- getNNmatrix(nn_idx, k = app_k, start = 0, prune = -1)$nn
                computeSNN(nn, k = app_k, prune = 0)
            })
        }
    }

    if (file.exists(stability_h5)) {
        h5_index <- rhdf5::h5ls(stability_h5)
        ecc_paths <- h5_index[h5_index$otype == "H5I_DATASET" & grepl("/ecc$", h5_index$group), ]
        if (nrow(ecc_paths) > 0) {
            selected <- utils::head(seq_len(nrow(ecc_paths)), 12)
            ecc_vectors <- lapply(selected, function(i) {
                rhdf5::h5read(stability_h5, file.path(ecc_paths$group[i], ecc_paths$name[i]))
            })

            results[[length(results) + 1]] <- run_case("app_fixture_stability_ecc_summary", {
                vapply(ecc_vectors, stats::median, numeric(1))
            })
        }
    }
}

results <- do.call(rbind, results)
print(results, row.names = FALSE)

output <- Sys.getenv("CA_BENCH_OUTPUT", "")
if (nzchar(output)) {
    utils::write.csv(results, output, row.names = FALSE)
}
