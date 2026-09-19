#!/usr/bin/env Rscript

fixture_h5 <- Sys.getenv("CA_COMPARE_FIXTURE", unset = "")
if (!nzchar(fixture_h5)) {
    stop("Set CA_COMPARE_FIXTURE to an expression.h5 fixture file.")
}
if (!file.exists(fixture_h5)) {
    stop("Missing fixture: ", fixture_h5)
}

env_int <- function(name, default) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (is.na(value) || !nzchar(value)) default else as.integer(value)
}

env_num_seq <- function(name, default) {
    value <- Sys.getenv(name, unset = default)
    as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

env_flag <- function(name, default = FALSE) {
    value <- Sys.getenv(name, unset = if (default) "true" else "false")
    tolower(value) %in% c("1", "true", "yes", "y")
}

safe_filename <- function(x) {
    gsub("[^A-Za-z0-9_.-]+", "_", x)
}

summarize_rprofmem <- function(memory_file) {
    lines <- readLines(memory_file, warn = FALSE)
    sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    sizes <- sizes[!is.na(sizes)]
    list(
        allocation_count = length(sizes),
        allocated_bytes = sum(sizes),
        largest_allocation_bytes = if (length(sizes) == 0) 0 else max(sizes)
    )
}

lib_path <- Sys.getenv("CA_COMPARE_LIB", unset = "")
if (nzchar(lib_path)) {
    .libPaths(c(lib_path, .libPaths()))
}

suppressPackageStartupMessages({
    library(ClustAssess)
    library(foreach)
})

if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package `rhdf5` is required to run this benchmark.")
}

read_fixture <- function(ngenes, ncells) {
    genes <- rhdf5::h5read(fixture_h5, "genes")
    cells <- rhdf5::h5read(fixture_h5, "cells")
    avg <- rhdf5::h5read(fixture_h5, "average_expression")
    idx <- order(avg, decreasing = TRUE)[seq_len(min(ngenes, length(avg)))]
    cidx <- seq_len(min(ncells, length(cells)))
    expr <- rhdf5::h5read(fixture_h5, "expression_matrix", index = list(idx, cidx))
    rownames(expr) <- genes[idx]
    colnames(expr) <- cells[cidx]
    expr
}

make_inputs <- function(expr) {
    vars <- matrixStats::rowVars(expr)
    list(
        features_sets = list(
            Most_Abundant = rownames(expr),
            Highly_Variable = rownames(expr)[order(vars, decreasing = TRUE)]
        ),
        steps = list(
            Most_Abundant = c(500L, 1000L, min(1500L, nrow(expr))),
            Highly_Variable = c(500L, 1000L, min(1500L, nrow(expr)))
        )
    )
}

register_workers <- function(workers) {
    if (workers <= 1L) {
        foreach::registerDoSEQ()
        return(NULL)
    }
    if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("Package `doParallel` is required for workers > 1.")
    }
    cl <- parallel::makeCluster(workers)
    doParallel::registerDoParallel(cl)
    cl
}

setup_tracing <- function() {
    if (!env_flag("CA_COMPARE_TRACE", FALSE)) {
        return(NULL)
    }

    `%||%` <- function(x, y) if (is.null(x)) y else x
    profile <- new.env(parent = emptyenv())
    profile$totals <- list()
    profile$counts <- list()
    profile$weighted_ecc_shape_totals <- list()
    profile$weighted_ecc_shape_counts <- list()
    profile$add <- function(name, elapsed) {
        profile$totals[[name]] <- (profile$totals[[name]] %||% 0) + elapsed
        profile$counts[[name]] <- (profile$counts[[name]] %||% 0L) + 1L
    }
    profile$add_weighted_ecc_shape <- function(shape, elapsed) {
        profile$weighted_ecc_shape_totals[[shape]] <-
            (profile$weighted_ecc_shape_totals[[shape]] %||% 0) + elapsed
        profile$weighted_ecc_shape_counts[[shape]] <-
            (profile$weighted_ecc_shape_counts[[shape]] %||% 0L) + 1L
    }

    trace_one <- function(fun) {
        if (!exists(fun, envir = asNamespace("ClustAssess"), inherits = FALSE)) {
            return(invisible(FALSE))
        }
        if (fun == "weighted_element_consistency") {
            trace(
                fun,
                tracer = substitute({
                    .ca_t0 <- proc.time()[["elapsed"]]
                    .ca_shape <- paste(
                        paste0("elements=", length(clustering_list[[1]])),
                        paste0("clusterings=", length(clustering_list)),
                        paste0("matrix=", isTRUE(calculate_sim_matrix)),
                        sep = " "
                    )
                }),
                exit = substitute({
                    .ca_elapsed <- proc.time()[["elapsed"]] - .ca_t0
                    PROFILE$add(FUN, .ca_elapsed)
                    PROFILE$add_weighted_ecc_shape(.ca_shape, .ca_elapsed)
                }, list(PROFILE = profile, FUN = fun)),
                where = asNamespace("ClustAssess"),
                print = FALSE
            )
        } else {
            trace(
                fun,
                tracer = substitute(.ca_t0 <- proc.time()[["elapsed"]]),
                exit = substitute(
                    PROFILE$add(FUN, proc.time()[["elapsed"]] - .ca_t0),
                    list(PROFILE = profile, FUN = fun)
                ),
                where = asNamespace("ClustAssess"),
                print = FALSE
            )
        }
        TRUE
    }

    for (fun in strsplit(
        Sys.getenv(
            "CA_COMPARE_TRACE_FUNCTIONS",
            "merge_resolutions,summarize_clustering_memberships,merge_partitions,weighted_element_consistency,merge_identical_partitions"
        ),
        ",",
        fixed = TRUE
    )[[1]]) {
        trace_one(trimws(fun))
    }

    profile
}

print_profile <- function(profile) {
    if (is.null(profile)) {
        return(invisible(NULL))
    }
    cat("function_profile\n")
    for (name in names(profile$totals)) {
        cat(
            name,
            "calls=", profile$counts[[name]],
            "elapsed=", sprintf("%.3f", profile$totals[[name]]),
            "\n",
            sep = " "
        )
    }
    if (length(profile$weighted_ecc_shape_totals) > 0) {
        cat("weighted_ecc_shape_profile\n")
        ordered_shapes <- names(sort(
            unlist(profile$weighted_ecc_shape_totals),
            decreasing = TRUE
        ))
        for (shape in ordered_shapes) {
            cat(
                shape,
                "calls=", profile$weighted_ecc_shape_counts[[shape]],
                "elapsed=", sprintf(
                    "%.3f",
                    profile$weighted_ecc_shape_totals[[shape]]
                ),
                "\n",
                sep = " "
            )
        }
    }
}

setup_rprof <- function(label) {
    if (!env_flag("CA_COMPARE_RPROF", FALSE)) {
        return(NULL)
    }

    profile_dir <- Sys.getenv(
        "CA_COMPARE_PROFILE_DIR",
        file.path(tempdir(), "clustassess-branch-profiles")
    )
    dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
    base <- file.path(profile_dir, safe_filename(label))
    profile_file <- paste0(base, ".Rprof")
    memory_file <- paste0(base, ".Rprofmem")

    utils::Rprof(
        profile_file,
        interval = as.numeric(Sys.getenv("CA_COMPARE_RPROF_INTERVAL", "0.01")),
        memory.profiling = TRUE
    )
    utils::Rprofmem(
        memory_file,
        threshold = as.numeric(Sys.getenv("CA_COMPARE_RPROFMEM_THRESHOLD", "0"))
    )

    list(
        profile_file = profile_file,
        memory_file = memory_file,
        total_file = paste0(base, "_by_total.csv"),
        self_file = paste0(base, "_by_self.csv")
    )
}

stop_rprof <- function(rprof) {
    if (is.null(rprof)) {
        return(invisible(NULL))
    }

    utils::Rprofmem(NULL)
    utils::Rprof(NULL)
    summary <- tryCatch(
        utils::summaryRprof(rprof$profile_file, memory = "both"),
        error = function(e) NULL
    )
    if (!is.null(summary)) {
        if (!is.null(summary$by.total) && nrow(summary$by.total) > 0) {
            utils::write.csv(summary$by.total, rprof$total_file)
        }
        if (!is.null(summary$by.self) && nrow(summary$by.self) > 0) {
            utils::write.csv(summary$by.self, rprof$self_file)
        }
    }

    memory <- summarize_rprofmem(rprof$memory_file)
    cat("rprof_files\n")
    cat("profile=", rprof$profile_file, "\n", sep = "")
    cat("by_total=", rprof$total_file, "\n", sep = "")
    cat("by_self=", rprof$self_file, "\n", sep = "")
    cat("memory=", rprof$memory_file, "\n", sep = "")
    cat(
        "rprofmem_summary allocations=", memory$allocation_count,
        " allocated_mb=", sprintf("%.1f", memory$allocated_bytes / 1024^2),
        " largest_mb=", sprintf("%.1f", memory$largest_allocation_bytes / 1024^2),
        "\n",
        sep = ""
    )
}

run_assessment <- function(expr, inputs, workers, mode) {
    cl <- register_workers(workers)
    on.exit({
        if (!is.null(cl)) {
            parallel::stopCluster(cl)
        }
        foreach::registerDoSEQ()
    }, add = TRUE)

    common_args <- list(
        expression_matrix = expr,
        n_repetitions = env_int("CA_COMPARE_REPETITIONS", 20L),
        n_neigh_sequence = env_num_seq("CA_COMPARE_N_NEIGH", "5,15,25,35,45"),
        resolution_sequence = env_num_seq("CA_COMPARE_RESOLUTION", paste(seq(0.5, 2.5, 0.1), collapse = ",")),
        features_sets = inputs$features_sets,
        steps = inputs$steps,
        seed_sequence = seq(from = 1L, by = 100L, length.out = env_int("CA_COMPARE_REPETITIONS", 20L)),
        n_top_configs = env_int("CA_COMPARE_TOP_CONFIGS", 2L),
        umap_arguments = list(
            n_neighbors = 30L,
            n_threads = 1L,
            n_sgd_threads = 1L,
            n_epochs = env_int("CA_COMPARE_UMAP_EPOCHS", 10L)
        ),
        verbose = FALSE,
        save_temp = FALSE
    )

    if (mode == "current_fast") {
        if (!exists("make_pca_matrix_processing", mode = "function")) {
            stop("current_fast mode requires make_pca_matrix_processing().")
        }
        common_args$matrix_processing <- make_pca_matrix_processing(
            actual_npcs = env_int("CA_COMPARE_NPCS", 30L),
            method = "irlba",
            seed = 1L
        )
        common_args$compute_umap_embeddings <- FALSE
        common_args$return_nn_idx <- FALSE
    } else if (mode != "strict") {
        stop("CA_COMPARE_MODE must be `strict` or `current_fast`.")
    }

    assessment_result <- NULL
    elapsed <- system.time(assessment_result <- suppressWarnings(do.call(
        automatic_stability_assessment,
        common_args
    )))

    result_file <- Sys.getenv("CA_COMPARE_SAVE_RESULT", unset = "")
    result_md5 <- NA_character_
    result_size_bytes <- as.numeric(utils::object.size(assessment_result))
    if (nzchar(result_file)) {
        saveRDS(assessment_result, result_file, version = 2)
        result_md5 <- unname(tools::md5sum(result_file))
    }

    data.frame(
        label = Sys.getenv("CA_COMPARE_LABEL", unset = "unnamed"),
        mode = mode,
        package_version = as.character(utils::packageVersion("ClustAssess")),
        workers = workers,
        genes = nrow(expr),
        cells = ncol(expr),
        elapsed = unname(elapsed[["elapsed"]]),
        user = unname(elapsed[["user.self"]]),
        sys = unname(elapsed[["sys.self"]]),
        result_size_bytes = result_size_bytes,
        result_file = if (nzchar(result_file)) result_file else NA_character_,
        result_md5 = result_md5,
        stringsAsFactors = FALSE
    )
}

expr <- read_fixture(
    ngenes = env_int("CA_COMPARE_GENES", 1500L),
    ncells = env_int("CA_COMPARE_CELLS", 600L)
)
inputs <- make_inputs(expr)
workers <- env_int("CA_COMPARE_WORKERS", 8L)
mode <- Sys.getenv("CA_COMPARE_MODE", unset = "strict")
profile <- setup_tracing()
label <- Sys.getenv("CA_COMPARE_LABEL", unset = "unnamed")

cat("fixture=", fixture_h5, "\n", sep = "")
cat("mode=", mode, " genes=", nrow(expr), " cells=", ncol(expr), "\n", sep = "")
cat("settings=repetitions=", env_int("CA_COMPARE_REPETITIONS", 20L),
    " neigh=", paste(env_num_seq("CA_COMPARE_N_NEIGH", "5,15,25,35,45"), collapse = ","),
    " resolutions=", length(env_num_seq("CA_COMPARE_RESOLUTION", paste(seq(0.5, 2.5, 0.1), collapse = ","))),
    " top_configs=", env_int("CA_COMPARE_TOP_CONFIGS", 2L),
    " workers=", workers, "\n", sep = "")

rprof <- setup_rprof(label)
assessment_summary <- tryCatch(
    run_assessment(expr, inputs, workers, mode),
    finally = stop_rprof(rprof)
)
rprof <- NULL
print(assessment_summary, row.names = FALSE)
print_profile(profile)
