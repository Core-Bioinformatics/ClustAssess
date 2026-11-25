test_that("native weighted_element_consistency and sim_matrix match pairwise R bit-for-bit", {
    ref_weighted_ecc <- function(clustering_list, weights = NULL, calculate_sim_matrix = FALSE) {
        n_clusterings <- length(clustering_list)
        if (calculate_sim_matrix) {
            ecs_sim_matrix <- matrix(0, nrow = n_clusterings, ncol = n_clusterings)
        }
        if (is.null(weights)) {
            weights <- rep(1, n_clusterings)
        }
        consistency <- 0
        total_weights <- sum(weights)
        norm_factor <- total_weights * (total_weights - 1) / 2
        for (i in seq_len(n_clusterings - 1)) {
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
        list(ecc = consistency, ecs_matrix = ecs_sim_matrix)
    }

    set.seed(2026)
    cases <- list(
        integers = lapply(1:6, function(i) sample.int(7L, 80L, replace = TRUE)),
        negative_large_int = lapply(1:5, function(i) sample(c(-100000L, 0L, 1000000L, 42L), 60L, replace = TRUE)),
        fractional_special = lapply(1:4, function(i) sample(c(-0.0, 0.0, 1.0, 1.0 + 1e-15, 2.5, Inf, -Inf), 50L, replace = TRUE)),
        factors_unused = lapply(1:5, function(i) factor(sample(letters[1:4], 60L, replace = TRUE), levels = letters[6:1])),
        mixed_num_fac = list(
            sample.int(4L, 50L, replace = TRUE),
            as.double(sample.int(5L, 50L, replace = TRUE)),
            factor(sample.int(3L, 50L, replace = TRUE), levels = 5:1),
            sample.int(6L, 50L, replace = TRUE)
        ),
        characters = lapply(1:5, function(i) sample(c("", "a", "b", "c", "alpha"), 70L, replace = TRUE)),
        named_ints = lapply(1:5, function(i) stats::setNames(sample.int(5L, 40L, replace = TRUE), paste0("c", 1:40))),
        singletons_sparse = lapply(1:4, function(i) sample.int(300L, 300L, replace = FALSE))
    )

    for (nm in names(cases)) {
        parts <- cases[[nm]]
        w_list <- list(NULL, seq_along(parts), stats::runif(length(parts), 0.5, 4.0))
        for (w in w_list) {
            for (calc_mat in c(FALSE, TRUE)) {
                expect_identical(
                    weighted_element_consistency(parts, weights = w, calculate_sim_matrix = calc_mat),
                    ref_weighted_ecc(parts, weights = w, calculate_sim_matrix = calc_mat),
                    info = paste(nm, "mat=", calc_mat)
                )
            }
        }
        # Serial element_sim_matrix parity
        expected_mat <- ref_weighted_ecc(parts, calculate_sim_matrix = TRUE)$ecs_matrix
        expect_identical(
            element_sim_matrix(parts, output_type = "matrix"),
            expected_mat,
            info = paste(nm, "sim_matrix")
        )
    }

    # Named weights with length-1 partitions preserve R's exact name semantics
    len1 <- list(1L, 1L)
    named_w <- c(a = 2, b = 3)
    expect_identical(
        weighted_element_consistency(len1, weights = named_w, calculate_sim_matrix = FALSE),
        ref_weighted_ecc(len1, weights = named_w, calculate_sim_matrix = FALSE)
    )
    expect_identical(
        weighted_element_consistency(len1, weights = named_w, calculate_sim_matrix = TRUE),
        ref_weighted_ecc(len1, weights = named_w, calculate_sim_matrix = TRUE)
    )
})

test_that("character weighted_element_consistency and element_sim_matrix leave absent .Random.seed untouched", {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        rm(".Random.seed", envir = .GlobalEnv)
    }
    on.exit({
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
        }
    }, add = TRUE)

    parts <- list(c("", "a", "a", "b"), c("a", "b", "b", ""), c("x", "x", "y", "y"))
    res_wec <- weighted_element_consistency(parts, weights = c(2, 1, 3), calculate_sim_matrix = TRUE)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    res_mat <- element_sim_matrix(parts)
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    expect_identical(res_wec$ecs_matrix, res_mat)
})

test_that("native get_highest_prune_param and get_highest_prune_param_embedding match igraph reference", {
    ref_hpp <- function(nn_matrix, n_neigh) {
        snn <- ClustAssess:::computeSNN(nn_matrix, n_neigh, 0)
        g <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected", weighted = TRUE)
        target <- igraph::components(g)$no
        vals <- sapply(0:n_neigh, function(i) i / (2 * n_neigh - i))
        start_n <- 1
        stop_n <- length(vals)
        prev_g <- g
        while (start_n <= stop_n) {
            middle <- as.integer((stop_n + start_n) / 2)
            cur_g <- igraph::delete_edges(prev_g, which(igraph::E(prev_g)$weight <= vals[middle]))
            if (igraph::components(cur_g)$no > target) {
                stop_n <- middle
                if (start_n == stop_n) start_n <- middle - 1
            } else {
                start_n <- middle + 1
                prev_g <- cur_g
                if (start_n == stop_n) break
            }
        }
        list(prune_value = vals[middle], adj_matrix = ClustAssess:::pruneSNN(snn, vals[middle]))
    }

    set.seed(101)
    # Disconnected 2-blob graph, connected graph, small k (best_c == 0 possible), and weighted fallback
    emb_blobs <- rbind(
        matrix(stats::rnorm(40 * 4, mean = 0), nrow = 40),
        matrix(stats::rnorm(40 * 4, mean = 30), nrow = 40)
    )
    emb_single <- matrix(stats::rnorm(60 * 5), nrow = 60)

    for (emb in list(emb_blobs, emb_single)) {
        for (k in c(2L, 5L, 12L)) {
            nn_idx <- RANN::nn2(emb, k = k)$nn.idx
            nn_mat <- getNNmatrix(nn_idx, k, 0, -1)$nn
            got <- get_highest_prune_param(nn_mat, k)
            exp <- ref_hpp(nn_mat, k)
            expect_identical(got$prune_value, exp$prune_value)
            expect_identical(got$adj_matrix, exp$adj_matrix)
            expect_identical(get_highest_prune_param_embedding(emb, k), exp$prune_value)
        }
    }

    # Non-binary weighted fallback path
    nn_weighted <- getNNmatrix(RANN::nn2(emb_single, k = 6)$nn.idx, 6, 0, -1)$nn
    nn_weighted@x <- stats::runif(length(nn_weighted@x), 0.5, 1.5)
    got_w <- get_highest_prune_param(nn_weighted, 6)
    exp_w <- ref_hpp(nn_weighted, 6)
    expect_identical(got_w$prune_value, exp_w$prune_value)
    expect_identical(got_w$adj_matrix, exp_w$adj_matrix)
})

test_that("native nn_conn_comps matches igraph connected components for PCA and UMAP", {
    foreach::registerDoSEQ()
    set.seed(77)
    emb <- rbind(
        matrix(stats::rnorm(25 * 4, mean = 0), nrow = 25),
        matrix(stats::rnorm(25 * 4, mean = 25), nrow = 25)
    )
    rownames(emb) <- paste0("cell_", seq_len(nrow(emb)))
    ks <- c(2L, 5L, 10L, 15L)

    nn_idx <- RANN::nn2(emb, k = max(ks))$nn.idx
    expected_pca <- stats::setNames(lapply(ks, function(k) {
        g <- igraph::graph_from_adjacency_matrix(getNNmatrix(nn_idx, k, 0, -1)$nn)
        igraph::components(g)$no
    }), as.character(ks))

    got <- get_nn_conn_comps(emb, n_neigh_sequence = ks, include_umap = FALSE)
    expect_identical(got$PCA, expected_pca)
})

test_that("with_nn2_idx_cache reuses exact PCA kNN index once and cleans up on error", {
    foreach::registerDoSEQ()
    set.seed(9)
    emb <- matrix(stats::rnorm(30 * 4), nrow = 30, dimnames = list(paste0("c", 1:30), paste0("PC_", 1:4)))
    calls <- 0L
    real_nn2 <- RANN::nn2
    testthat::local_mocked_bindings(
        nn2 = function(...) {
            calls <<- calls + 1L
            real_nn2(...)
        },
        .package = "RANN"
    )

    ClustAssess:::with_nn2_idx_cache({
        first <- ClustAssess:::parallel_nn2_idx(emb, 5)
        second <- ClustAssess:::parallel_nn2_idx(emb, 5)
        third <- ClustAssess:::parallel_nn2_idx(emb, 5)
    })
    fourth <- ClustAssess:::parallel_nn2_idx(emb, 5)
    expect_identical(first, second)
    expect_identical(first, third)
    expect_identical(first, fourth)
    # Populated once inside scope, then recomputed after scope exits
    expect_identical(calls, 2L)
    expect_false(isTRUE(ClustAssess:::.nn2_idx_cache$active))
    expect_length(ClustAssess:::.nn2_idx_cache$entries, 0L)

    expect_error(
        ClustAssess:::with_nn2_idx_cache({
            ClustAssess:::parallel_nn2_idx(emb, 5)
            stop("boom")
        }),
        "boom"
    )
    expect_false(isTRUE(ClustAssess:::.nn2_idx_cache$active))
    expect_length(ClustAssess:::.nn2_idx_cache$entries, 0L)
})

test_that("automatic_stability_assessment honors save_temp = FALSE when temp_file is supplied", {
    foreach::registerDoSEQ()
    set.seed(11)
    expr <- matrix(stats::rpois(30 * 24, 5) + 0, nrow = 30, ncol = 24,
                   dimnames = list(paste0("g", 1:30), paste0("c", 1:24)))
    tmp_file <- tempfile(fileext = ".rds")
    on.exit(unlink(tmp_file), add = TRUE)

    suppressWarnings(automatic_stability_assessment(
        expression_matrix = expr,
        n_repetitions = 2,
        n_neigh_sequence = c(3, 5),
        resolution_sequence = c(0.6, 0.8),
        features_sets = list(all = rownames(expr)),
        steps = list(all = c(10, 20)),
        n_top_configs = 1,
        umap_arguments = list(n_neighbors = 5, n_epochs = 5, init = "random"),
        algorithms_clustering_assessment = 1,
        verbose = FALSE,
        temp_file = tmp_file,
        save_temp = FALSE,
        compute_umap_embeddings = FALSE
    ))
    expect_false(file.exists(tmp_file))
})

test_that("write_objects handles chunk_size = 1 and nrow < chunk_size without errors or chunk warnings", {
    testthat::skip_if_not_installed("rhdf5")
    tmp <- tempfile("wo_chunks_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

    set.seed(5)
    for (case in list(list(ng = 5L, cs = 1L), list(ng = 12L, cs = 100L), list(ng = 101L, cs = 100L))) {
        expr <- matrix(stats::rpois(case$ng * 8L, 3) + 0, nrow = case$ng, ncol = 8L,
                       dimnames = list(paste0("g", seq_len(case$ng)), paste0("c", 1:8)))
        meta <- data.frame(UMAP_1 = stats::rnorm(8), UMAP_2 = stats::rnorm(8), row.names = colnames(expr))
        subd <- file.path(tmp, paste0("g", case$ng, "_cs", case$cs))
        dir.create(subd)
        expect_no_warning(
            suppressMessages(utils::capture.output(
                write_objects(NULL, expr, meta, subd, chunk_size = case$cs, gene_variance_threshold = -1)
            ))
        )
        expr_h5 <- file.path(subd, "expression.h5")
        mat_h5 <- rhdf5::h5read(expr_h5, "expression_matrix")
        rank_h5 <- rhdf5::h5read(expr_h5, "rank_matrix")
        expected_ranks <- t(apply(mat_h5, 1, function(x) as.integer(rank(x, ties.method = "min"))))
        expect_identical(rank_h5, expected_ranks)
    }
})

test_that("wilcox_test produces valid two-sided p-values matching rank-sum normal approximation", {
    set.seed(19)
    m <- matrix(sample(1:5, 10 * 20, replace = TRUE), nrow = 10, ncol = 20)
    pvals <- ClustAssess:::wilcox_test(m, 10L, max(m))
    expect_length(pvals, 10L)
    expect_true(all(is.finite(pvals) & pvals >= 0 & pvals <= 1))
    # Compare against stats::wilcox.test normal approximation (with continuity correction)
    for (i in seq_len(nrow(m))) {
        ref_p <- suppressWarnings(stats::wilcox.test(m[i, 1:10], m[i, 11:20], exact = FALSE, correct = TRUE)$p.value)
        expect_equal(pvals[i], ref_p, tolerance = 1e-12)
    }
})

test_that("is_package_installed matches installed.packages without scanning every library entry", {
    expect_true(ClustAssess:::is_package_installed("stats"))
    expect_true(ClustAssess:::is_package_installed("Matrix"))
    expect_false(ClustAssess:::is_package_installed("definitely_not_an_installed_package_xyz"))
})

test_that("chunked parallel clustering and NN stability match sequential execution with PSOCK workers", {
    testthat::skip_if_not_installed("doParallel")
    if (!nzchar(system.file(package = "ClustAssess")) ||
        !file.exists(file.path(system.file(package = "ClustAssess"), "Meta", "package.rds"))) {
        testthat::skip("ClustAssess is loaded from source, not an installed package")
    }

    set.seed(42)
    emb <- rbind(
        matrix(stats::rnorm(20 * 4, mean = 0), nrow = 20),
        matrix(stats::rnorm(20 * 4, mean = 8), nrow = 20)
    )
    rownames(emb) <- paste0("cell_", seq_len(nrow(emb)))
    colnames(emb) <- paste0("PC_", 1:4)
    adj <- get_highest_prune_param(
        getNNmatrix(RANN::nn2(emb, k = 8)$nn.idx, 8, 0, -1)$nn,
        8
    )$adj_matrix
    rownames(adj) <- rownames(emb)
    colnames(adj) <- rownames(emb)

    seeds <- c(1L, 101L, 201L, 301L, 401L, 501L)
    res_seq <- c(0.5, 1.0, 1.5)

    foreach::registerDoSEQ()
    seq_nn <- assess_nn_stability(emb, c(5, 8), seed_sequence = seeds, graph_reduction_type = "PCA")
    seq_cl <- assess_clustering_stability(adj, resolution = res_seq, seed_sequence = seeds,
                                          clustering_algorithm = c(1, 2, 3, 6), verbose = FALSE)

    cl <- parallel::makeCluster(2)
    on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
    }, add = TRUE)
    doParallel::registerDoParallel(cl)

    par_nn <- assess_nn_stability(emb, c(5, 8), seed_sequence = seeds, graph_reduction_type = "PCA")
    par_cl <- assess_clustering_stability(adj, resolution = res_seq, seed_sequence = seeds,
                                          clustering_algorithm = c(1, 2, 3, 6), verbose = FALSE)

    expect_identical(par_nn, seq_nn)
    expect_identical(par_cl, seq_cl)
})


test_that("weighted_element_consistency matches 1.2.0 on non-logical flags and bytes-encoded labels", {
    parts <- list(c(1L, 1L, 2L, 2L, 3L, 3L), c(1L, 2L, 2L, 3L, 3L, 1L), c(2L, 2L, 1L, 1L, 3L, 3L))
    expected_mat <- weighted_element_consistency(parts, calculate_sim_matrix = TRUE)
    # Numeric truthy flag (1L) must behave identically to 1.2.0 if (calculate_sim_matrix)
    expect_identical(weighted_element_consistency(parts, calculate_sim_matrix = 1L), expected_mat)
    expect_error(weighted_element_consistency(parts, calculate_sim_matrix = NA))
    expect_error(weighted_element_consistency(parts, calculate_sim_matrix = c(TRUE, FALSE)))
    expect_error(weighted_element_consistency(parts, calculate_sim_matrix = logical(0)))
    expect_error(weighted_element_consistency(parts, calculate_sim_matrix = "TRUE"))

    # Bytes-encoded strings must fall back cleanly instead of failing translation
    b1 <- c("a", "\xff", "a", "\xff")
    Encoding(b1) <- "bytes"
    b2 <- c("x", "x", "y", "y")
    res_b <- weighted_element_consistency(list(b1, b2), calculate_sim_matrix = TRUE)
    expect_true(is.list(res_b))
    expect_identical(dim(res_b$ecs_matrix), c(2L, 2L))
    expect_identical(element_sim_matrix(list(b1, b2)), res_b$ecs_matrix)
})

test_that("get_highest_prune_param terminates on negative SNN weights and preserves fractional/zero n_neigh", {
    # Negative SNN weights from weighted matrix(2, 2, 2) must terminate with a clear error (never hang)
    x <- methods::as(methods::as(matrix(2, 2, 2), "generalMatrix"), "CsparseMatrix")
    expect_error(get_highest_prune_param(x, 2L), "No valid pruning threshold")

    # Fractional n_neigh = 2.5 and n_neigh = 0 use exact 1.2.0 R fallback
    a <- methods::as(Matrix::sparseMatrix(i = c(1:3, 1:2), j = c(1:3, 2:3), x = 1, dims = c(3, 3)), "dgCMatrix")
    res_frac <- get_highest_prune_param(a, 2.5)
    expect_equal(res_frac$prune_value, 1 / (2 * 2.5 - 1))
    res_zero <- get_highest_prune_param(a, 0L)
    expect_true(is.nan(res_zero$prune_value))
})

test_that("nn_conn_comps_from_idx_cpp safely handles NA and NaN entries without overflow", {
    expect_identical(
        ClustAssess:::nn_conn_comps_from_idx_cpp(matrix(c(NA_integer_, 2L), ncol = 1L), 1L),
        2
    )
    expect_identical(
        ClustAssess:::nn_conn_comps_from_idx_cpp(matrix(c(NaN, 2), ncol = 1L), 1L),
        2
    )
})

test_that("wilcox_test does not overflow signed 32-bit int for n1 >= 46341 and validates bounds", {
    # Equal ranks across n1 = 46341 and n2 = 46341 should give p = 1 (not 0 from int overflow)
    n1 <- 46341L
    ranks <- matrix(rep(1L, 2L * n1), nrow = 1L)
    expect_equal(ClustAssess:::wilcox_test(ranks, n1, 1L), 1)

    expect_error(ClustAssess:::wilcox_test(matrix(c(1L, 5L), nrow = 1L), 1L, 2L), "rank_values")
    expect_error(ClustAssess:::wilcox_test(matrix(c(1L, 2L), nrow = 1L), -1L, 2L), "n1")
})

test_that("is_package_installed matches installed.packages on empty, NA, NULL, and vector inputs", {
    expect_identical(ClustAssess:::is_package_installed(""), FALSE)
    expect_identical(ClustAssess:::is_package_installed(NA_character_), FALSE)
    expect_identical(ClustAssess:::is_package_installed(character(0)), logical(0))
    expect_identical(ClustAssess:::is_package_installed(NULL), logical(0))
    expect_identical(
        ClustAssess:::is_package_installed(c("stats", "definitely_not_installed_xyz", "")),
        c(TRUE, FALSE, FALSE)
    )
})

test_that("parallel_nn2_idx cache distinguishes fractional k from integer k", {
    foreach::registerDoSEQ()
    emb <- matrix(c(0, 1, 2, 0, 1, 2), nrow = 3L, ncol = 2L)
    ClustAssess:::with_nn2_idx_cache({
        ok <- ClustAssess:::parallel_nn2_idx(emb, 3)
        expect_identical(dim(ok), c(3L, 3L))
        # k = 3.5 must not falsely hit the k = 3 cache entry (RANN::nn2 errors when k > ND)
        expect_error(ClustAssess:::parallel_nn2_idx(emb, 3.5))
    })
})


test_that("ECS rejects NA element names consistently and matches fallback on classed weights", {
    a <- stats::setNames(c(1L, 1L, 2L, 2L), c("c1", NA, "c3", "c4"))
    b <- stats::setNames(c(1L, 2L, 2L, 1L), c("c1", NA, "c3", "c4"))
    expect_error(element_sim_elscore(a, b))
    expect_error(weighted_element_consistency(list(a, b), calculate_sim_matrix = TRUE))
    expect_error(element_sim_matrix(list(a, b)))

    # Classed numeric weights must behave consistently regardless of calculate_sim_matrix type
    p1 <- c(1L, 1L, 2L, 2L, 3L, 3L)
    p2 <- c(1L, 2L, 2L, 3L, 3L, 1L)
    w_classed <- structure(c(1, 3), class = "review_num")
    expect_identical(
        weighted_element_consistency(list(p1, p2), weights = w_classed, calculate_sim_matrix = TRUE),
        weighted_element_consistency(list(p1, p2), weights = w_classed, calculate_sim_matrix = 1L)
    )
})

test_that("direct pairwise disjointECS and disjointECSaverage sparse path matches dense bit-for-bit", {
    set.seed(321)
    # 400 singletons -> 400 * 400 = 160,000 > max(4 * 400, 65536), exercising sparse path!
    n <- 400L
    p1 <- sample.int(n, n, replace = FALSE)
    p2 <- sample.int(n, n, replace = FALSE)
    # Also test partial ties with high K (300 clusters across 400 elements)
    q1 <- c(sample.int(260L, 260L, replace = FALSE), sample.int(260L, 140L, replace = TRUE))
    q2 <- c(sample.int(260L, 260L, replace = FALSE), sample.int(260L, 140L, replace = TRUE))

    for (pair in list(list(p1, p2), list(q1, q2))) {
        a <- pair[[1]]; b <- pair[[2]]
        # Reference using myContTable (dense)
        ct <- ClustAssess:::myContTable(a, b, min(a), min(b))
        rs <- rowSums(ct); cs <- colSums(ct)
        uniq <- ct / outer(rs, cs, pmax)
        ref_el <- uniq[cbind(a - min(a) + 1L, b - min(b) + 1L)]
        ref_avg <- 0
        for (i in seq_len(nrow(ct))) {
            for (j in rev(seq_len(ncol(ct)))) {
                m <- max(rs[i], cs[j])
                if (m > 0) ref_avg <- ref_avg + ct[i, j] * ct[i, j] / (as.double(n) * m)
            }
        }
        expect_identical(ClustAssess:::disjointECS(a, b), ref_el)
        expect_identical(ClustAssess:::disjointECSaverage(a, b), ref_avg)
    }
})

test_that("get_highest_prune_param rejects non-finite weights and embedding/NN builders validate bounds", {
    m_na <- methods::as(Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(1, NA_real_, 1), dims = c(3, 3)), "dgCMatrix")
    expect_error(get_highest_prune_param(m_na, 2L), "missing or non-finite")
    m_inf <- methods::as(Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(1, Inf, 1), dims = c(3, 3)), "dgCMatrix")
    expect_error(get_highest_prune_param(m_inf, 2L), "missing or non-finite")

    emb <- matrix(stats::rnorm(12), nrow = 4, ncol = 3)
    expect_error(get_highest_prune_param_embedding(emb, 2.5), "whole number")
    expect_error(get_highest_prune_param_embedding(emb, 0), "whole number")
    expect_error(get_highest_prune_param_embedding(emb, 5), "whole number")
    expect_error(getNNmatrix(matrix(c(1, 99, 1, 2), nrow = 2), 2, 0, -1), "1-based indices")
    expect_error(getNNmatrix(matrix(c(1, 2, 1, 2), nrow = 2), 2, -1, -1), "start")
})

test_that("parallel_nn2_idx handles ncores > nrow(embedding) and bypasses cache for SharedObject", {
    emb_small <- matrix(c(0, 1, 2, 0, 1, 2), nrow = 3L, ncol = 2L)
    # Mock getDoParWorkers() = 4 > nrow(emb_small) = 3 while executing sequentially
    foreach::registerDoSEQ()
    testthat::local_mocked_bindings(
        getDoParWorkers = function() 4L,
        .package = "foreach"
    )
    res_small <- ClustAssess:::parallel_nn2_idx(emb_small, 2L)
    expect_identical(res_small, RANN::nn2(emb_small, k = 2L)$nn.idx)

    if (requireNamespace("SharedObject", quietly = TRUE)) {
        sh_emb <- SharedObject::share(matrix(as.numeric(1:12), nrow = 6L, ncol = 2L))
        on.exit(sh_emb <- SharedObject::unshare(sh_emb), add = TRUE)
        SharedObject::setCopyOnWrite(sh_emb, FALSE)
        ClustAssess:::with_nn2_idx_cache({
            r1 <- ClustAssess:::parallel_nn2_idx(sh_emb, 3L)
            sh_emb[1, 1] <- 1000
            r2 <- ClustAssess:::parallel_nn2_idx(sh_emb, 3L)
            expect_identical(r2, RANN::nn2(as.matrix(sh_emb), k = 3L)$nn.idx)
            expect_false(identical(r1, r2))
        })
    }
})

test_that("automatic_stability_assessment keeps per-set n_top_configs and handles single-step sets", {
    foreach::registerDoSEQ()
    set.seed(23)
    expr <- matrix(stats::rpois(40 * 24, 5) + 0, nrow = 40, ncol = 24,
                   dimnames = list(paste0("g", 1:40), paste0("c", 1:24)))
    genes <- rownames(expr)
    # Single-step feature set first (1 eligible step) followed by multi-step set (4 steps) with n_top_configs = 2
    res <- suppressWarnings(automatic_stability_assessment(
        expression_matrix = expr,
        n_repetitions = 2,
        n_neigh_sequence = c(3, 5),
        resolution_sequence = c(0.6, 0.8),
        features_sets = list(small = genes[1:10], big = genes),
        steps = list(small = 10L, big = c(10L, 20L, 30L, 40L)),
        n_top_configs = 2L,
        umap_arguments = list(n_neighbors = 5, n_epochs = 5, init = "random"),
        algorithms_clustering_assessment = 1,
        verbose = FALSE,
        save_temp = FALSE,
        compute_umap_embeddings = FALSE
    ))
    expect_identical(setdiff(names(res$small), "feature_list"), "10")
    expect_gte(length(setdiff(names(res$big), "feature_list")), 2L)
})
