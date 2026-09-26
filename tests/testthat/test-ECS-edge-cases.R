test_that("ECS handles exact matches with relabelled and non-contiguous labels", {
    x <- c(10, 10, 20, 20, 30, 30)
    y <- c("b", "b", "a", "a", "c", "c")
    z <- factor(c(-1, -1, 7, 7, 99, 99))

    expect_equal(element_sim_elscore(x, y), rep(1, length(x)))
    expect_equal(element_sim_elscore(x, z), rep(1, length(x)))
    expect_equal(element_sim(x, y), 1)
    expect_equal(element_sim(x, z), 1)
})

test_that("average ECS ignores unused rows and columns in label ranges", {
    x <- c(1L, 3L, 1L, 3L)

    expect_equal(element_sim_elscore(x, x), rep(1, length(x)))
    expect_equal(element_sim(x, x), 1)
    expect_equal(
        element_sim_matrix(list(x, x)),
        matrix(1, nrow = 2, ncol = 2)
    )
})

test_that("ECS handles crossed partitions with known element scores", {
    x <- c(1, 1, 2, 2)
    y <- c(1, 2, 1, 2)

    expect_equal(element_sim_elscore(x, y), rep(0.5, 4))
    expect_equal(element_sim(x, y), 0.5)
})

test_that("ECS and ECC preserve distinct fractional and precision-sensitive labels", {
    y <- c(1L, 2L, 1L, 2L)
    crossed <- c(1L, 1L, 2L, 2L)
    for (labels in list(c(1.2, 1.7), c(1, 1 + 1e-15), c(-1.7, -1.2))) {
        x <- labels[y]
        expect_equal(element_sim(x, y), 1)
        expect_equal(element_sim_elscore(x, y), rep(1, 4))
        expect_equal(element_sim(x, crossed), 0.5)
        expect_equal(element_sim_elscore(x, as.character(y)), rep(1, 4))
        expect_equal(element_consistency(list(x, y, x)), rep(1, 4))

        for (n_partitions in c(2L, 5L)) {
            memberships <- c(list(x), rep(list(y), n_partitions - 1L))
            weighted <- weighted_element_consistency(
                memberships,
                weights = seq_len(n_partitions),
                calculate_sim_matrix = TRUE
            )
            expect_equal(weighted$ecc, rep(1, 4))
            expect_equal(weighted$ecs_matrix, matrix(1, n_partitions, n_partitions))
            expect_equal(element_sim_matrix(memberships), matrix(1, n_partitions, n_partitions))
        }

        merged <- merge_partitions(list(x, y), ecs_thresh = 0.99, order_logic = "none")
        expect_length(merged$partitions, 1)
        expect_equal(merged$partitions[[1]]$freq, 2)
    }
})

test_that("flat scoring uses compact ids for negative and widely spaced labels", {
    y <- c(1L, 2L, 1L, 2L)
    crossed <- c(1L, 1L, 2L, 2L)
    label_sets <- list(
        c(-1000000000L, 1000000000L),
        c(1L, .Machine$integer.max),
        c(-1e100, 1e100)
    )
    for (labels in label_sets) {
        x <- setNames(labels[y], letters[1:4])
        expect_equal(element_sim(x, y), 1)
        expect_equal(element_sim_elscore(x, y), setNames(rep(1, 4), names(x)))
        expect_equal(element_sim(x, crossed), 0.5)

        for (n_partitions in c(2L, 5L)) {
            memberships <- c(list(x), rep(list(setNames(y, names(x))), n_partitions - 1L))
            expect_equal(
                unname(weighted_element_consistency(memberships)),
                rep(1, 4)
            )
            expect_equal(element_sim_matrix(memberships), matrix(1, n_partitions, n_partitions))
        }
    }
})

test_that("flat label encoding preserves safe ids and compacts other labels", {
    for (dense in list(c(2L, 1L, 2L, 1L), c(2, 1, 2, 1))) {
        expect_identical(ClustAssess:::encode_flat_disjoint_membership(dense), as.integer(dense))
    }
    for (x in list(c(2L, 4L, 2L, 4L), c(-2L, -1L, -2L, -1L),
                   c(1.2, 1.7, 1.2, 1.7), factor(c("b", "a", "b", "a")))) {
        expect_identical(ClustAssess:::encode_flat_disjoint_membership(x), c(1L, 2L, 1L, 2L))
    }
    encoded <- ClustAssess:::encode_flat_disjoint_memberships(
        rep(list(c(1.2, 1.7, 1.2, 1.7)), 5)
    )
    expect_identical(encoded, matrix(rep(c(1L, 2L, 1L, 2L), 5), ncol = 5))
})

test_that("flat element agreement preserves precise and empty-string labels", {
    inputs <- list(
        rep(c(1, 1 + 1e-15), 2),
        c(1.2, 1.7, 1.2, 1.7),
        c("", "cluster", "", "cluster"),
        factor(c("", "cluster", "", "cluster"))
    )
    comparisons <- list(c(1L, 2L, 1L, 2L), c("a", "a", "b", "b"))
    for (reference in inputs) {
        names(reference) <- letters[1:4]
        expected <- Reduce(`+`, lapply(comparisons, function(x) {
            unname(element_sim_elscore(reference, x))
        })) / length(comparisons)
        actual <- element_agreement(reference, comparisons)
        expect_equal(actual, expected)
        expect_null(names(actual))
        expect_equal(element_agreement(comparisons[[1]], list(reference)), rep(1, 4))
    }
    expect_error(element_agreement(inputs[[1]], list(c(1, 2))), "same length")
})

test_that("flat membership APIs reject missing labels before early returns", {
    missing_inputs <- list(
        c(1, NA_real_, 2, 2), c(1, NaN, 2, 2), c(1L, NA_integer_, 2L, 2L),
        c("a", NA_character_, "b", "b"), factor(c("a", NA, "b", "b"))
    )
    valid <- c(1L, 1L, 2L, 2L)
    for (bad in missing_inputs) {
        calls <- list(
            function() element_sim(bad, valid),
            function() element_sim(valid, bad),
            function() element_sim_elscore(bad, valid),
            function() element_sim_matrix(list(bad)),
            function() element_sim_matrix(list(valid, bad)),
            function() weighted_element_consistency(list(bad)),
            function() weighted_element_consistency(list(valid, bad)),
            function() weighted_element_consistency(c(rep(list(valid), 4), list(bad))),
            function() element_consistency(list(bad)),
            function() element_consistency(list(bad, bad)),
            function() merge_partitions(list(bad)),
            function() merge_partitions(list(bad, bad)),
            function() merge_partitions(list(list(mb = bad, freq = 2))),
            function() merge_partitions(list(bad), ecs_thresh = 0.99),
            function() element_agreement(bad, list(valid)),
            function() element_agreement(valid, list(bad)),
            function() ClustAssess:::create_clustering(bad)
        )
        for (call in calls) {
            expect_error(call(), "Flat cluster labels must not contain missing values")
        }
    }
})

test_that("empty strings remain valid flat cluster labels", {
    x <- c("", "b", "", "b")
    y <- c(1L, 2L, 1L, 2L)
    expect_equal(element_sim(x, y), 1)
    expect_equal(element_sim_elscore(x, y), rep(1, 4))
    expect_equal(element_sim_matrix(list(x, y)), matrix(1, 2, 2))
    expect_equal(weighted_element_consistency(rep(list(x), 5)), rep(1, 4))
    expect_equal(element_consistency(list(x, y)), rep(1, 4))
    expect_length(merge_partitions(list(x, y))$partitions, 1)
    expect_equal(element_agreement(x, list(y)), rep(1, 4))
})

test_that("flat element agreement works on installed PSOCK workers", {
    package_path <- system.file(package = "ClustAssess")
    skip_if(!nzchar(package_path), "ClustAssess is not installed in a library path")
    skip_if(
        !file.exists(file.path(package_path, "Meta", "package.rds")),
        "ClustAssess is loaded from source, not an installed package"
    )
    cl <- parallel::makePSOCKcluster(2)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    reference <- rep(c(1, 1 + 1e-15), 2)
    comparisons <- list(c("", "b", "", "b"), c(1L, 1L, 2L, 2L))
    results <- parallel::clusterCall(cl, function(reference, comparisons) {
        list(
            agreement = ClustAssess::element_agreement(reference, comparisons),
            missing = tryCatch(
                ClustAssess::element_agreement(c(1, NA, 2, 2), comparisons),
                error = conditionMessage
            )
        )
    }, reference, comparisons)
    for (result in results) {
        expect_equal(result$agreement, rep(0.75, 4))
        expect_match(result$missing, "Flat cluster labels must not contain missing values")
    }
})

test_that("average ECS remains correct when cluster-size products exceed integer range", {
    for (n_cells in c(50000L, 65536L)) {
        x <- rep(1L, n_cells)
        y <- rep(1:2, each = n_cells / 2)

        expect_equal(element_sim(x, x), 1)
        expect_equal(element_sim(x, y), mean(element_sim_elscore(x, y)))
        expect_equal(
            element_sim_matrix(list(x, y)),
            matrix(c(1, 0.5, 0.5, 1), nrow = 2)
        )
    }
})

test_that("threshold merging handles large nearly identical partitions", {
    x <- rep(1L, 50000)
    y <- x
    y[length(y)] <- 2L

    merged <- merge_partitions(list(x, y), ecs_thresh = 0.99, order_logic = "none")

    expect_length(merged$partitions, 1)
    expect_equal(merged$partitions[[1]]$freq, 2)
    expect_equal(merged$ecc, rep(1, length(x)))
})

test_that("weighted ECC matches ECC on expanded repeated partitions", {
    x <- c(1, 1, 2, 2)
    y <- c(1, 2, 1, 2)

    weighted <- weighted_element_consistency(
        clustering_list = list(x, y),
        weights = c(2, 1)
    )
    expanded <- element_consistency(list(x, x, y))

    expect_equal(weighted, expanded)
})

test_that("weighted ECC requires at least two clusterings", {
    x <- c(1, 1, 2, 2)

    expect_error(
        weighted_element_consistency(list(x)),
        "At least two clusterings are required"
    )
})

manual_weighted_element_consistency <- function(clustering_list,
                                                weights,
                                                calculate_sim_matrix = FALSE) {
    n_clusterings <- length(clustering_list)
    norm_factor <- sum(weights) * (sum(weights) - 1) / 2
    ecc <- rep(0, length(clustering_list[[1]]))

    if (calculate_sim_matrix) {
        ecs_matrix <- diag(1, n_clusterings)
    }

    for (i in seq_len(n_clusterings)) {
        ecc <- ecc + weights[i] * (weights[i] - 1) / 2 / norm_factor

        if (i == n_clusterings) {
            next
        }

        for (j in seq(from = i + 1, to = n_clusterings)) {
            ecs <- element_sim_elscore(clustering_list[[i]], clustering_list[[j]])
            ecc <- ecc + ecs * weights[i] * weights[j] / norm_factor

            if (calculate_sim_matrix) {
                ecs_matrix[i, j] <- mean(ecs)
                ecs_matrix[j, i] <- ecs_matrix[i, j]
            }
        }
    }

    if (calculate_sim_matrix) {
        return(list(ecc = ecc, ecs_matrix = ecs_matrix))
    }

    ecc
}

test_that("weighted ECC matches pairwise reference for mixed flat labels", {
    clustering_list <- list(
        c(100, 100, 200, 200, 300, 300, 300, 400),
        factor(c("x", "x", "y", "z", "y", "z", "z", "w")),
        c("left", "right", "left", "right", "right", "left", "right", "left"),
        c(1, 1, 1, 2, 2, 3, 3, 3),
        c(-2, -2, 4, 4, 9, 9, 9, 12)
    )
    clustering_list <- lapply(clustering_list, function(x) {
        names(x) <- paste0("cell", seq_along(x))
        x
    })
    weights <- c(2, 3, 1, 4, 2)

    expected <- manual_weighted_element_consistency(clustering_list, weights)
    actual <- weighted_element_consistency(clustering_list, weights)

    expect_equal(actual, expected)
    expect_equal(names(actual), names(clustering_list[[1]]))
})

test_that("weighted ECC matrix output matches pairwise reference", {
    clustering_list <- list(
        c(10, 10, 20, 20, 30, 30),
        c("a", "b", "a", "b", "a", "b"),
        factor(c("u", "u", "u", "v", "v", "v")),
        c(100, 100, 200, 200, 300, 300),
        c("left", "right", "left", "right", "left", "right")
    )
    clustering_list <- lapply(clustering_list, function(x) {
        names(x) <- paste0("cell", seq_along(x))
        x
    })
    weights <- c(2, 1, 3, 2, 1)

    expected <- manual_weighted_element_consistency(
        clustering_list,
        weights,
        calculate_sim_matrix = TRUE
    )
    actual <- weighted_element_consistency(
        clustering_list,
        weights,
        calculate_sim_matrix = TRUE
    )

    expect_equal(actual$ecc, expected$ecc)
    expect_equal(actual$ecs_matrix, expected$ecs_matrix)
})

test_that("weighted ECC preserves mixed named and unnamed membership behavior", {
    unnamed <- c(1, 1, 2, 2)
    named_a <- setNames(c(1, 2, 1, 2), paste0("cell", 1:4))
    named_b <- setNames(c(1, 2, 1, 2), paste0("other", 1:4))
    good <- list(unnamed, named_a, named_a, unnamed, named_a)
    bad <- list(unnamed, named_a, named_a, unnamed, named_b)
    all_named_bad <- list(named_a, named_a, named_a, named_a, named_b)

    expect_equal(
        weighted_element_consistency(good, weights = rep(1, 5)),
        manual_weighted_element_consistency(good, weights = rep(1, 5))
    )
    expect_error(
        weighted_element_consistency(bad, weights = rep(1, 5)),
        "Not all elements"
    )
    expect_error(
        weighted_element_consistency(all_named_bad, weights = rep(1, 5)),
        "Not all elements"
    )
})

test_that("flat weighted ECC matches the pairwise reference across list sizes", {
    set.seed(11)
    element_names <- paste0("cell", 1:128)

    for (n_clusterings in c(2L, 3L, 5L, 20L)) {
        clustering_list <- lapply(seq_len(n_clusterings), function(i) {
            setNames(sample.int(7L, length(element_names), replace = TRUE), element_names)
        })
        weights <- sample.int(4L, n_clusterings, replace = TRUE)

        expect_equal(
            weighted_element_consistency(clustering_list, weights),
            manual_weighted_element_consistency(clustering_list, weights)
        )
    }
})

test_that("merge_partitions merges relabelled identical memberships", {
    x <- c(1, 1, 2, 2)
    relabelled_x <- c("a", "a", "b", "b")
    y <- c(1, 2, 1, 2)

    merged <- merge_partitions(
        list(x, relabelled_x, y),
        ecs_thresh = 1,
        order_logic = "freq"
    )

    expect_equal(length(merged$partitions), 2)
    expect_equal(sum(vapply(merged$partitions, `[[`, numeric(1), "freq")), 3)
    expect_equal(merged$partitions[[1]]$freq, 2)
    expect_equal(merged$partitions[[2]]$freq, 1)
    expect_equal(merged$ecc, rep(2 / 3, length(x)))
})

test_that("merge_partitions merges non-contiguous numeric relabellings", {
    x <- c(1L, 3L, 1L, 3L)
    relabelled_x <- c(10L, 20L, 10L, 20L)

    merged <- merge_partitions(
        list(x, relabelled_x),
        ecs_thresh = 1,
        order_logic = "none"
    )

    expect_equal(length(merged$partitions), 1)
    expect_equal(merged$partitions[[1]]$freq, 2)
    expect_equal(merged$ecc, rep(1, length(x)))
})

test_that("duplicate keys preserve numeric label precision", {
    x <- rep(c(1, 1 + 1e-15), 2)
    relabelled_x <- rep(1:2, 2)

    expect_true(ClustAssess:::are_identical_memberships(x, relabelled_x))
    expect_identical(
        ClustAssess:::canonical_membership_key(x),
        ClustAssess:::canonical_membership_key(relabelled_x)
    )

    merged <- merge_partitions(list(x, relabelled_x), order_logic = "none")

    expect_length(merged$partitions, 1)
    expect_equal(merged$partitions[[1]]$freq, 2)
    expect_equal(merged$ecc, rep(1, length(x)))
})

test_that("merge_partitions returns a symmetric ECS matrix for unique partitions", {
    x <- c(1, 1, 2, 2)
    y <- c(1, 2, 1, 2)
    z <- c(1, 1, 1, 2)

    merged <- merge_partitions(
        list(x, y, z),
        ecs_thresh = 1,
        order_logic = "none",
        return_ecs_matrix = TRUE
    )

    expect_equal(dim(merged$ecs_matrix), c(3L, 3L))
    expect_equal(merged$ecs_matrix, t(merged$ecs_matrix))
    expect_equal(diag(merged$ecs_matrix), rep(1, 3))
})

test_that("merge_partitions reorders the full ECS matrix after frequency tie-breaking", {
    a <- c(1L, 1L, 1L, 1L)
    b <- c(1L, 1L, 2L, 2L)
    c <- c(1L, 1L, 2L, 3L)
    merged <- merge_partitions(
        list(a, a, b, b, c),
        order_logic = "freq",
        check_ties = TRUE,
        return_ecs_matrix = TRUE
    )

    # B wins the frequency tie on average agreement. Its similarity to C
    # differs from A's, so permuting only the tied top-left block is wrong.
    returned_memberships <- lapply(merged$partitions, `[[`, "mb")
    expect_equal(returned_memberships, list(b, a, c))
    expect_equal(merged$ecs_matrix, element_sim_matrix(returned_memberships))
    expect_equal(merged$ecs_matrix[1:2, 3], c(0.75, 0.375))
    expect_equal(merged$ecs_matrix[3, 1:2], c(0.75, 0.375))
})

test_that("element_sim_matrix matches pairwise ECS for flat labels", {
    clustering_list <- list(
        c(10, 10, 30, 30, 50, 50),
        factor(c("left", "left", "right", "right", "right", "left")),
        c("a", "b", "a", "b", "b", "a")
    )

    expected <- matrix(NA_real_, nrow = length(clustering_list), ncol = length(clustering_list))
    for (i in seq_along(clustering_list)) {
        for (j in seq_along(clustering_list)) {
            expected[i, j] <- if (i == j) {
                1
            } else {
                element_sim(clustering_list[[i]], clustering_list[[j]])
            }
        }
    }

    expect_equal(element_sim_matrix(clustering_list), expected)
})

test_that("element_sim_matrix data.frame output preserves matrix order", {
    clustering_list <- list(
        c(1, 1, 2, 2),
        c(1, 2, 1, 2),
        c(1, 1, 1, 2)
    )

    ecs_matrix <- element_sim_matrix(clustering_list, output_type = "matrix")
    ecs_data <- element_sim_matrix(clustering_list, output_type = "data.frame")

    expect_equal(ecs_data$i.clustering, rep(seq_along(clustering_list), length(clustering_list)))
    expect_equal(ecs_data$j.clustering, rep(seq_along(clustering_list), each = length(clustering_list)))
    expect_equal(ecs_data$element_sim, c(ecs_matrix))
})

test_that("element_sim_matrix keeps pairwise input validation", {
    named_x <- c(a = 1, b = 1, c = 2)
    named_y <- c(a = 1, c = 1, b = 2)

    expect_error(
        element_sim_matrix(list(c(1, 1, 2), c(1, 2))),
        "same length"
    )
    expect_error(
        element_sim_matrix(list(named_x, named_y)),
        "same"
    )
})

test_that("merge_partitions avg_agreement uses merged partition frequencies", {
    x <- c(1, 1, 2, 2)
    relabelled_x <- c("x", "x", "y", "y")
    y <- c(1, 2, 1, 2)
    z <- c(1, 1, 1, 2)

    merged <- merge_partitions(
        list(x, relabelled_x, y, z),
        ecs_thresh = 1,
        order_logic = "avg_agreement",
        check_ties = FALSE
    )

    frequencies <- vapply(merged$partitions, `[[`, numeric(1), "freq")
    average_agreement <- vapply(merged$partitions, `[[`, numeric(1), "avg_agreement")

    expected_x_agreement <- (
        1 +
            element_sim(x, y) +
            element_sim(x, z)
    ) / 3

    expect_equal(frequencies[[1]], 2)
    expect_equal(merged$partitions[[1]]$mb, x)
    expect_equal(average_agreement[[1]], expected_x_agreement)
    expect_true(all(average_agreement[-1] <= average_agreement[[1]]))
})

test_that("merge_partitions freq ordering keeps ECC stable without requesting ECS matrix", {
    x <- c(1, 1, 2, 2)
    relabelled_x <- c("x", "x", "y", "y")
    y <- c(1, 2, 1, 2)
    z <- c(1, 1, 1, 2)

    no_matrix <- merge_partitions(
        list(x, relabelled_x, y, z),
        ecs_thresh = 1,
        order_logic = "freq",
        check_ties = FALSE,
        return_ecs_matrix = FALSE
    )
    with_matrix <- merge_partitions(
        list(x, relabelled_x, y, z),
        ecs_thresh = 1,
        order_logic = "freq",
        check_ties = FALSE,
        return_ecs_matrix = TRUE
    )

    expect_null(no_matrix$ecs_matrix)
    expect_false(any(vapply(no_matrix$partitions, function(x) {
        !is.null(x$avg_agreement)
    }, logical(1))))
    expect_true(all(vapply(with_matrix$partitions, function(x) {
        !is.null(x$avg_agreement)
    }, logical(1))))
    expect_equal(no_matrix$ecc, with_matrix$ecc)
    expect_equal(
        vapply(no_matrix$partitions, `[[`, numeric(1), "freq"),
        vapply(with_matrix$partitions, `[[`, numeric(1), "freq")
    )
})

test_that("merge_partitions preserves explicit frequencies when merging duplicates", {
    x <- c(10, 10, 20, 20)
    relabelled_x <- c("left", "left", "right", "right")
    y <- c(1, 2, 1, 2)
    partition_list <- list(
        list(mb = x, freq = 2),
        list(mb = relabelled_x, freq = 3),
        list(mb = y, freq = 4)
    )

    merged <- merge_partitions(
        partition_list,
        ecs_thresh = 1,
        order_logic = "freq",
        check_ties = FALSE
    )

    expect_equal(length(merged$partitions), 2)
    expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), c(5, 4))
    expect_equal(merged$partitions[[1]]$mb, x)
    expect_equal(merged$partitions[[2]]$mb, y)
})

test_that("merge_partitions handles large relabelled memberships without long-key failures", {
    x <- rep(seq_len(50), each = 300)
    relabelled_x <- paste0("cluster_", x)
    y <- rep(seq_len(50), length.out = length(x))

    merged <- merge_partitions(
        list(x, relabelled_x, y),
        ecs_thresh = 1,
        order_logic = "freq",
        check_ties = FALSE
    )

    expect_equal(length(merged$partitions), 2)
    expect_equal(vapply(merged$partitions, `[[`, numeric(1), "freq"), c(2, 1))
    expect_equal(merged$partitions[[1]]$mb, x)
})

test_that("merge_partitions forwards options for pre-merged and nested inputs", {
    x <- c(1, 1, 2, 2)
    y <- c(1, 2, 1, 2)
    z <- c(1, 1, 1, 2)

    premerged <- list(partitions = list(
        list(mb = x, freq = 2),
        list(mb = y, freq = 1),
        list(mb = z, freq = 1)
    ))
    direct <- merge_partitions(
        premerged$partitions,
        order_logic = "freq",
        check_ties = FALSE,
        return_ecs_matrix = TRUE
    )
    forwarded <- merge_partitions(
        premerged,
        order_logic = "freq",
        check_ties = FALSE,
        return_ecs_matrix = TRUE
    )

    expect_equal(forwarded$ecs_matrix, direct$ecs_matrix)
    expect_equal(forwarded$ecc, direct$ecc)

    nested <- list(
        first = premerged$partitions,
        second = list(
            list(mb = y, freq = 1),
            list(mb = z, freq = 1)
        )
    )
    nested_result <- merge_partitions(
        nested,
        order_logic = "none",
        return_ecs_matrix = TRUE
    )

    expect_named(nested_result, c("first", "second"))
    expect_true(!is.null(nested_result$first$ecs_matrix))
    expect_true(!is.null(nested_result$second$ecs_matrix))
})

test_that("merge_partitions ECS-threshold path is stable with a registered PSOCK backend", {
    skip_if_not_installed("doParallel")

    x <- setNames(c(1, 1, 1, 2, 2, 2), paste0("cell", 1:6))
    y <- setNames(c(1, 1, 2, 2, 2, 3), names(x))
    z <- setNames(c(1, 2, 2, 2, 3, 3), names(x))
    partition_list <- list(x, y, z)

    expected <- merge_partitions(
        partition_list,
        ecs_thresh = 0.5,
        order_logic = "none",
        return_ecs_matrix = TRUE
    )

    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
    }, add = TRUE)

    actual <- merge_partitions(
        partition_list,
        ecs_thresh = 0.5,
        order_logic = "none",
        return_ecs_matrix = TRUE
    )

    expect_equal(actual, expected)
})
