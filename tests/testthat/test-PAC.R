test_that("consensus_cluster is deterministic with full sampling on separated data", {
    x <- rbind(
        c(0, 0),
        c(0, 0.1),
        c(5, 5),
        c(5.1, 5),
        c(10, 10),
        c(10.1, 10)
    )

    set.seed(42)
    pac <- consensus_cluster(
        x,
        k_min = 2,
        k_max = 3,
        n_reps = 3,
        p_sample = 1,
        p_feature = 1,
        lower_lim = 0.1,
        upper_lim = 0.9
    )

    expect_named(
        pac,
        c(
            "dist_method", "linkage", "iteration", "pac", "lower_lim",
            "upper_lim", "n.clusters", "p_sample", "p_feature"
        )
    )
    expect_equal(nrow(pac), 6L)
    expect_equal(pac$iteration, rep(seq_len(3), each = 2))
    expect_equal(pac$n.clusters, rep(2:3, times = 3))
    expect_equal(pac$pac, rep(0, 6))
})

test_that("consensus_cluster caps k_max to sampled cell count", {
    x <- matrix(seq_len(12), ncol = 3)

    set.seed(1)
    pac <- consensus_cluster(
        x,
        k_min = 2,
        k_max = 10,
        n_reps = 2,
        p_sample = 1,
        p_feature = 1
    )

    expect_lte(max(pac$n.clusters), nrow(x))
})

test_that("PAC plotting helpers return ggplot objects", {
    testthat::skip_if_not_installed("ggplot2")

    x <- rbind(
        c(0, 0),
        c(0, 0.1),
        c(5, 5),
        c(5.1, 5),
        c(10, 10),
        c(10.1, 10)
    )

    set.seed(42)
    pac <- consensus_cluster(
        x,
        k_min = 2,
        k_max = 3,
        n_reps = 3,
        p_sample = 1,
        p_feature = 1
    )

    expect_s3_class(suppressWarnings(pac_convergence(pac, 2)), "ggplot")
    expect_s3_class(suppressWarnings(pac_landscape(pac)), "ggplot")
})
