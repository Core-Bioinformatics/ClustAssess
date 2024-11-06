Sys.unsetenv("RETICULATE_PYTHON")
library(reticulate)
virtualenv_create("ClustAssess-env", packages = c("numpy", "clusim", "ClustAssessPy"))
use_python("~/.virtualenvs/ClustAssess-env/bin/python")
use_virtualenv("ClustAssess-env", required = TRUE)

test_that("ECS is between 0 and 1", {
    set.seed(1234)
    for (i in seq_len(100)) {
        x <- sample.int(10, size = 100, replace = TRUE)
        y <- sample.int(10, size = 100, replace = TRUE)
        ecs <- element_sim_elscore(x, y)
        max_ecs <- max(ecs)
        min_ecs <- min(ecs)

        expect_true(0 <= min_ecs && max_ecs <= 1)
    }
})

test_that("ECS produces expected output for known cases", {

})

test_that("ECS provides same results for different types of inputs (character, factor, numeric etc)", {
    set.seed(1234)
    for (i in seq_len(100)) {
        x <- sample.int(10, size = 100, replace = TRUE)
        y <- sample.int(10, size = 100, replace = TRUE)
        ecs1 <- element_sim_elscore(x, y)
        ecs2 <- element_sim_elscore(as.factor(x), as.factor(y))
        ecs3 <- element_sim_elscore(as.character(x), as.character(y))
        ecs4 <- element_sim_elscore(as.numeric(x), as.numeric(y))

        expect_equal(ecs1, ecs2)
        expect_equal(ecs1, ecs3)
        expect_equal(ecs1, ecs4)
    }
})

test_that("`element_sim` produces the average of `element_sim_elscore`", {
    set.seed(1234)
    for (i in seq_len(100)) {
        x <- sample.int(10, size = 100, replace = TRUE)
        y <- sample.int(10, size = 100, replace = TRUE)
        ecs1 <- element_sim_elscore(x, y)
        ecs2 <- element_sim(x, y)

        expect_equal(mean(ecs1), ecs2)
    }
})

test_that("ECS produces consistent results with clusim", {
    csim <- import("clusim")

    set.seed(1234)
    for (i in seq_len(100)) {
        x <- sample.int(10, size = 100, replace = TRUE)
        y <- sample.int(10, size = 100, replace = TRUE)
        ecs1 <- element_sim_elscore(x, y)

        cl1 <- csim$clustering$Clustering()
        cl1$from_membership_list(x)
        cl2 <- csim$clustering$Clustering()
        cl2$from_membership_list(y)
        ecs2 <- csim$sim$element_sim_elscore(cl1, cl2)[[1]]$tolist()

        expect_equal(ecs1, ecs2)
    }

})

test_that("ECS produces consistent results with ClustAssessPy", {
    np <- import("numpy")
    cpy <- import("ClustAssessPy")

    set.seed(1234)
    for (i in seq_len(100)) {
        x <- sample.int(10, size = 100, replace = TRUE)
        y <- sample.int(10, size = 100, replace = TRUE)
        ecs1 <- element_sim_elscore(x, y)
        ecs2 <- py_to_r(cpy$element_sim_elscore(np$array(x), np$array(y))$tolist())

        expect_equal(ecs1, ecs2)
    }

})

test_that("ECC gives the average of the pairwise ECS", {
    set.seed(1234)
    clustering_list <- lapply(seq_len(100), function(i) {
        sample.int(10, size = 100, replace = TRUE)
    })

    expected_ecc <- element_consistency(clustering_list)
    ecc <- rep(0, 100)
    n_comps <- 0
    for (i in seq_along(clustering_list)) {
        for (j in seq_along(clustering_list)) {
            if (i == j) {
                next
            }
            n_comps <- n_comps + 1
            ecs <- element_sim_elscore(clustering_list[[i]], clustering_list[[j]])
            ecc <- ecc + ecs
        }
    }

    ecc <- ecc / n_comps
    expect_equal(ecc, expected_ecc)
})

test_that("ECS matrix is calculated correctly", {
    set.seed(1234)
    clustering_list <- lapply(seq_len(100), function(i) {
        sample.int(10, size = 100, replace = TRUE)
    })

    suppressWarnings(ecs_matrix <- element_sim_matrix(clustering_list))
    expect_equal(ecs_matrix, t(ecs_matrix))

    for (i in seq_along(clustering_list)) {
        for (j in seq_along(clustering_list)) {
            if (i > j) {
                next
            }
            ecs <- element_sim_elscore(clustering_list[[i]], clustering_list[[j]])
            expect_equal(ecs_matrix[i, j], mean(ecs))
        }
    }
})

test_that("merging partitions keeps all partitions intact", {
    set.seed(1234)
    clustering_list <- lapply(seq_len(100), function(i) {
        sample.int(10, size = 100, replace = TRUE)
    })

    merged_clustering <- (merge_partitions(clustering_list))$partitions
    expect_equal(length(clustering_list), sum(sapply(merged_clustering, function(x) x$freq )))

    for (i in seq_along(merged_clustering)) {
        expect_true(!is.na(match(list(merged_clustering[[i]]$mb), clustering_list)))
    }
})

test_that("merging partitions provides correct results for ECC and ECS matrix", {
    set.seed(1234)
    clustering_list <- lapply(seq_len(100), function(i) {
        sample.int(10, size = 100, replace = TRUE)
    })

    merged_clustering <- merge_partitions(clustering_list, return_ecs_matrix = TRUE)
    partition_list <- lapply(merged_clustering$partitions, function(x) x$mb)

    ecc <- element_consistency(partition_list)
    ecs_matrix <- element_sim_matrix(partition_list)

    expect_equal(merged_clustering$ecc, ecc)
    expect_equal(merged_clustering$ecs_matrix, ecs_matrix)
})

test_that("merging partitions sorts the partitions correctly", {
    set.seed(1234)
    clustering_list <- lapply(seq_len(100), function(i) {
        sample.int(10, size = 100, replace = TRUE)
    })

    for (order_key in c("freq", "avg_agreement")) {
        merged_clustering <- (merge_partitions(clustering_list, order_logic = order_key))$partitions

        for (i in seq_along(merged_clustering)) {
            if (i == 1) {
                next
            }
            expect_true(merged_clustering[[i - 1]][[order_key]] >= merged_clustering[[i]][[order_key]])
        }
    }
})
