###### UI ######

ui_graph_clustering_choice <- function(id) {
  ns <- shiny::NS(id)


  shiny::fluidRow(
    shiny::splitLayout(
      cellWidths = c("40px", "90%"),
      shinyWidgets::circleButton(ns("info_choice"),
        icon = shiny::icon("info"),
        size = "sm",
        status = "success"
      ),
      shiny::h2("Fixing the clustering method"),
    ),
    shiny::uiOutput(ns("clustering_options")),
    shiny::actionButton(ns("fix_cluster_button"),
      "Fix the method!",
      style = "font-size:20px;",
      class = "btn-danger"
    ),
    style = "padding:50px; font-size:20px;"
  )
}

#' Graph clustering - ui side
#'
#' @description to be completed
#'
#' @export
ui_graph_clustering <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Graph Clustering",
    shiny::uiOutput(ns("selected_conf")),
    shiny::uiOutput(ns("test_output")),
    shiny::uiOutput(ns("header1")),
    shiny::fluidRow(
      shiny::column(
        8,
        shiny::tabsetPanel(
          id = ns("boxplot_tabset"),
          shiny::tabPanel(
            "Vary by k",
            shiny::plotOutput(ns("boxplot_per_k"))
          ),
          shiny::tabPanel(
            "Vary by resolution",
            shiny::plotOutput(ns("boxplot_per_res"))
          ),
        )
      ),
      shiny::column(
        4,
        shiny::plotOutput(ns("umap_config"))
      )
    ),
    shiny::uiOutput(ns("header2")),
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::plotOutput(ns("boxplot_overall_resolution"))
      ),
      shiny::column(
        6,
        shiny::plotOutput(ns("boxplot_overall_k"))
      )
    ),
    shiny::uiOutput(ns("header3")),
    shiny::plotOutput(ns("k_resolution")),
    shiny::uiOutput(ns("header4")),
    shiny::plotOutput(ns("k_stability")),
    ui_graph_clustering_choice(ns("cluster_method_choice"))
  )
}

###### SERVER ######
server_graph_clustering_choice <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      clustering_options <- names(pkg_env$stab_obj$clustering_stability[[1]])
      output$clustering_options <- shiny::renderUI({
        ns <- session$ns
        shiny::radioButtons(
          inputId = ns("radio_cluster_method"),
          label = "Choose the clustering method for the downstream analysis:",
          choices = clustering_options,
          selected = clustering_options[1],
          width = "100%"
        )
      })

      # shiny::observe(dr_choice_info(session)) %>% shiny::bindEvent(input$info_choice)
      # shiny::observe({
      #   print(names(input))
      # }) %>% shiny::bindEvent(input)

      user_choice <- shiny::reactive(input$radio_cluster_method) %>% shiny::bindEvent(input$fix_cluster_button)

      # shiny::observe({
      #   print(paste("Schimbat user_choice la", user_choice()))
      # }) %>% shiny::bindEvent(user_choice())

      return(user_choice)
    }
  )
}

#' Graph clustering - server side
#'
#' @description to be completed
#'
#' @export
server_graph_clustering <- function(id, feature_choice) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      isolated_fchoice <- shiny::isolate(feature_choice())
      # print(isolated_fchoice)
      ftype <- isolated_fchoice$chosen_feature_type
      fsize <- isolated_fchoice$chosen_set_size

      output$test_output <- shiny::renderText(isolated_fchoice[[1]])
      
      # print(rhdf5::h5ls("stability.h5", recursive = 2))
      # print(paste(ftype, fsize, sep = "/"))

      temp_list <- rhdf5::h5read("stability.h5", paste(ftype, fsize, sep = "/"))
      temp_list$pca <- NULL
      temp_list$nn_stability <- NULL
      temp_list$nn_conn_comps <- NULL
      add_env_variable("stab_obj", temp_list)
      rm(temp_list)
      gc()

      # print(names(temp_list))
      

      stable_config <- pkg_env$stab_obj$stable_config
      # print(stable_config)
      # print(names(pkg_env))
      # print(names(pkg_env$stab_obj))

      # shiny::observe({
        output$selected_conf <- shiny::renderUI({
          shiny::fluidRow(
            shiny::div(glue::glue("Feature type: {stable_config[[1]]}")),
            div(glue::glue("Feature set size: {stable_config[[2]]}")),
            div(glue::glue("Base embedding: {stable_config[[3]]}")),
            div(glue::glue("Graph type: {stable_config[[4]]}")),
            div(glue::glue("Number of neighbours: {stable_config[[5]]}"))
          )
        })

        output$header1 <- shiny::renderUI({
          shiny::h3("Boxplot distribution of the resolution-wise stability of the clustering methods")
        })

        output$header2 <- shiny::renderUI({
          shiny::h3("Boxplot distribution of the overall stability")
        })

        output$header3 <- shiny::renderUI({
          shiny::h3("Correspondence between the resolution value and the number of clusters")
        })

        output$header4 <- shiny::renderUI({
          shiny::h3("The stability of the number of clusters")
        })
      # }) %>% shiny::bindEvent(feature_choice)

      output$boxplot_per_res <- shiny::renderPlot({
        shiny_plot_clustering_per_value_stability(pkg_env$stab_obj$clustering_stability,
          value_type = "resolution"
        )
      })

      output$boxplot_per_k <- shiny::renderPlot({
        shiny_plot_clustering_per_value_stability(pkg_env$stab_obj$clustering_stability,
          value_type = "k"
        )
      })

      output$boxplot_overall_resolution <- shiny::renderPlot({
        shiny_plot_clustering_overall_stability(pkg_env$stab_obj$clustering_stability,
          value_type = "resolution"
        )
      })

      output$boxplot_overall_k <- shiny::renderPlot({
        shiny_plot_clustering_overall_stability(pkg_env$stab_obj$clustering_stability,
          value_type = "k"
        )
      })

      output$k_resolution <- shiny::renderPlot({
        shiny_plot_k_resolution_corresp(pkg_env$stab_obj$clustering_stability,
          colour_information = "ecc"
        )
      })

      output$k_stability <- shiny::renderPlot({
        shiny_plot_k_n_partitions(pkg_env$stab_obj$clustering_stability,
          colour_information = "ecc"
        )
      })

      clustering_choice <- server_graph_clustering_choice("cluster_method_choice") #, clustering_importance)




      return(clustering_choice)

    }
  )
}

shiny_plot_clustering_per_value_stability <- function(clust_object,
                                                value_type = c("k", "resolution")) {
  value_type <- value_type[value_type %in% c("k", "resolution")]
  # TODO add empty boxplots for the missing values for k at least to help differentiate

  if (length(value_type) > 1) {
    value_type <- value_type[1]
  }

  if (length(value_type) == 0) {
    stop("`value_type` should contain either `k` or `resolution`")
  }

  clust_object <- clust_object[[paste0("split_by_", value_type)]]


  ecc_vals <- lapply(clust_object, function(by_alg) {
    lapply(by_alg, function(by_res_value) {
      as.numeric(by_res_value$ecc)
    })
  })

  melted_df <- reshape2::melt(ecc_vals)
  colnames(melted_df) <- c("ecc", value_type, "method")
  melted_df$method <- factor(melted_df$method)
  unique_vals <- stringr::str_sort(unique(melted_df[[value_type]]), numeric = TRUE)
  melted_df[[value_type]] <- factor(melted_df[[value_type]], levels = unique_vals)

  ggplot2::ggplot(
    melted_df,
    ggplot2::aes(x = .data[[value_type]], y = .data$ecc, fill = .data$method)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Clustering stability per ", value_type))
}

shiny_plot_clustering_overall_stability <- function(clust_object,
                                              value_type = c("k", "resolution"),
                                              summary_function = median) {
  value_type <- value_type[value_type %in% c("k", "resolution")]

  if (length(value_type) > 1) {
    value_type <- value_type[1]
  }

  if (length(value_type) == 0) {
    stop("`value_type` should contain either `k` or `resolution`")
  }

  ecc_vals <- lapply(
    clust_object[[paste0("split_by_", value_type)]],
    function(by_alg) {
      lapply(by_alg, function(by_value) {
        summary_function(by_value$ecc)
      })
    }
  )

  melted_df <- reshape2::melt(ecc_vals)
  colnames(melted_df) <- c("ecc", value_type, "method")
  melted_df$method <- factor(melted_df$method)
  melted_df[[value_type]] <- factor(melted_df[[value_type]])

  ggplot2::ggplot(
    melted_df,
    ggplot2::aes(x = .data$method, y = .data$ecc, fill = .data$method)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Overall clustering stability grouped by ", value_type))
}

shiny_plot_k_resolution_corresp <- function(clust_object,
                                      colour_information = c("ecc", "freq_k"),
                                      dodge_width = 0.3,
                                      pt_size_range = c(1.5, 4)) {
  # TODO check the colors and the vertical lines, try to help the user
  if (length(colour_information) > 1) {
    colour_information <- colour_information[1]
  }

  if (!(colour_information %in% c("ecc", "freq_k"))) {
    stop("colour_information can be either `ecc` or `freq_k`")
  }

  clust_object <- clust_object$split_by_resolution

  # use the names of the fields from the list
  res_object_names <- names(clust_object)

  # create a dataframe that contains the number of cases when,
  # for a given resolution, a number of clusters was obtained
  for (i in seq_along(clust_object)) {
    res_object <- clust_object[[i]]

    n_runs <- sum(sapply(res_object[[names(res_object)[1]]]$clusters, function(x) {
      x$total_freq
    }))

    list_appereances <- lapply(res_object, function(x) {
      lapply(x$clusters, function(y) {
        as.integer(y$first_freq)
      })
    })
    temp_appereances <- reshape2::melt(list_appereances)
    colnames(temp_appereances) <- c(
      "freq_partition",
      "number_clusters",
      "resolution_value"
    )

    temp_appereances[["freq_k"]] <- unlist(lapply(res_object, function(x) {
      lapply(x$clusters, function(y) {
        as.integer(y$total_freq)
      })
    }))

    temp_appereances[["configuration"]] <- rep(res_object_names[i], nrow(temp_appereances))
    temp_appereances$freq_partition <- temp_appereances$freq_partition / temp_appereances$freq_k
    temp_appereances$freq_k <- temp_appereances$freq_k / n_runs
    temp_appereances$ecc <- unlist(lapply(res_object, function(x) {
      sapply(x$clusters, function(k) {
        mean(k$ecc)
      })
    }))

    if (i == 1) {
      appearances_df <- temp_appereances
    } else {
      appearances_df <- rbind(appearances_df, temp_appereances)
    }
  }

  appearances_df[["configuration"]] <- factor(appearances_df[["configuration"]])
  appearances_df[["number_clusters"]] <- factor(as.numeric(appearances_df[["number_clusters"]]))

  main_plot <- ggplot2::ggplot(
    appearances_df,
    ggplot2::aes(
      y = .data$number_clusters,
      x = .data$resolution_value,
      size = .data$freq_partition,
      fill = .data[[colour_information]],
      shape = .data$configuration,
      group = .data$configuration
    )
  ) +
    ggplot2::geom_hline(
      yintercept = unique(appearances_df$number_clusters),
      linetype = "dashed",
      color = "#e3e3e3"
    ) +
    ggplot2::geom_vline(
      xintercept = unique(appearances_df$resolution_value),
      linetype = "dashed",
      color = "#e3e3e3"
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = dodge_width)) +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_viridis_c(guide = "colorbar") +
    ggplot2::scale_shape_manual(
      name = "Clustering method",
      values = 21:24,
      guide = "legend"
    ) +
    ggplot2::labs(
      x = "resolution",
      y = "k"
    ) +
    ggplot2::scale_size_continuous(range = pt_size_range, guide = "legend") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(size = max(pt_size_range)))
    )

  return(main_plot)
}

shiny_plot_k_n_partitions <- function(clust_object,
                                colour_information = c("ecc", "freq_part"),
                                dodge_width = 0.3,
                                pt_size_range = c(1.5, 4),
                                y_step = 5) {
  if (length(colour_information) > 1) {
    colour_information <- colour_information[1]
  }

  if (!(colour_information %in% c("ecc", "freq_part"))) {
    stop("colour_information can be either `ecc` or `freq_part`")
  }

  clust_object <- clust_object$split_by_k

  # use the names of the fields from the list
  object_names <- names(clust_object)

  max_n_part <- 0

  # creates a dataframe that contains, for each configuration
  # the number of different partitions with a given number of clusters
  for (i in seq_along(clust_object)) {
    partition_object <- clust_object[[i]]

    unique_parts_temp <- reshape2::melt(lapply(partition_object, function(x) {
      as.integer(x$n_partitions)
    }))
    colnames(unique_parts_temp) <- c("n.partitions", "n.clusters")

    unique_parts_temp[["configuration"]] <- rep(object_names[i], nrow(unique_parts_temp))

    unique_parts_temp[["first.occ"]] <- as.numeric(lapply(partition_object, function(x) {
      as.integer(x$first_freq)
    }))
    unique_parts_temp[["total.occ"]] <- as.numeric(lapply(partition_object, function(x) {
      as.integer(x$total_freq)
    }))
    unique_parts_temp[["freq_part"]] <- unique_parts_temp$first.occ / unique_parts_temp$total.occ
    unique_parts_temp[["ecc"]] <- sapply(partition_object, function(x) {
      mean(x$ecc)
    })
    overall_total_occ <- sum(unique_parts_temp$total.occ)
    unique_parts_temp[["frequency_k"]] <- unique_parts_temp$total.occ / overall_total_occ

    max_n_part <- max(c(max(
      unique_parts_temp$n.partitions
    ), max_n_part))

    if (i == 1) {
      unique_parts <- unique_parts_temp
    } else {
      unique_parts <- rbind(unique_parts, unique_parts_temp)
    }
  }

  unique_parts$configuration <- factor(unique_parts$configuration)
  unique_parts$n.clusters <- factor(unique_parts$n.clusters,
    levels = stringr::str_sort(unique(
      unique_parts$n.clusters
    ), numeric = TRUE)
  )

  main_plot <- ggplot2::ggplot(
    unique_parts,
    ggplot2::aes(
      x = .data$n.clusters,
      y = .data$n.partitions,
      shape = .data$configuration,
      size = .data$frequency_k,
      fill = .data[[colour_information]],
      group = .data$configuration
    )
  ) +
    ggplot2::scale_y_continuous(breaks = seq(
      from = 0,
      to = max_n_part,
      by = y_step
    )) +
    ggplot2::scale_size_continuous(range = pt_size_range, guide = "legend") +
    ggplot2::geom_hline(
      yintercept = seq(from = 0, to = max_n_part, by = y_step),
      linetype = "dashed",
      color = "#C3C3d3"

    ) +
    ggplot2::geom_vline(
      xintercept = unique(unique_parts$n.clusters),
      linetype = "dashed",
      color = "#c3c3c3d3"
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(dodge_width)) +
    ggplot2::theme_classic() +
    ggplot2::scale_shape_manual(name = "Clustering method", values = 21:24, guide = "legend") +
    ggplot2::scale_fill_viridis_c(
      name = ifelse(colour_information == "ecc",
        "ECC",
        "partition\nfrequency"
      ),
      guide = "colorbar"
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("k") +
    ggplot2::ylab("# partitions") +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(size = max(pt_size_range)))
    )


  return(main_plot)
}
