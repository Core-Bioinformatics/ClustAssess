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
            shiny::h2("Fixing the clustering method")
        ),
        shiny::radioButtons(
            inputId = ns("radio_cluster_method"),
            label = "Choose the clustering method for the downstream analysis:",
            choices = "",
            width = "100%"
        ),
        shinyWidgets::pickerInput(
            inputId = ns("select_n_clusters"),
            label = "Select the number of clusters for the downstream analysis",
            choices = "",
            inline = FALSE,
            options = list(
                `actions-box` = TRUE,
                title = "Select/deselect clusters",
                size = 10,
                width = "90%",
                `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
        ),
        shiny::actionButton(ns("fix_cluster_button"),
            "Fix the method!",
            style = "font-size:20px;",
            class = "btn-danger"
        ),
        style = "padding:50px; font-size:20px;"
    )
}

ui_graph_clustering_per_value_boxplot <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_boxplot_distro_res"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("Boxplot distribution of the resolution and ncluster-wise stability of the clustering methods"),
            class = "first-element-tab"
        ),
        shiny::tabsetPanel(
            id = ns("boxplot_tabset"),
            shiny::splitLayout(
                cellWidths = c("40px", "40px"),
                boxplot_settings(ns),
                gear_download(ns, "boxplot_vary", "boxplot_vary")
            ),
            shiny::tabPanel(
                "Vary by k",
                shiny::plotOutput(ns("boxplot_per_k"), height = "auto")
            ),
            shiny::tabPanel(
                "Vary by resolution",
                shiny::plotOutput(ns("boxplot_per_res"), height = "auto")
            )
        )
    )
}

ui_graph_clustering_per_value_umap <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        # shiny::uiOutput(ns("select_cluster_method_generator")),
        shiny::selectInput(
            inputId = ns("select_method"),
            label = "Select the clustering method",
            choices = ""
        ),
        # shiny::uiOutput(ns("select_n_clusters_generator")),
        shiny::selectInput(
            inputId = ns("select_k"),
            label = "Select the number of clusters (k)",
            choices = ""
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "300px"),
            gear_umaps(ns, "clustering_umap", TRUE, "lowest"),
            shinyWidgets::pickerInput(
                inputId = ns("select_clusters"),
                choices = "",
                inline = FALSE,
                # width = "100%",
                # width = "30%",
                options = list(
                    `actions-box` = TRUE,
                    title = "Select/deselect clusters",
                    # actionsBox = TRUE,
                    size = 10,
                    width = "90%",
                    `selected-text-format` = "count > 3"
                ),
                multiple = TRUE
            )
        ),
        shiny::fluidRow(
            # shiny::uiOutput(ns("umap_k_generator")),
            shiny::column(
                6,
                shiny::plotOutput(ns("umap_k"), height = "auto"),
                shiny::plotOutput(ns("umap_k_legend"), height = "auto"),
            ),
            shiny::column(
                6,
                shiny::plotOutput(ns("umap_ecc"), height = "auto"),
                shiny::plotOutput(ns("umap_ecc_legend"), height = "auto")
            )
            # shiny::uiOutput(ns("umap_ecc_generator"))
        )
    )
}

ui_graph_clustering_overall_boxplot <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_boxplot_distro_overall"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("Boxplot distribution of the overall stability"),
        ),
        shiny::splitLayout(
            gear_download(ns, "boxplot_overall_resolution", "boxplot_overall_resolution"),
            gear_download(ns, "boxplot_overall_k", "boxplot_overall_k")
        ),
        shiny::splitLayout(
            shiny::plotOutput(ns("boxplot_overall_resolution"), height = "auto"),
            shiny::plotOutput(ns("boxplot_overall_k"), height = "auto")
        )
    )
}

ui_graph_clustering_k_stab <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_k_corresp"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("Correspondence between the resolution value and the number of clusters")
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            shinyWidgets::dropdownButton(
                shiny::checkboxGroupInput(
                    inputId = ns("res_select_groups"),
                    label = "Select groups",
                    width = "100%",
                    choices = ""
                ),
                shiny::sliderInput(
                    inputId = ns("res_point_range"),
                    label = "Size point range",
                    min = 0.1, max = 5.00, value = c(1.35, 2.35)
                ),
                shiny::sliderInput(
                    inputId = ns("res_distance_factor"),
                    label = "Distance between groups",
                    min = 0.01, max = 1.00, value = 0.7
                ),
                shiny::sliderInput(
                    inputId = ns("res_text_size"),
                    label = "Text size",
                    min = 0.50, max = 10.00, value = 1.5
                ),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "k_resolution", "corresp_k_res")
        ),
        shiny::plotOutput(ns("k_resolution"), height = "auto"),
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_k_stab"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("The stability of the number of clusters")
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            shinyWidgets::dropdownButton(
                shiny::checkboxGroupInput(
                    inputId = ns("k_select_groups"),
                    label = "Select groups",
                    width = "100%",
                    choices = ""
                ),
                shiny::sliderInput(
                    inputId = ns("k_point_range"),
                    label = "Size point range",
                    min = 0.1, max = 5.00, value = c(1.35, 2.35)
                ),
                shiny::sliderInput(
                    inputId = ns("k_distance_factor"),
                    label = "Distance between groups",
                    min = 0.01, max = 1.00, value = 0.7
                ),
                shiny::sliderInput(
                    inputId = ns("k_ecc_thresh"),
                    label = "Low ECC threshold",
                    min = 0.00, max = 1.00, value = 0.90, step = 0.005
                ),
                shiny::sliderInput(
                    inputId = ns("k_occ_thresh"),
                    label = "Low frequency threshold",
                    min = 0.00, max = 1.00, value = 0.00
                ),
                shiny::sliderInput(
                    inputId = ns("k_nparts_thresh"),
                    label = "High # partitions threshold",
                    min = 0.00, max = 1.00, value = 0.00
                ),
                shiny::sliderInput(
                    inputId = ns("k_text_size"),
                    label = "Text size",
                    min = 0.50, max = 10.00, value = 1.5
                ),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "k_stability", "k_stability")
        ),
        shiny::plotOutput(ns("k_stability"), height = "auto"),
    )
}

#' UI - Graph clustering module
#'
#' @description Creates the UI interface for the graph clustering module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to identify the UI elements.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
ui_graph_clustering <- function(id) {
    ns <- shiny::NS(id)

    shiny::tabPanel(
        "Graph Clustering",
        shinyWidgets::circleButton(ns("info_title"),
            icon = shiny::icon("info"),
            size = "sm",
            status = "info",
            class = "page-info"
        ),
        shiny::actionButton(ns("show_config"), "Show config", type = "info", class = "btn-info show_config", disabled = TRUE),
        ui_graph_clustering_per_value_boxplot(ns("per_value_boxplot")),
        ui_graph_clustering_per_value_umap(ns("per_value_umap")),
        ui_graph_clustering_overall_boxplot(ns("overall_boxplot")),
        ui_graph_clustering_k_stab(ns("k_stab")),
        ui_graph_clustering_choice(ns("cluster_method_choice")),
    )
}

###### SERVER ######
server_graph_clustering_choice <- function(id, parent_session) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            clustering_options <- names(pkg_env$stab_obj$structure_list)
            shiny::updateRadioButtons(
                session = session,
                inputId = "radio_cluster_method",
                choices = clustering_options,
                selected = clustering_options[1],
            )

            shiny::observe({
                shiny::req(input$radio_cluster_method)
                k_values <- pkg_env$stab_obj$structure_list[[input$radio_cluster_method]]
                shinyWidgets::updatePickerInput(
                    session = session,
                    inputId = "select_n_clusters",
                    choices = k_values,
                    selected = k_values
                )
            }) %>% shiny::bindEvent(input$radio_cluster_method)


            user_choice <- shiny::reactive(
                list(
                    method_name = input$radio_cluster_method,
                    n_clusters = input$select_n_clusters
                )
            ) %>% shiny::bindEvent(input$fix_cluster_button)

            shiny::observe({
                shiny::req(input$select_n_clusters)
                shiny::showTab("tabset_id", "Comparison", select = TRUE, session = parent_session)
            }) %>% shiny::bindEvent(input$fix_cluster_button, ignoreInit = TRUE)

            shiny::observe(gclust_choice_info(session)) %>% shiny::bindEvent(input$info_choice, ignoreInit = TRUE)

            return(user_choice)
        }
    )
}

server_graph_clustering_per_value_boxplot <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            plt_height <- shiny::reactive(
                floor(pkg_env$height_ratio * pkg_env$dimension()[2])
            )

            plt_width <- shiny::reactive(
                pkg_env$dimension()[1]
            )

            output$boxplot_per_res <- shiny::renderPlot(
                height = function() {
                    plt_height()
                },
                {
                    shiny::req(pkg_env$current_tab == "Graph Clustering")
                    ecc_list <- pkg_env$stab_obj$ecc_by_res
                    split_names <- lapply(names(ecc_list), function(x) {
                        strsplit(x, ";")[[1]]
                    })
                    k_or_res_values <- as.numeric(sapply(split_names, function(x) {
                        x[1]
                    }))
                    cl_method <- sapply(split_names, function(x) {
                        x[2]
                    })

                    grouped_boxplot_list(
                        groups_list = ecc_list,
                        groups_values = cl_method,
                        x_values = k_or_res_values,
                        plt_height = plt_height(),
                        plt_width = plt_width(),
                        display_legend = TRUE,
                        xlab_text = "resolution",
                        boxplot_width = input$boxplot_width,
                        space_inter = input$space_inter_groups,
                        space_intra = input$space_intra_groups,
                        text_size = input$text_size
                    )
                }
            )

            output$boxplot_per_k <- shiny::renderPlot(
                height = function() {
                    plt_height()
                },
                {
                    # NOTE you need this check becuasee of the `suspendWhenHidden` option
                    shiny::req(pkg_env$current_tab == "Graph Clustering")
                    ecc_list <- pkg_env$stab_obj$ecc_by_k
                    split_names <- lapply(names(ecc_list), function(x) {
                        strsplit(x, ";")[[1]]
                    })
                    k_or_res_values <- as.numeric(sapply(split_names, function(x) {
                        x[1]
                    }))
                    cl_method <- sapply(split_names, function(x) {
                        x[2]
                    })


                    grouped_boxplot_list(
                        groups_list = ecc_list,
                        groups_values = cl_method,
                        x_values = k_or_res_values,
                        plt_height = plt_height(), #* 0.95,
                        plt_width = plt_width(),
                        display_legend = TRUE, # FIXME fix this (works atm)
                        xlab_text = "k",
                        boxplot_width = input$boxplot_width,
                        space_inter = input$space_inter_groups,
                        space_intra = input$space_intra_groups,
                        text_size = input$text_size
                    )
                }
            )

            shiny::outputOptions(output, "boxplot_per_k", suspendWhenHidden = FALSE)
            shiny::outputOptions(output, "boxplot_per_res", suspendWhenHidden = FALSE)

            output$download_boxplot_vary <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_boxplot_vary, ".", tolower(input$filetype_boxplot_vary))
                },
                content = function(file) {
                    shiny::req(pkg_env$current_tab == "Graph Clustering")
                    var_value <- ifelse(strsplit(input$boxplot_tabset, " ")[[1]][3] == "k", "k", "res")
                    ecc_list <- pkg_env$stab_obj[[paste0("ecc_by_", var_value)]]
                    split_names <- lapply(names(ecc_list), function(x) {
                        strsplit(x, ";")[[1]]
                    })
                    k_or_res_values <- as.numeric(sapply(split_names, function(x) {
                        x[1]
                    }))
                    cl_method <- sapply(split_names, function(x) {
                        x[2]
                    })


                    filetypes[[input$filetype_boxplot_vary]](file, width = input$width_boxplot_vary, height = input$height_boxplot_vary)

                    grouped_boxplot_list(
                        groups_list = ecc_list,
                        groups_values = cl_method,
                        x_values = k_or_res_values,
                        plt_height = input$height_boxplot_vary * ppi,
                        plt_width = input$width_boxplot_vary * ppi,
                        display_legend = TRUE, # FIXME fix this (works atm)
                        xlab_text = var_value,
                        boxplot_width = input$boxplot_width,
                        space_inter = input$space_inter_groups,
                        space_intra = input$space_intra_groups,
                        text_size = input$text_size
                    )


                    grDevices::dev.off()
                }
            )

            shiny::observe(gclust_distro_res_boxplot_info(session)) %>% shiny::bindEvent(input$info_boxplot_distro_res, ignoreInit = TRUE)
        }
    )
}

server_graph_clustering_overall_boxplot <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            plt_height <- shiny::reactive(
                floor(pkg_env$height_ratio * pkg_env$dimension()[2])
            )

            plt_width <- shiny::reactive(
                pkg_env$dimension()[1]
            )

            output$boxplot_overall_resolution <- shiny::renderPlot(
                height = function() {
                    floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2))
                },
                {
                    # shiny::req(pkg_env$lock_k())
                    # shiny_plot_clustering_overall_stability(pkg_env$stab_obj$clustering_stability,
                    #   value_type = "resolution"
                    # )
                    ftype <- pkg_env$lock_stable$feature_set
                    fsize <- pkg_env$lock_stable$n_features

                    grouped_boxplot_dataframe(
                        dataframe = pkg_env$stab_obj$summary_res, # rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_resolution", "summary", sep = "/")),
                        y_column = "ecc",
                        x_column = "cl_method",
                        plt_height = plt_height() - 1,
                        plt_width = plt_width(),
                        plot_title = "Stability by resolution"
                    )
                }
            )

            output$boxplot_overall_k <- shiny::renderPlot(
                height = function() {
                    floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2))
                },
                {
                    grouped_boxplot_dataframe(
                        dataframe = pkg_env$stab_obj$summary_k, # rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "summary", sep = "/")),
                        y_column = "ecc",
                        x_column = "cl_method",
                        plt_height = plt_height() - 1,
                        plt_width = plt_width(),
                        plot_title = "Stability by k"
                    )
                }
            )


            output$download_boxplot_overall_resolution <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_boxplot_overall_resolution, ".", tolower(input$filetype_boxplot_overall_resolution))
                },
                content = function(file) {
                    filetypes[[input$filetype_boxplot_overall_resolution]](file, width = input$width_boxplot_overall_resolution, height = input$height_boxplot_overall_resolution)

                    grouped_boxplot_dataframe(
                        dataframe = pkg_env$stab_obj$summary_res,
                        y_column = "ecc",
                        x_column = "cl_method",
                        plt_height = input$height_boxplot_overall_resolution * ppi,
                        plt_width = input$width_boxplot_overall_resolution * ppi,
                        plot_title = "Stability by k"
                    )
                    grDevices::dev.off()
                }
            )


            output$download_boxplot_overall_k <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_boxplot_overall_k, ".", tolower(input$filetype_boxplot_overall_k))
                },
                content = function(file) {
                    filetypes[[input$filetype_boxplot_overall_k]](file, width = input$width_boxplot_overall_k, height = input$height_boxplot_overall_k)

                    grouped_boxplot_dataframe(
                        dataframe = pkg_env$stab_obj$summary_k,
                        y_column = "ecc",
                        x_column = "cl_method",
                        plt_height = input$height_boxplot_overall_k * ppi,
                        plt_width = input$width_boxplot_overall_k * ppi,
                        plot_title = "Stability by k"
                    )
                    grDevices::dev.off()
                }
            )

            shiny::observe(gclust_distro_overall_boxplot_info(session)) %>% shiny::bindEvent(input$info_boxplot_distro_overall, ignoreInit = TRUE)
        }
    )
}

server_graph_clustering_per_value_umap <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            changed_method <- shiny::reactiveVal(FALSE)
            changed_k <- shiny::reactiveVal(FALSE)
            k_legend_height <- shiny::reactiveVal(0)
            ecc_legend_height <- shiny::reactiveVal(0)

            shiny::updateSelectInput(
                session = session,
                inputId = "select_method",
                choices = names(pkg_env$stab_obj$structure_list),
                selected = names(pkg_env$stab_obj$structure_list)[1],
            )

            shiny::observe({
                shiny::req(input$select_method != "")
                changed_method(TRUE)
                changed_k(TRUE)
                shiny::updateSelectInput(
                    session = session,
                    inputId = "select_k",
                    choices = pkg_env$stab_obj$structure_list[[input$select_method]],
                    selected = pkg_env$stab_obj$structure_list[[input$select_method]][1]
                )
            }) %>% shiny::bindEvent(input$select_method)

            shiny::observe({
                changed_method()
                input$select_k

                shiny::isolate({
                    shiny::req(input$select_k != "")
                    changed_k(TRUE)
                    shinyWidgets::updatePickerInput(
                        session = session,
                        inputId = "select_clusters",
                        choices = seq_len(as.integer(input$select_k)),
                        selected = seq_len(as.integer(input$select_k))
                    )
                })
            })

            plt_height <- shiny::reactive(
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.5))
            )



            output$umap_k <- shiny::renderPlot(
                height = function() {
                    plt_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    input$select_clusters
                    input$clustering_umap_labels
                    input$clustering_umap_text_size
                    input$clustering_umap_axis_size
                    input$clustering_umap_legend_size
                    input$clustering_umap_pt_size
                    input$clustering_umap_pt_type
                    input$select_method
                    input$select_k
                    plt_height()

                    shiny::isolate({
                        shiny::req(input$select_method != "", cancelOutput = TRUE)
                        if (changed_method()) {
                            shiny::req(input$select_k != "", input$select_k == pkg_env$stab_obj$structure_list[[input$select_method]][1], cancelOutput = TRUE)
                            changed_method(FALSE)
                        }
                        unique_values <- as.character(seq_len(as.integer(input$select_k)))


                        if (changed_k()) {
                            matched_elems <- match(input$select_clusters, unique_values)
                            matched_elems <- matched_elems[!is.na(matched_elems)]
                            shiny::req(
                                length(matched_elems) == length(unique_values),
                                length(matched_elems) == length(input$select_clusters),
                                cancelOutput = TRUE
                            )
                            changed_k(FALSE)
                        }

                        predicted_width <- graphics::strwidth(c(" ", unique_values), units = "inches", cex = input$clustering_umap_legend_size) * ppi
                        space_width <- predicted_width[1]
                        predicted_width <- predicted_width[2:length(predicted_width)]

                        number_columns <- min(
                            max(
                                plt_height() %/% (6 * space_width + max(predicted_width)),
                                1
                            ),
                            length(unique_values)
                        )
                        number_rows <- ceiling(length(unique_values) / number_columns)

                        text_height <- graphics::strheight(
                            paste(
                                rep("TEXT", number_rows + 1),
                                collapse = "\n"
                            ),
                            units = "inches",
                            cex = input$clustering_umap_legend_size
                        ) * ppi
                        k_legend_height(text_height)

                        color_plot2(
                            embedding = pkg_env$stab_obj$umap,
                            color_info = factor(rhdf5::h5read("stability.h5", paste(
                                pkg_env$lock_stable$feature_set,
                                pkg_env$lock_stable$n_features,
                                "clustering_stability",
                                "split_by_k",
                                "mbs",
                                input$select_method,
                                input$select_k,
                                sep = "/"
                            ))),
                            color_values = rhdf5::h5read("stability.h5", paste0("colors/", input$select_k)),
                            unique_values = seq_len(as.integer(input$select_k)),
                            plt_height = plt_height(),
                            plt_width = plt_height(),
                            pch = ifelse(input$clustering_umap_pt_type == "Pixel", ".", 19),
                            pt_size = input$clustering_umap_pt_size,
                            text_size = input$clustering_umap_text_size,
                            axis_size = input$clustering_umap_axis_size,
                            labels = input$clustering_umap_labels,
                            groups_highlight = input$select_clusters
                        )
                    })
                }
            )

            shiny::observe({
                shiny::req(k_legend_height() > 0, input$select_method != "")


                output$umap_k_legend <- shiny::renderPlot(
                    height = function() {
                        k_legend_height()
                    },
                    width = function() {
                        plt_height()
                    },
                    {
                        input$select_k
                        input$select_method
                        input$select_clusters
                        input$clustering_umap_legend_size
                        plt_height()

                        shiny::isolate({
                            unique_values <- as.character(seq_len(as.integer(input$select_k)))

                            if (!is.null(input$select_clusters)) {
                                unique_values <- as.integer(input$select_clusters)
                            }

                            only_legend_plot(
                                unique_values = unique_values,
                                color_values = rhdf5::h5read("stability.h5", paste0("colors/", input$select_k))[unique_values],
                                color_info = NULL,
                                plt_width = plt_height(),
                                text_size = input$clustering_umap_legend_size
                            )
                        })
                    }
                )
            })

            output$umap_ecc <- shiny::renderPlot(
                height = function() {
                    plt_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    input$select_clusters
                    input$clustering_umap_labels
                    input$clustering_umap_axis_size
                    input$clustering_umap_legend_size
                    input$clustering_umap_pt_size
                    input$clustering_umap_pt_order
                    input$clustering_umap_pt_type
                    input$select_method
                    input$select_k
                    plt_height()

                    shiny::isolate({
                        shiny::req(input$select_method != "", cancelOutput = TRUE)
                        if (changed_method()) {
                            shiny::req(input$select_k != "", input$select_k == pkg_env$stab_obj$structure_list[[input$select_method]][1], cancelOutput = TRUE)
                        }
                        unique_values <- as.character(seq_len(as.integer(input$select_k)))


                        if (changed_k()) {
                            matched_elems <- match(input$select_clusters, unique_values)
                            matched_elems <- matched_elems[!is.na(matched_elems)]
                            shiny::req(
                                length(matched_elems) == length(unique_values),
                                length(matched_elems) == length(input$select_clusters),
                                cancelOutput = TRUE
                            )
                        }
                        formatted_k <- sprintf("%06d", as.integer(input$select_k))

                        old_par <- graphics::par(mai = c(0.1, 0, 0.1, 0))
                        text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$clustering_umap_legend_size)
                        graphics::par(old_par)
                        ecc_legend_height(text_height * ppi)

                        color_plot2(
                            embedding = pkg_env$stab_obj$umap,
                            color_info = rhdf5::h5read("stability.h5", paste(
                                pkg_env$lock_stable$feature_set,
                                pkg_env$lock_stable$n_features,
                                "clustering_stability",
                                "split_by_k",
                                "ecc",
                                paste(formatted_k, input$select_method, sep = ";"),
                                sep = "/"
                            ))[
                                rhdf5::h5read("stability.h5", paste(
                                    pkg_env$lock_stable$feature_set,
                                    pkg_env$lock_stable$n_features,
                                    "clustering_stability",
                                    "split_by_k",
                                    "ecc_order",
                                    paste(formatted_k, input$select_method, sep = ";"),
                                    sep = "/"
                                ))
                            ],
                            color_values = NULL,
                            unique_values = NULL,
                            plt_height = plt_height(),
                            plt_width = plt_height(),
                            display_legend = FALSE,
                            sort_cells = input$clustering_umap_pt_order,
                            pch = ifelse(input$clustering_umap_pt_type == "Pixel", ".", 19),
                            pt_size = input$clustering_umap_pt_size,
                            legend_text_size = input$clustering_umap_legend_size,
                            axis_size = input$clustering_umap_axis_size
                        )
                    })
                }
            )

            shiny::observe({
                shiny::req(ecc_legend_height() > 0, input$select_method != "")


                output$umap_ecc_legend <- shiny::renderPlot(
                    height = function() {
                        ecc_legend_height()
                    },
                    width = function() {
                        plt_height()
                    },
                    {
                        input$select_k
                        input$select_method
                        input$select_clusters
                        input$clustering_umap_legend_size
                        plt_height()

                        shiny::isolate({
                            unique_values <- as.character(seq_len(as.integer(input$select_k)))

                            if (!is.null(input$select_clusters)) {
                                unique_values <- as.integer(input$select_clusters)
                            }
                            formatted_k <- sprintf("%06d", as.integer(input$select_k))

                            only_legend_plot(
                                unique_values = NULL,
                                color_values = NULL,
                                color_info = rhdf5::h5read("stability.h5", paste(
                                    pkg_env$lock_stable$feature_set,
                                    pkg_env$lock_stable$n_features,
                                    "clustering_stability",
                                    "split_by_k",
                                    "ecc",
                                    paste(formatted_k, input$select_method, sep = ";"),
                                    sep = "/"
                                ))[
                                    rhdf5::h5read("stability.h5", paste(
                                        pkg_env$lock_stable$feature_set,
                                        pkg_env$lock_stable$n_features,
                                        "clustering_stability",
                                        "split_by_k",
                                        "ecc_order",
                                        paste(formatted_k, input$select_method, sep = ";"),
                                        sep = "/"
                                    ))
                                ],
                                plt_width = plt_height(),
                                text_size = input$clustering_umap_legend_size
                            )
                        })
                    }
                )
            })
        }
    )
}

server_graph_clustering_k_stab <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            # render_plot_by_height("k_resolution", session)
            # render_plot_by_height("k_stability", session)

            shiny::updateSliderInput(
                session = session,
                inputId = "k_occ_thresh",
                max = max(pkg_env$stab_obj$summary_k$total_freq),
                min = min(pkg_env$stab_obj$summary_k$total_freq) - 1,
                step = 1
            )

            shiny::updateSliderInput(
                session = session,
                inputId = "k_nparts_thresh",
                max = max(pkg_env$stab_obj$summary_k$n_partitions) + 1,
                min = min(pkg_env$stab_obj$summary_k$n_partitions),
                value = max(pkg_env$stab_obj$summary_k$n_partitions) + 1,
                step = 1
            )

            shiny::updateCheckboxGroupInput(
                session = session,
                inputId = "k_select_groups",
                choices = unique(pkg_env$stab_obj$summary_k$cl_method),
                selected = unique(pkg_env$stab_obj$summary_k$cl_method)
            )

            shiny::updateCheckboxGroupInput(
                session = session,
                inputId = "res_select_groups",
                choices = unique(pkg_env$stab_obj$summary_k$cl_method),
                selected = unique(pkg_env$stab_obj$summary_k$cl_method)
            )

            plt_height <- shiny::reactive(
                floor(pkg_env$height_ratio * pkg_env$dimension()[2])
            )

            plt_width <- shiny::reactive(
                pkg_env$dimension()[1]
            )

            output$k_resolution <- shiny::renderPlot(
                {
                    # shiny::req(pkg_env$lock_k())
                    shiny_plot_k_res(
                        summary_df = pkg_env$stab_obj$summary_res,
                        distance_factor = input$res_distance_factor,
                        filtered_cl_methods = input$res_select_groups,
                        text_size = input$res_text_size,
                        pt_size_range = input$res_point_range
                    )
                },
                height = function() {
                    plt_height()
                }
            )

            output$k_stability <- shiny::renderPlot(
                {
                    # shiny::req(pkg_env$lock_k())
                    shiny::req(input$k_nparts_thresh > 0, input$k_select_groups)
                    shiny_plot_k_n_partitions(
                        summary_df = pkg_env$stab_obj$summary_k,
                        distance_factor = input$k_distance_factor,
                        filtered_cl_methods = input$k_select_groups,
                        text_size = input$k_text_size,
                        pt_size_range = input$k_point_range,
                        threshold_ecc = input$k_ecc_thresh,
                        threshold_occurences = input$k_occ_thresh,
                        threshold_nparts = input$k_nparts_thresh,
                        plt_height = plt_height(),
                        display_legend = TRUE
                    )
                },
                height = function() {
                    plt_height()
                }
            )

            shiny::observe(gclust_k_corresp_info(session)) %>% shiny::bindEvent(input$info_k_corresp, ignoreInit = TRUE)
            shiny::observe(gclust_k_stab_info(session)) %>% shiny::bindEvent(input$info_k_stab, ignoreInit = TRUE)
        }
    )
}

#' Server - Graph clustering module
#'
#' @description Creates the backend interface for the graph clustering module
#' inside the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to acess the UI elements.
#' @param feature_choice A reactive object that contains the chosen configuration
#' from the Dimensionality Reduction tab.
#' @param parent_session The session of the parent module, used to control the
#' tabs of the application.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
server_graph_clustering <- function(id, feature_choice, parent_session) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            isolated_fchoice <- shiny::isolate(feature_choice())
            ftype <- isolated_fchoice$chosen_feature_type
            fsize <- isolated_fchoice$chosen_set_size

            print(paste(Sys.time(), "Loading the stability object"))
            if (exists("stable_config")) {
                stable_config <- NULL
            }
            add_env_variable("stab_obj", list(
                ecc_by_k = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "ecc", sep = "/")),
                ecc_by_res = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_resolution", "ecc", sep = "/")),
                structure_list = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "structure_list", sep = "/")),
                summary_k = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "summary", sep = "/")),
                summary_res = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_resolution", "summary", sep = "/")),
                umap = rhdf5::h5read("stability.h5", paste(ftype, fsize, "umap", sep = "/"))
            ))
            stable_config <- rhdf5::h5read("stability.h5", paste(ftype, fsize, "stable_config", sep = "/"))
            add_env_variable("lock_choice", shiny::reactive(input$"cluster_method_choice-radio_cluster_method"))
            add_env_variable("lock_k", shiny::reactive(input$"per_value_umap-select_k"))
            add_env_variable("lock_stable", stable_config)
            add_env_variable("current_tab", "Graph Clustering")
            print(paste(Sys.time(), "Finished loading"))

            shiny::observe(
                shinyjs::enable("show_config")
            ) %>% shiny::bindEvent(input$"cluster_method_choice-radio_cluster_method", ignoreInit = TRUE, once = TRUE)

            clustering_choice <- server_graph_clustering_choice("cluster_method_choice", parent_session)

            server_graph_clustering_per_value_umap("per_value_umap")
            server_graph_clustering_per_value_boxplot("per_value_boxplot")
            server_graph_clustering_overall_boxplot("overall_boxplot")
            server_graph_clustering_k_stab("k_stab")

            shiny::observe({
                shiny::showModal(
                    stable_config_info(stable_config),
                    session = session
                )
            }) %>% shiny::bindEvent(input$show_config, ignoreInit = TRUE)

            shiny::observe(gclust_info(session)) %>% shiny::bindEvent(input$info_title, ignoreInit = TRUE)

            return(clustering_choice)
        }
    )
}

boxplot_settings <- function(ns) {
    shinyWidgets::dropdownButton(
        shiny::sliderInput(
            inputId = ns("boxplot_width"),
            label = "Boxplot width",
            min = 0.00, max = 1.00, value = 0.50
        ),
        shiny::sliderInput(
            inputId = ns("space_inter_groups"),
            label = "Space between groups",
            min = 1, max = 20, value = 1, step = 1
        ),
        shiny::sliderInput(
            inputId = ns("space_intra_groups"),
            label = "Space between boxplots inside group",
            min = 1, max = 20, value = 1, step = 1
        ),
        shiny::sliderInput(
            inputId = ns("text_size"),
            label = "Text size",
            min = 0.50, max = 10.00, value = 1
        ),
        circle = TRUE,
        status = "success",
        size = "sm",
        icon = shiny::icon("cog")
    )
}

stable_config_info <- function(stable_config) {
    shiny::modalDialog(
        shiny::HTML(paste0("<div><b>Feature type:</b> ", stable_config$feature_set, "</div>")),
        shiny::HTML(paste0("<div><b>Feature set size:</b> ", stable_config$n_features, "</div>")),
        shiny::HTML(paste0("<div><b>Number of PCs:</b> ", stable_config$n_pcs, "</div>")),
        shiny::HTML(paste0("<div><b>Graph base embedding:</b> ", stable_config$base_embedding, "</div>")),
        shiny::HTML(paste0("<div><b>Graph type:</b> ", stable_config$graph_type, "</div>")),
        shiny::HTML(paste0("<div><b>Graph pruning value:</b> ", round(stable_config$prune_param, 6), "</div>")),
        shiny::HTML(paste0("<div><b>Number of neighbours:</b> ", stable_config$n_neighbours, "</div>")),
        shiny::br(),
        shiny::em("Note: The stable configuration will be updated when changing the choices from the previous tabs."),
        title = "Current stable configuration",
        easyClose = TRUE
    )
}

##### GGPLOTS #####

shiny_ggplot_clustering_per_value_stability <- function(clust_object,
                                                        value_type = c("k", "resolution"),
                                                        width = 0.2,
                                                        dodge_width = 0.3,
                                                        text_size = 5) {
    print(Sys.time(), "Start shiny plot per value clust")
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
    lims <- as.numeric(c(unique_vals[1], unique_vals[length(unique_vals)]))
    lims <- lims - min(lims)
    melted_df[[value_type]] <- factor(melted_df[[value_type]], levels = unique_vals)

    print(Sys.time(), "Stop shiny plot per value clust")
    ggplot2::ggplot(
        melted_df,
        ggplot2::aes(x = .data[[value_type]], y = .data$ecc, fill = .data$method)
    ) +
        ggplot2::coord_cartesian(ylim = c(0, 1), xlim = lims) +
        ggplot2::geom_boxplot(
            position = ggplot2::position_dodge(width = dodge_width),
            width = width,
            outlier.shape = NA,
            outlier.size = 0.1
        ) +
        # ggplot2::geom_violin(
        #   position = ggplot2::position_dodge(width = dodge_width),
        #   width = width
        # ) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste0("Clustering stability per ", value_type)) +
        ggplot2::theme(
            legend.position = "bottom",
            text = ggplot2::element_text(size = text_size)
        )
}

shiny_ggplot_clustering_overall_stability <- function(clust_object,
                                                      value_type = c("k", "resolution"),
                                                      summary_function = stats::median) {
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
        ggplot2::ggtitle(paste0("Overall clustering stability grouped by ", value_type)) +
        ggplot2::theme(legend.position = "bottom")
}

##### PLOTS ######
shiny_plot_clustering_per_value_stability <- function(ecc_list,
                                                      value_type = c("k", "resolution"),
                                                      width = 0.5,
                                                      space_inter_groups = 1,
                                                      space_intra_groups = 1,
                                                      text_size = 1) {
    split_names <- lapply(names(ecc_list), function(x) {
        strsplit(x, ";")[[1]]
    })
    k_or_res_values <- sapply(split_names, function(x) {
        x[1]
    })
    unique_values <- as.numeric(unique(k_or_res_values))
    cl_method <- sapply(split_names, function(x) {
        x[2]
    })
    unique_cl_methods <- unique(cl_method)
    n_methods <- length(unique_cl_methods)
    cl_method <- match(cl_method, unique_cl_methods)

    at_values <- rep(0, length(cl_method))
    text_coords <- rep(0, length(unique_values))
    abline_coords <- rep(0, length(unique_values) - 1)

    count_diff <- -1
    prev_value <- -1
    updated_count <- FALSE
    start_index <- 1

    for (i in seq_along(cl_method)) {
        if (k_or_res_values[i] != prev_value) {
            if (prev_value != -1) {
                updated_count <- TRUE
            }
            prev_value <- k_or_res_values[i]
            count_diff <- count_diff + 1
            # abline_coords[count_diff+1] <- n_methods * space_intra + (n_methods * space_intra + space_inter) * count_diff + (space_inter) / 2
        }

        at_values[i] <- count_diff * (n_methods * space_intra_groups + space_inter_groups) + cl_method[i] + (cl_method[i] - 1) * (space_intra_groups - 1)

        if (updated_count) {
            abline_coords[count_diff] <- (at_values[i] + at_values[i - 1]) / 2
            updated_count <- FALSE
            text_coords[count_diff] <- mean(at_values[start_index:(i - 1)])
            start_index <- i
        }
    }

    text_coords[count_diff + 1] <- mean(at_values[start_index:length(cl_method)])


    graphics::boxplot(
        ecc_list,
        outline = FALSE,
        at = at_values,
        col = rhdf5::h5read("stability.h5", paste0("colors/", n_methods))[cl_method],
        boxwex = width * (n_methods + (space_intra_groups - 1) * (n_methods - 1)) / n_methods,
        xaxt = "n",
        xlab = NA,
        cex.axis = text_size,
        cex.lab = text_size
    )
    graphics::abline(v = abline_coords, lty = "dashed", col = "grey")
    graphics::title(xlab = "k", ylab = "ecc", cex.lab = text_size)
    graphics::axis(side = 1, at = text_coords, labels = unique_values, las = 2, cex.axis = text_size)
}

shiny_plot_clustering_per_value_stability_old <- function(clust_object,
                                                          value_type = c("k", "resolution"),
                                                          width = 0.5,
                                                          space_inter_groups = 1,
                                                          space_intra_groups = 1,
                                                          text_size = 5) {
    value_type <- value_type[value_type %in% c("k", "resolution")]

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
    colnames(melted_df) <- c("ecc", "x_value", "method")
    melted_df$method <- factor(melted_df$method)
    unique_vals <- stringr::str_sort(unique(melted_df[["x_value"]]), numeric = TRUE)
    lims <- as.numeric(c(unique_vals[1], unique_vals[length(unique_vals)]))
    lims <- lims - min(lims)
    melted_df[["x_value"]] <- factor(melted_df[["x_value"]], levels = unique_vals)

    method_k_df <- reshape2::melt(lapply(clust_object, function(by_alg) {
        names(by_alg)
    })) %>% dplyr::arrange(.data$value)
    method_k_df$L1 <- as.numeric(factor(method_k_df$L1))
    max_levels <- max(method_k_df$L1)
    offset <- floor(max_levels / 2) + max_levels %/% 2
    cluster_values <- unique(method_k_df$value)

    at_values <- c()
    name_values <- rep("", length(cluster_values) * max_levels)
    abline_coords <- rep(0, length(cluster_values) - 1)

    for (i in seq_along(cluster_values)) {
        name_values[(i - 1) * max_levels + offset] <- cluster_values[i]
        at_values <- c(at_values, seq_len(max_levels) * space_intra_groups + (max_levels * space_intra_groups + space_inter_groups) * (i - 1))

        if (i == 1) {
            next
        }

        abline_coords[i - 1] <- max_levels * space_intra_groups + (max_levels * space_intra_groups + space_inter_groups) * (i - 2) + (space_inter_groups + space_intra_groups) / 2
    }
    middle_values <- at_values[seq(from = offset, by = max_levels, to = length(name_values))]

    colorsi <- grDevices::rainbow(max_levels, s = 0.5)


    return({
        graphics::boxplot(
            ecc ~ method + x_value,
            data = melted_df,
            col = colorsi,
            pch = ".",
            outline = FALSE,
            at = at_values,
            xaxt = "n",
            boxwex = width * (max_levels + (space_intra_groups - 1) * (max_levels - 1)) / max_levels,
            xlab = NA,
            cex.axis = text_size,
            cex.lab = text_size
        )
        graphics::axis(side = 1, at = middle_values, labels = cluster_values, las = 2, cex.axis = text_size)
        graphics::title(xlab = value_type, cex.lab = text_size)
        graphics::abline(v = abline_coords, lty = "dashed", col = "grey")
    })
}

shiny_plot_clustering_overall_stability <- function(clust_object,
                                                    value_type = c("k", "resolution"),
                                                    summary_function = stats::median) {
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

    colorsi <- grDevices::rainbow(nlevels(melted_df$method), s = 0.5)
    return({
        graphics::boxplot(
            ecc ~ method,
            data = melted_df,
            col = colorsi,
            pch = "."
        )
    })
}

shiny_plot_k_res <- function(summary_df,
                             distance_factor = 0.6,
                             text_size = 1,
                             pt_size_range = c(0.1, 1.5),
                             filtered_cl_methods = NULL,
                             threshold_ecc = 0,
                             threshold_freq = 0) {
    available_pchs <- 21:24
    if (is.null(filtered_cl_methods)) {
        filtered_cl_methods <- c("Louvain", "Louvain.refined", "SLM", "Leiden")
    }

    summary_df$freq_partition <- summary_df$first_freq / summary_df$total_freq

    summary_df <- summary_df %>% dplyr::filter(
        .data$ecc > threshold_ecc &
            .data$cl_method %in% filtered_cl_methods &
            .data$freq_partition > threshold_freq
    )

    if (nrow(summary_df) == 0) {
        return()
    }

    cl_method <- summary_df$cl_method
    unique_cl_method <- unique(cl_method)
    cl_method <- match(cl_method, unique_cl_method)
    pchs <- available_pchs[factor(cl_method)]
    n_methods <- length(unique_cl_method)
    cl_method <- cl_method - mean(seq_len(n_methods))
    mask <- (cl_method < 0)
    cl_method[mask] <- floor(cl_method[mask])
    cl_method[!mask] <- ceiling(cl_method[!mask])
    max_distance <- ifelse(n_methods > 1, 0.9 / (n_methods - n_methods %% 2), 0)

    color_values <- viridis::viridis(min(50, nrow(summary_df)))
    color_info <- cut(summary_df$ecc, breaks = min(50, nrow(summary_df)))
    color_info <- color_values[color_info]

    x_axis <- summary_df$res
    unique_x_axis <- unique(x_axis)
    x_axis <- match(x_axis, unique_x_axis)
    x_axis <- x_axis + cl_method * max_distance * distance_factor

    actual_min <- min(summary_df$freq_partition)
    actual_max <- max(summary_df$freq_partition)

    if (actual_max > actual_min) {
        pt_sizes <- (pt_size_range[2] - pt_size_range[1]) * (summary_df$freq_partition - actual_min) / (actual_max - actual_min) + pt_size_range[1]
    } else {
        pt_sizes <- pt_size_range[2]
    }

    plot(
        x = x_axis,
        y = summary_df$k,
        pch = pchs,
        bg = color_info,
        xaxt = "n",
        yaxt = "n",
        xlab = "res",
        ylab = "k",
        cex = pt_sizes,
        cex.axis = text_size,
        cex.lab = text_size
    )
    graphics::abline(v = seq(from = 0.5, by = 1, length.out = length(unique_x_axis)), lty = "dashed", col = "grey")
    graphics::axis(side = 1, at = seq_along(unique_x_axis), labels = unique_x_axis, las = 2, cex.axis = text_size)
    graphics::axis(side = 2, at = unique(summary_df$k), cex.axis = text_size)
}

shiny_plot_k_n_partitions <- function(summary_df,
                                      plt_height,
                                      distance_factor = 0.6,
                                      text_size = 1,
                                      pt_size_range = c(0.1, 1.5),
                                      filtered_cl_methods = NULL,
                                      threshold_ecc = 0,
                                      threshold_occurences = 0,
                                      display_legend = FALSE,
                                      threshold_nparts = NULL) {
    available_pchs <- 21:24
    if (is.null(filtered_cl_methods)) {
        filtered_cl_methods <- c("Louvain", "Louvain.refined", "SLM", "Leiden")
    }

    if (is.null(threshold_nparts)) {
        threshold_nparts <- max(summary_df$n_partitions) + 1
    }

    summary_df <- summary_df %>% dplyr::filter(
        .data$ecc > threshold_ecc &
            .data$cl_method %in% filtered_cl_methods &
            .data$n_partitions < threshold_nparts &
            .data$total_freq > threshold_occurences
    )

    if (nrow(summary_df) == 0) {
        return()
    }

    cl_method <- summary_df$cl_method
    unique_cl_method <- unique(cl_method)
    cl_method <- match(cl_method, unique_cl_method)
    pchs <- available_pchs[factor(cl_method)]
    n_methods <- length(unique_cl_method)
    cl_method <- cl_method - mean(seq_len(n_methods))
    mask <- (cl_method < 0)
    cl_method[mask] <- floor(cl_method[mask])
    cl_method[!mask] <- ceiling(cl_method[!mask])
    max_distance <- ifelse(n_methods > 1, 0.9 / (n_methods - n_methods %% 2), 0)

    color_values <- viridis::viridis(min(50, nrow(summary_df)))

    if (max(summary_df$ecc) - min(summary_df$ecc) < 1e-3) {
        color_info <- factor(rep(summary_df$ecc[1], nrow(summary_df)))
    } else {
        color_info <- cut(summary_df$ecc, breaks = min(50, nrow(summary_df)), dig.lab = 3)
    }

    color_info <- color_values[color_info]

    x_axis <- summary_df$k
    unique_x_axis <- unique(x_axis)
    x_axis <- match(x_axis, unique_x_axis)
    x_axis <- x_axis + cl_method * max_distance * distance_factor
    actual_min <- min(summary_df$total_freq)
    actual_max <- max(summary_df$total_freq)

    if (actual_max > actual_min) {
        pt_sizes <- (pt_size_range[2] - pt_size_range[1]) * (summary_df$total_freq - actual_min) / (actual_max - actual_min) + pt_size_range[1]
    } else {
        pt_sizes <- pt_size_range[2]
    }

    if (display_legend) {
        plt_height <- plt_height / ppi

        predicted_height <- graphics::strheight(paste(rep("TEXT", 4), collapse = "\n"), units = "inches", cex = text_size) + graphics::strheight("TE\nXT", units = "inches", cex = pt_size_range[2] / 1.5)
        layout.matrix <- matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, ncol = 3, byrow = TRUE)

        graphics::layout(
            mat = layout.matrix,
            heights = c(
                graphics::lcm((plt_height - predicted_height - 0.1) * 2.54),
                graphics::lcm(predicted_height * 2.54)
            )
        )
    }

    plot(
        x = x_axis,
        y = summary_df$n_partitions,
        pch = pchs,
        bg = color_info,
        xaxt = "n",
        xlab = "k",
        ylab = "# different partitions",
        cex = pt_sizes,
        cex.axis = text_size,
        cex.lab = text_size
    )
    graphics::abline(v = seq(from = 0.5, by = 1, length.out = length(unique_x_axis)), lty = "dashed", col = "grey")
    graphics::axis(side = 1, at = seq_along(unique_x_axis), labels = unique_x_axis, las = 2, cex.axis = text_size)

    if (display_legend) {
        # point shape
        plot(seq_len(n_methods) - 1, rep(1, n_methods), axes = FALSE, bty = "n", ylab = "", xlab = "", cex = pt_size_range[2], pch = available_pchs[seq_len(n_methods)], main = "Clustering\nmethods", ylim = c(-1, 1.5), cex.main = text_size)
        graphics::text(y = -0.15, x = seq_len(n_methods) - 1, labels = unique_cl_method, cex = text_size, xpd = TRUE)

        # point size
        nfreqs <- sum(summary_df$total_freq)
        chosen_numbers <- seq(from = actual_min, to = actual_max, length.out = 5)
        if (actual_max > actual_min) {
            pt_sizes <- (pt_size_range[2] - pt_size_range[1]) * (chosen_numbers - actual_min) / (actual_max - actual_min) + pt_size_range[1]
        } else {
            pt_sizes <- rep(pt_size_range[2], 5)
        }
        plot(0:4, rep(1, 5), axes = FALSE, bty = "n", ylab = "", xlab = "", pch = 19, cex = pt_sizes, main = "k frequency", ylim = c(-1, 1.5), cex.main = text_size)
        graphics::text(y = -0.15, x = 0:4, labels = round(chosen_numbers / nfreqs, digits = 3), cex = text_size, xpd = TRUE)

        # point colour
        unique_values <- c(min(summary_df$ecc), max(summary_df$ecc))
        legend_image <- grDevices::as.raster(matrix(color_values, nrow = 1))
        plot(c(0, 1), c(-1, 1), type = "n", axes = F, bty = "n", ylab = "", xlab = "", main = "Average ECC", cex.main = text_size)
        graphics::text(y = -0.5, x = seq(0, 1, l = 5), labels = round(seq(from = unique_values[1], to = unique_values[2], length.out = 5), digits = 3), cex = text_size)
        graphics::rasterImage(legend_image, 0, 0, 1, 1)
    }
}




shiny_ggplot_k_resolution_corresp <- function(clust_object,
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
        ) +
        ggplot2::theme(legend.position = "bottom")

    return(main_plot)
}

shiny_ggplot_k_n_partitions <- function(clust_object,
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
        ) +
        ggplot2::theme(legend.position = "bottom")



    return(main_plot)
}
