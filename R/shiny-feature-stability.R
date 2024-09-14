proportion_widths <- 55

#### UI ####

ui_dimensionality_stability <- function(id) {
    ns <- shiny::NS(id)


    shiny::tagList(
        shiny::uiOutput(ns("stepchoosing")),
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_ecc_res"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("ECC per individual resolution values"),
            class = "first-element-tab"
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            shinyWidgets::dropdownButton(
                shinyWidgets::sliderTextInput(
                    inputId = ns("by_step_resolution"),
                    label = "Resolution",
                    choices = c("")
                ),
                gear_overall(ns, "by_step_res"),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "by_step_res", "by_step_res")
        ),
        shiny::splitLayout(
            cellWidths = paste0(c(proportion_widths - 3, 100 - proportion_widths), "%"),
            shiny::tagList(
                shiny::plotOutput(ns("boxplot_ecc"), height = "auto"),
                shiny::splitLayout(
                    cellWidths = c("40px", "90%"),
                    shinyWidgets::circleButton(ns("info_ecs_incr_res"),
                        icon = shiny::icon("info"),
                        size = "sm",
                        status = "success"
                    ),
                    shiny::h2("Incremental ECS per individual resolution values")
                ),
                shiny::splitLayout(
                    cellWidths = c("40px", "40px"),
                    shinyWidgets::dropdownButton(
                        shinyWidgets::sliderTextInput(
                            inputId = ns("incremental_resolution"),
                            label = "Resolution",
                            choices = c("")
                        ),
                        gear_overall(ns, "incremental_res"),
                        circle = TRUE,
                        status = "success",
                        size = "sm",
                        icon = shiny::icon("cog")
                    ),
                    gear_download(ns, "incremental_resolution", "incremental_resolution")
                ),
                shiny::plotOutput(ns("boxplot_incr"), height = "auto")
            ),
            shiny::div(
                style = "width: 90%; height: auto;",
                shiny::selectInput(
                    inputId = ns("select_ftype_umap"),
                    label = "Select feature type",
                    choices = ""
                ),
                shiny::selectInput(
                    inputId = ns("select_fsize_umap"),
                    label = "Select the feature set size",
                    choices = ""
                ),
                gear_umaps(ns, "ecc_res_umap", FALSE, "lowest"),
                shiny::plotOutput(ns("umap_ecc"), height = "auto", width = "100%"),
                shiny::plotOutput(ns("umap_ecc_legend"), height = "auto", width = "100%"),
                shiny::tableOutput(ns("table_ecc_info"))
            )
        ),
        shiny::hr(style = "border-top:3px solid;"),
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_ecc_overall"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("Overall stability")
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            shinyWidgets::dropdownButton(
                gear_overall(ns, "by_step"),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "by_step", "by_step")
        ),
        shiny::plotOutput(ns("overall_boxplot_ecc"), height = "auto"),
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_ecs_incremental_overall"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("Overall incremental stability")
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            shinyWidgets::dropdownButton(
                gear_overall(ns, "incremental"),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "incremental", "incremental")
        ),
        shiny::plotOutput(ns("overall_boxplot_incremental"), height = "auto")
    )
}

ui_dimensionality_distribution_plots <- function(id, draw_line) {
    ns <- shiny::NS(id)
    style <- ifelse(draw_line, "border-right:5px solid", "")

    shinyWidgets::panel(
        style = style,
        shiny::selectizeInput(ns("feature_type"), "Feature names", NULL),
        shiny::selectizeInput(ns("feature_steps"), "Feature set size", NULL),
        shiny::hr(),
        shiny::selectizeInput(
            inputId = ns("gene_expr"),
            choices = NULL,
            label = "Gene name(s)",
            width = "95%",
            multiple = TRUE,
            options = list(
                plugins = list("remove_button", "drag_drop")
            )
        ),
        shiny::splitLayout(
            shiny::numericInput(
                inputId = ns("expr_threshold"),
                label = "Gene expression threshold",
                min = 0, max = 10, value = 0, step = 0.01,
                width = "95%"
            ),
            shiny::numericInput(
                inputId = ns("relaxation"),
                label = "#genes not expressed",
                min = 0, max = 10, value = 0, step = 1,
                width = "95%"
            )
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            gear_umaps(ns, "gene", FALSE, "highest")
        ),
        shiny::plotOutput(ns("umap_gene"), height = "auto"),
        shiny::plotOutput(ns("umap_gene_legend"), height = "auto"),
        shiny::splitLayout(
            shiny::selectizeInput(
                inputId = ns("metadata"),
                label = "Metadata",
                choices = NULL
            ),
            shiny::verticalLayout(
                shiny::tags$b("Select the highlighted groups"),
                shinyWidgets::pickerInput(
                    inputId = ns("select_groups"),
                    choices = "",
                    inline = FALSE,
                    options = list(
                        `actions-box` = TRUE,
                        title = "Select/deselect groups",
                        size = 10,
                        width = "90%",
                        `selected-text-format` = "count > 3"
                    ),
                    multiple = TRUE
                )
            )
        ),
        gear_umaps(ns, "metadata", default_order = "highest"),
        shiny::plotOutput(ns("umap_metadata"), height = "auto"),
        shiny::plotOutput(ns("umap_metadata_legend"), height = "auto")
    )
}

ui_dimensionality_choice <- function(id) {
    ns <- shiny::NS(id)

    shiny::fluidRow(
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_choice"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("Fixing the feature configuration"),
        ),
        shiny::uiOutput(ns("radio_ft")),
        shiny::uiOutput(ns("radio_fs")),
        shiny::uiOutput(ns("choice")),
        shiny::actionButton(ns("fix_feature_button"),
            "Fix the configuration!",
            style = "font-size:20px;",
            class = "btn-danger"
        ),
        style = "padding:50px; font-size:20px;"
    )
}

#' UI - Dimensionality reduction module
#'
#' @description Creates the UI interface for the dimensionality reduction
#' module inside the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to identify the UI elements.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
ui_dimensionality_reduction <- function(id) {
    ns <- shiny::NS(id)

    shiny::tabPanel(
        "Dimensionality Reduction",
        shinyWidgets::circleButton(ns("info_title"),
            icon = shiny::icon("info"),
            size = "sm",
            status = "info",
            class = "page-info"
        ),
        ui_dimensionality_stability(ns("stability")),
        shiny::splitLayout(
            cellWidths = c("40px", "90%"),
            shinyWidgets::circleButton(ns("info_comparison"),
                icon = shiny::icon("info"),
                size = "sm",
                status = "success"
            ),
            shiny::h2("Pairwise comparison of gene and metadata distribution"),
        ),
        shiny::splitLayout(
            cellWidths = c("48%", "48%"),
            ui_dimensionality_distribution_plots(ns("distribution_left"), TRUE),
            ui_dimensionality_distribution_plots(ns("distribution_right"), FALSE),
        ),
        ui_dimensionality_choice(ns("feature_choice")),
        style = "margin-bottom:30px;"
    )
}

#### SERVER ####

server_dimensionality_stability <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            # TODO less urgent add option to zoom in on the ecc umap

            shiny::updateSelectInput(
                session = session,
                inputId = "select_ftype_umap",
                choices = names(pkg_env$feature_ordering$original)
            )

            shiny::observe({
                shiny::req(input$select_ftype_umap, input$select_ftype_umap != "")
                shiny::updateSelectInput(
                    session = session,
                    inputId = "select_fsize_umap",
                    choices = pkg_env$feature_ordering$original[[input$select_ftype_umap]]
                )
            })

            plt_height <- shiny::reactive({
                shiny::req(pkg_env$dimension())
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * (1 - proportion_widths / 100)))
            })

            legend_height <- shiny::reactive({
                grDevices::pdf(file = NULL, width = plt_height(), height = plt_height())
                graphics::par(mai = c(0.1, 0, 0.1, 0))
                text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$by_step_res_text_size)
                grDevices::dev.off()
                return((0.2 + text_height) * ppi)
            })

            ecc_value <- shiny::reactive({
                shiny::req(
                    input$select_fsize_umap,
                    input$select_fsize_umap != "",
                    input$by_step_resolution != "",
                    input$by_step_resolution != "0"
                )

                fsize_index <- which(input$select_fsize_umap == pkg_env$feature_ordering$original[[input$select_ftype_umap]])

                (rhdf5::h5read(
                    "stability.h5",
                    paste(
                        "feature_stability",
                        "by_steps",
                        input$by_step_resolution,
                        sep = "/"
                    )
                ) %>% dplyr::filter(.data$fsize == fsize_index & .data$ftype == input$select_ftype_umap))$ecc
            })

            output$umap_ecc <- shiny::renderPlot(
                height = function() {
                    shiny::req(plt_height())
                },
                width = function() {
                    shiny::req(plt_height())
                },
                {
                    shiny::req(
                        ecc_value(),
                        plt_height()
                    )

                    color_plot2(
                        embedding = rhdf5::h5read(
                            "stability.h5",
                            paste(
                                "feature_stability",
                                "embedding_list",
                                input$select_ftype_umap,
                                input$select_fsize_umap,
                                sep = "/"
                            )
                        ),
                        color_info = ecc_value(),
                        color_values = NULL,
                        plt_height = plt_height(),
                        plt_width = plt_height(),
                        axis_size = input$ecc_res_umap_axis_size,
                        pt_size = input$ecc_res_umap_pt_size,
                        sort_cells = input$ecc_res_umap_pt_order,
                        pch = ifelse(input$ecc_res_umap_pt_type == "Pixel", ".", 19),
                        legend_text_size = input$ecc_res_umap_legend_size,
                        text_size = input$ecc_res_umap_text_size
                    )
                }
            )

            output$umap_ecc_legend <- shiny::renderPlot(
                height = function() {
                    legend_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    shiny::req(ecc_value())

                    only_legend_plot(
                        plt_width = plt_height(),
                        text_size = input$ecc_res_umap_legend_size,
                        color_values = NULL,
                        unique_values = NULL,
                        color_info = ecc_value()
                    )
                }
            )

            output$table_ecc_info <- shiny::renderTable({
                shiny::req(ecc_value())

                data.frame(
                    "Stat" = c("Min", "Q1", "Median", "Q3", "Max"),
                    "ECC value" = as.character(trunc(stats::fivenum(ecc_value()) * 1e4) / 1e4)
                )
            })

            output$boxplot_ecc <- shiny::renderPlot(
                {
                    shiny::req(input$by_step_resolution != "", input$by_step_resolution != "0")

                    shiny_plot_feature_stability_boxplot(
                        resval = input$by_step_resolution,
                        feature_ordering = pkg_env$feature_ordering,
                        plt_height = floor(pkg_env$height_ratio * pkg_env$dimension()[2]),
                        plt_width = pkg_env$dimension()[1] * (proportion_widths - 5) / 100,
                        text_size = input$by_step_res_text_size,
                        width = input$by_step_res_boxplot_width,
                        space_intra_groups = input$by_step_res_intra_distance,
                        plt_title = paste("res", input$by_step_resolution),
                        space_inter_groups = input$by_step_res_inter_distance
                    )
                },
                height = function() {
                    floor(pkg_env$height_ratio * pkg_env$dimension()[2])
                }
            )

            output$boxplot_incr <- shiny::renderPlot(
                {
                    shiny::req(input$by_step_resolution != "", input$by_step_resolution != "0")

                    shiny_plot_feature_stability_incremental(
                        resval = input$incremental_resolution,
                        feature_ordering = pkg_env$feature_ordering,
                        plt_height = floor(pkg_env$height_ratio * pkg_env$dimension()[2]),
                        plt_width = pkg_env$dimension()[1] * (proportion_widths - 5) / 100,
                        text_size = input$incremental_res_text_size,
                        plt_title = paste("res", input$incremental_resolution),
                        width = input$incremental_res_boxplot_width,
                        space_intra_groups = input$incremental_res_intra_distance,
                        space_inter_groups = input$incremental_res_inter_distance
                    )
                },
                height = function() {
                    floor(pkg_env$height_ratio * pkg_env$dimension()[2])
                }
            )

            output$overall_boxplot_ecc <- shiny::renderPlot(
                {
                    shiny_plot_feature_stability_boxplot(
                        resval = "overall",
                        plt_height = floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2)),
                        plt_width = pkg_env$dimension()[1],
                        feature_ordering = pkg_env$feature_ordering,
                        text_size = input$by_step_text_size,
                        width = input$by_step_boxplot_width,
                        space_intra_groups = input$by_step_intra_distance,
                        space_inter_groups = input$by_step_inter_distance
                    )
                },
                height = function() {
                    floor(pkg_env$height_ratio * pkg_env$dimension()[2])
                }
            )

            output$overall_boxplot_incremental <- shiny::renderPlot(
                {
                    shiny_plot_feature_stability_incremental(
                        resval = "overall",
                        plt_height = floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2)),
                        plt_width = pkg_env$dimension()[1],
                        feature_ordering = pkg_env$feature_ordering,
                        text_size = input$incremental_text_size,
                        width = input$incremental_boxplot_width,
                        space_intra_groups = input$incremental_intra_distance,
                        space_inter_groups = input$incremental_inter_distance
                    )
                },
                height = function() {
                    floor(pkg_env$height_ratio * pkg_env$dimension()[2])
                }
            )

            output$download_by_step_res <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_by_step_res, ".", tolower(input$filetype_by_step_res))
                },
                content = function(file) {
                    shiny::req(input$by_step_resolution != "", input$by_step_resolution != "0")
                    filetypes[[input$filetype_by_step_res]](file, width = input$width_by_step_res, height = input$height_by_step_res)

                    shiny_plot_feature_stability_boxplot(
                        resval = input$by_step_resolution,
                        feature_ordering = pkg_env$feature_ordering,
                        plt_height = input$height_by_step_res * ppi,
                        plt_width = input$width_by_step_res * ppi,
                        text_size = input$by_step_res_text_size,
                        width = input$by_step_res_boxplot_width,
                        space_intra_groups = input$by_step_res_intra_distance,
                        plt_title = paste("res", input$by_step_resolution),
                        space_inter_groups = input$by_step_res_inter_distance
                    )
                    grDevices::dev.off()
                }
            )

            output$download_incremental_resolution <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_incremental_resolution, ".", tolower(input$filetype_incremental_resolution))
                },
                content = function(file) {
                    shiny::req(input$by_step_resolution != "", input$by_step_resolution != "0")
                    filetypes[[input$filetype_incremental_resolution]](file, width = input$width_incremental_resolution, height = input$height_incremental_resolution)

                    shiny_plot_feature_stability_incremental(
                        resval = input$incremental_resolution,
                        feature_ordering = pkg_env$feature_ordering,
                        plt_height = input$height_incremental_resolution * ppi,
                        plt_width = input$width_incremental_resolution * ppi,
                        text_size = input$incremental_res_text_size,
                        plt_title = paste("res", input$incremental_resolution),
                        width = input$incremental_res_boxplot_width,
                        space_intra_groups = input$incremental_res_intra_distance,
                        space_inter_groups = input$incremental_res_inter_distance
                    )
                    grDevices::dev.off()
                }
            )

            output$download_by_step <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_by_step, ".", tolower(input$filetype_by_step))
                },
                content = function(file) {
                    filetypes[[input$filetype_by_step]](file, width = input$width_by_step, height = input$height_by_step)

                    shiny_plot_feature_stability_boxplot(
                        resval = "overall",
                        feature_ordering = pkg_env$feature_ordering,
                        plt_height = input$height_by_step * ppi,
                        plt_width = input$width_by_step * ppi,
                        text_size = input$by_step_text_size,
                        width = input$by_step_boxplot_width,
                        space_intra_groups = input$by_step_intra_distance,
                        space_inter_groups = input$by_step_inter_distance
                    )
                    grDevices::dev.off()
                }
            )

            output$download_incremental <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_incremental, ".", tolower(input$filetype_incremental))
                },
                content = function(file) {
                    filetypes[[input$filetype_incremental]](file, width = input$width_incremental, height = input$height_incremental)

                    shiny_plot_feature_stability_incremental(
                        resval = "overall",
                        feature_ordering = pkg_env$feature_ordering,
                        plt_height = input$height_incremental * ppi,
                        plt_width = input$width_incremental * ppi,
                        text_size = input$incremental_text_size,
                        width = input$incremental_boxplot_width,
                        space_intra_groups = input$incremental_intra_distance,
                        space_inter_groups = input$incremental_inter_distance
                    )
                    grDevices::dev.off()
                }
            )

            shiny::observe(dr_individual_ecc_info(session)) %>% shiny::bindEvent(input$info_ecc_res, ignoreInit = TRUE)
            shiny::observe(dr_individual_incremental_info(session)) %>% shiny::bindEvent(input$info_ecs_incr_res, ignoreInit = TRUE)
            shiny::observe(dr_overall_ecc_info(session)) %>% shiny::bindEvent(input$info_ecc_overall, ignoreInit = TRUE)
            shiny::observe(dr_overall_incremental_info(session)) %>% shiny::bindEvent(input$info_ecs_incremental_overall, ignoreInit = TRUE)
        }
    )
}

server_dimensionality_distribution <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            changed_metadata <- shiny::reactiveVal(FALSE)
            metadata_legend_height <- shiny::reactiveVal(0)
            gene_legend_height <- shiny::reactiveVal(0)

            shiny::observeEvent(input$feature_type, {
                shiny::updateSelectizeInput(session,
                    inputId = "feature_steps",
                    choices = pkg_env$feature_ordering$stable[[input$feature_type]],
                    server = FALSE
                )
            })

            expr_matrix <- shiny::reactive({
                if ("genes" %in% names(pkg_env)) {
                    index_gene <- pkg_env$genes[input$gene_expr]
                    index_gene <- index_gene[!is.na(index_gene)] # not necesarry most probably

                    return(rhdf5::h5read("expression.h5", "expression_matrix", index = list(index_gene, NULL)))
                }

                # for backward-compatibility purposes
                index_interest <- pkg_env$genes_of_interest[input$gene_expr]
                index_interest <- index_interest[!is.na(index_interest)]

                index_others <- pkg_env$genes_others[input$gene_expr]
                index_others <- index_others[!is.na(index_others)]

                rbind(
                    rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(index_interest, NULL)),
                    rhdf5::h5read("expression.h5", "matrix_others", index = list(index_others, NULL))
                )
            }) %>% shiny::bindEvent(input$gene_expr)

            max_level_expr <- shiny::reactive(max(expr_matrix()))

            shiny::observe({
                shiny::updateSliderInput(session,
                    inputId = "expr_threshold",
                    max = round(max_level_expr(), 3),
                    step = round(max_level_expr() / 10, 3)
                )

                shiny::updateNumericInput(
                    session,
                    inputId = "relaxation",
                    max = length(input$gene_expr) - 1
                )

                if (length(input$gene_expr) > 1) {
                    shinyjs::show("relaxation")
                } else {
                    shinyjs::hide("relaxation")
                }
            }) %>% shiny::bindEvent(input$gene_expr)


            shiny::observe({
                mtd_names <- pkg_env$metadata_unique[[input$metadata]]
                if (is.null(mtd_names)) {
                    shinyjs::hide(id = "select_groups")
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = "select_groups",
                        choices = NULL,
                        selected = NULL
                    )
                } else {
                    shinyjs::show(id = "select_groups")
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = "select_groups",
                        choices = mtd_names,
                        selected = mtd_names
                    )
                }

                changed_metadata(TRUE)
            }) %>% shiny::bindEvent(input$metadata)

            plt_height <- shiny::reactive(
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.43))
            )

            output$umap_gene <- shiny::renderPlot(
                height = function() {
                    plt_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    shiny::req(input$gene_expr, input$feature_type, input$feature_steps, cancelOutput = TRUE)
                    plt_height()
                    relaxation <- input$relaxation
                    expr_threshold <- input$expr_threshold
                    input$gene_pt_type
                    input$gene_legend_size
                    input$gene_axis_size
                    input$gene_pt_size
                    input$gene_pt_order

                    shiny::isolate({
                        if (is.na(expr_threshold) || is.null(expr_threshold)) {
                            expr_threshold <- 0
                        }

                        if (is.na(relaxation) || is.null(relaxation)) {
                            relaxation <- 0
                        }

                        unique_values <- NULL
                        used_matrix <- expr_matrix()
                        color_values <- function(n) {
                            grDevices::colorRampPalette(c("grey85", RColorBrewer::brewer.pal(9, "OrRd")))(n)
                        }
                        if (length(input$gene_expr) > 1) {
                            unique_values <- c("other", "cells above threshold")
                            color_values <- c("FALSE" = "lightgray", "TRUE" = "red")
                            used_matrix <- matrixStats::colSums2(used_matrix > expr_threshold) >= (length(input$gene_expr) - relaxation)
                        } else if (expr_threshold > 0) {
                            unique_values <- c("other", "cells above threshold")
                            color_values <- c("FALSE" = "lightgray", "TRUE" = "red")
                            used_matrix <- used_matrix > expr_threshold
                        }

                        old_par <- graphics::par(mai = c(0.1, 0, 0.1, 0))
                        text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$gene_legend_size)
                        graphics::par(old_par)
                        gene_legend_height(text_height * ppi)

                        color_plot2(
                            embedding = rhdf5::h5read("stability.h5", paste(
                                "feature_stability",
                                "embedding_list",
                                input$feature_type,
                                input$feature_steps,
                                sep = "/"
                            )),
                            color_info = used_matrix,
                            plt_height = plt_height(),
                            plt_width = plt_height(),
                            display_legend = FALSE,
                            unique_values = unique_values,
                            color_values = color_values,
                            pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
                            pt_size = input$gene_pt_size,
                            axis_size = input$gene_axis_size,
                            sort_cells = input$gene_pt_order,
                            legend_text_size = input$gene_legend_size,
                            text_size = input$gene_legend_size
                        )
                    })
                }
            )

            output$umap_metadata <- shiny::renderPlot(
                height = function() {
                    plt_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    current_metadata <- input$metadata
                    groups <- input$select_groups
                    input$metadata_pt_size
                    input$metadata_axis_size
                    input$metadata_text_size
                    input$metadata_labels
                    input$metadata_pt_type
                    input$metadata_pt_order
                    input$metadata_legend_size
                    input$feature_type
                    input$feature_steps

                    shiny::isolate({
                        shiny::req(input$feature_steps, current_metadata, cancelOutput = TRUE)
                        if (changed_metadata() && !is.null(pkg_env$metadata_unique[[current_metadata]])) {
                            matched_elems <- match(groups, pkg_env$metadata_unique[[current_metadata]])
                            matched_elems <- matched_elems[!is.na(matched_elems)]
                            shiny::req(length(matched_elems) == length(pkg_env$metadata_unique[[current_metadata]]), cancelOutput = TRUE)
                            changed_metadata(FALSE)
                        }

                        if (is.null(pkg_env$metadata_unique[[current_metadata]])) {
                            old_par <- graphics::par(mai = c(0.1, 0, 0.1, 0))
                            text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$metadata_legend_size)
                            groups <- NULL
                        } else {
                            old_par <- graphics::par(mar = c(0, 0, 0, 0))
                            predicted_width <- graphics::strwidth(c(" ", pkg_env$metadata_unique[[current_metadata]]), units = "inches", cex = input$metadata_legend_size) * ppi
                            space_width <- predicted_width[1]
                            predicted_width <- predicted_width[2:length(predicted_width)]

                            number_columns <- min(
                                max(
                                    plt_height() %/% (6 * space_width + max(predicted_width)),
                                    1
                                ),
                                length(pkg_env$metadata_unique[[current_metadata]])
                            )
                            number_rows <- ceiling(length(pkg_env$metadata_unique[[current_metadata]]) / number_columns)

                            text_height <- graphics::strheight(
                                paste(
                                    rep("TEXT", number_rows + 1),
                                    collapse = "\n"
                                ),
                                units = "inches",
                                cex = input$metadata_legend_size
                            )
                        }
                        graphics::par(old_par)
                        metadata_legend_height(text_height * ppi)

                        metadata_plot(
                            embedding = rhdf5::h5read("stability.h5", paste(
                                "feature_stability",
                                "embedding_list",
                                input$feature_type,
                                input$feature_steps,
                                sep = "/"
                            )),
                            metadata_name = input$metadata,
                            plt_height = plt_height(),
                            plt_width = plt_height(),
                            pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
                            pt_size = input$metadata_pt_size,
                            sort_cells = input$metadata_pt_order,
                            text_size = input$metadata_text_size,
                            axis_size = input$metadata_axis_size,
                            labels = input$metadata_labels,
                            groups_highlight = groups
                        )
                    })
                }
            )

            shiny::observe({
                shiny::req(input$metadata, metadata_legend_height() > 0)
                output$umap_metadata_legend <- shiny::renderPlot(
                    height = function() {
                        metadata_legend_height()
                    },
                    width = function() {
                        plt_height()
                    },
                    {
                        plt_height()
                        input$select_groups
                        input$metadata_legend_size
                        current_metadata <- input$metadata

                        shiny::isolate({
                            if (!is.null(pkg_env$metadata_unique[[current_metadata]])) {
                                matched_elems <- match(input$select_groups, pkg_env$metadata_unique[[current_metadata]])
                                matched_elems <- matched_elems[!is.na(matched_elems)]
                                if (changed_metadata()) {
                                    shiny::req(
                                        length(matched_elems) == length(pkg_env$metadata_unique[[current_metadata]]),
                                        length(matched_elems) == length(input$select_groups),
                                        cancelOutput = TRUE
                                    )
                                }

                                # color_values <- pkg_env$metadata_colors[[current_metadata]][matched_elems]
                                unique_values <- pkg_env$metadata_unique[[current_metadata]][matched_elems]
                                color_values <- pkg_env$discrete_colors[[as.character(length(unique_values))]]
                            } else {
                                color_values <- NULL
                                unique_values <- NULL
                            }


                            only_legend_plot(
                                unique_values = unique_values,
                                color_values = color_values,
                                color_info = pkg_env$metadata[[current_metadata]],
                                plt_width = plt_height(),
                                text_size = input$metadata_legend_size
                            )
                        })
                    }
                )
            })

            shiny::observe({
                shiny::req(input$gene_expr, gene_legend_height() > 0)
                output$umap_gene_legend <- shiny::renderPlot(
                    height = function() {
                        gene_legend_height()
                    },
                    width = function() {
                        plt_height()
                    },
                    {
                        plt_height()
                        input$select_groups
                        expr_threshold <- input$expr_threshold
                        relaxation <- input$relaxation
                        input$gene_expr
                        input$gene_legend_size

                        shiny::isolate({
                            if (is.na(expr_threshold) || is.null(expr_threshold)) {
                                expr_threshold <- 0
                            }

                            if (is.na(relaxation) || is.null(relaxation)) {
                                relaxation <- 0
                            }

                            unique_values <- NULL
                            used_matrix <- expr_matrix()
                            color_values <- function(n) {
                                grDevices::colorRampPalette(c("grey85", RColorBrewer::brewer.pal(9, "OrRd")))(n)
                            }
                            if (length(input$gene_expr) > 1) {
                                unique_values <- c("other", "cells above threshold")
                                color_values <- c("#e3e3e3", "red")
                                used_matrix <- matrixStats::colSums2(used_matrix > expr_threshold) >= (length(input$gene_expr) - relaxation)
                            } else if (expr_threshold > 0) {
                                unique_values <- c("other", "cells above threshold")
                                color_values <- c("#e3e3e3", "red")
                                used_matrix <- used_matrix > expr_threshold
                            }


                            only_legend_plot(
                                unique_values = unique_values,
                                color_values = color_values,
                                color_info = used_matrix,
                                plt_width = plt_height(),
                                text_size = input$gene_legend_size
                            )
                        })
                    }
                )
            })
        }
    )
}

server_dimensionality_choice <- function(id, parent_session) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            output$radio_ft <- shiny::renderUI({
                ns <- session$ns
                shiny::radioButtons(
                    ns("radio_feature_type"),
                    label = "Choose the feature type for the downstream analysis:",
                    choices = names(pkg_env$feature_ordering$original),
                    selected = names(pkg_env$feature_ordering$original)[1],
                    width = "100%",
                )
            })

            output$radio_fs <- shiny::renderUI({
                ns <- session$ns
                shiny::radioButtons(
                    ns("radio_feature_size"),
                    label = "Choose the feature type for the downstream analysis:",
                    choices = pkg_env$feature_ordering$stable[[1]],
                    selected = pkg_env$feature_ordering$stable[[1]][1],
                    width = "100%",
                )
            })

            shiny::observeEvent(input$radio_feature_type, {
                ftype <- input$radio_feature_type
                shiny::updateRadioButtons(
                    session,
                    label = glue::glue("Choose the size of the feature set {ftype} for the downstream analysis: (We recommend {pkg_env$feature_ordering$stable[[ftype]][1]})"),
                    inputId = "radio_feature_size",
                    choices = pkg_env$feature_ordering$stable[[ftype]],
                    selected = pkg_env$feature_ordering$stable[[ftype]][1]
                )
            })

            shiny::observe(dr_choice_info(session)) %>% shiny::bindEvent(input$info_choice, ignoreInit = TRUE)

            user_choice <- shiny::reactive(list(
                chosen_feature_type = input$radio_feature_type,
                chosen_set_size = input$radio_feature_size
            )) %>% shiny::bindEvent(input$fix_feature_button)

            shiny::observe({
                shiny::showTab("tabset_id", "Graph Construction", select = FALSE, session = parent_session)
                shiny::showTab("tabset_id", "Graph Clustering", select = TRUE, session = parent_session)
            }) %>% shiny::bindEvent(input$fix_feature_button, ignoreInit = TRUE)

            return(user_choice)
        }
    )
}

#' Server - Dimensionality reduction module
#'
#' @description Creates the backend interface for the dimensionality
#' reduction module inside the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to acess the UI elements.
#' @param parent_session The session of the parent module, used to control the
#' tabs of the application.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
server_dimensionality_reduction <- function(id, parent_session) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            update_sliders(session)

            print(paste(Sys.time(), "Loading the stability object"))
            print(paste(Sys.time(), "Finished loading"))

            shiny::observe(dr_title_info(session)) %>% shiny::bindEvent(input$info_title, ignoreInit = TRUE)
            shiny::observe(dr_comparison_info(session)) %>% shiny::bindEvent(input$info_comparison, ignoreInit = TRUE)
            server_dimensionality_stability("stability")
            server_dimensionality_distribution("distribution_left")
            server_dimensionality_distribution("distribution_right")
            feature_choice <- server_dimensionality_choice("feature_choice", parent_session)

            feature_choice
        }
    )
}

update_sliders <- function(session) {
    shinyWidgets::updateSliderTextInput(
        session,
        "stability-by_step_resolution",
        choices = pkg_env$feature_ordering$resolution,
        selected = pkg_env$feature_ordering$resolution[1]
    )

    shinyWidgets::updateSliderTextInput(
        session,
        "stability-incremental_resolution",
        choices = pkg_env$feature_ordering$resolution,
        selected = pkg_env$feature_ordering$resolution[1]
    )

    for (panels in c("left", "right")) {
        shiny::updateSelectizeInput(
            session,
            inputId = glue::glue("distribution_{panels}-feature_type"),
            choices = names(pkg_env$feature_ordering$stable),
            selected = names(pkg_env$feature_ordering$stable)[1],
            server = FALSE
        )

        if ("genes" %in% names(pkg_env)) {
            gene_choices <- names(pkg_env$genes)
        } else { # for backward-compatibility purposes
            gene_choices <- c(names(pkg_env$genes_of_interest), names(pkg_env$genes_others))
        }

        shiny::updateSelectizeInput(
            session,
            inputId = glue::glue("distribution_{panels}-gene_expr"),
            choices = gene_choices,
            selected = gene_choices[1],
            server = TRUE,
            options = list(
                maxOptions = 7,
                create = TRUE,
                persist = TRUE
            )
        )

        shiny::updateSelectizeInput(
            session,
            inputId = glue::glue("distribution_{panels}-metadata"),
            server = FALSE,
            choices = colnames(pkg_env$metadata),
            selected = colnames(pkg_env$metadata)[1]
        )
    }
}

##### PLOTS #####

shiny_plot_feature_stability_boxplot <- function(resval,
                                                 feature_ordering,
                                                 plt_height,
                                                 plt_width,
                                                 plt_title = "",
                                                 text_size = 1,
                                                 width = 0.5,
                                                 space_intra_groups = 1,
                                                 space_inter_groups = 1,
                                                 text_offset = 0.01) {
    # convert pixels to inches
    plt_height <- plt_height / 96
    plt_width <- plt_width / 96
    # calculate space needed for the legend
    fgroups <- names(feature_ordering$original)
    predicted_width <- graphics::strwidth(fgroups, units = "inches", cex = text_size)
    predicted_height <- graphics::strheight(fgroups[1], units = "inches", cex = text_size)
    number_columns <- plt_width %/% (1.5 * max(predicted_width))
    number_rows <- ceiling(length(fgroups) / number_columns)
    current_margins <- graphics::par("mai")
    old_margin <- current_margins[1]
    current_margins[1] <- current_margins[1] + (number_rows + 2) * predicted_height[1] * 1.2
    if (text_size > 1.5) {
        current_margins[2] <- current_margins[2] + 0.5 * predicted_height[1]
    }
    graphics::par(mai = current_margins)

    # col <- rhdf5::h5read("stability.h5", "feature_stability/colours")
    n_groups <- length(fgroups)
    col <- pkg_env$discrete_colors[[as.character(n_groups)]]
    n_fsizes <- max(sapply(feature_ordering$original, function(x) {
        length(x)
    }))

    at_values <- c()
    name_values <- rep("", n_fsizes * n_groups)
    abline_coords <- rep(0, n_groups - 1)

    for (i in seq_len(n_fsizes)) {
        for (j in seq_along(fgroups)) {
            if (length(feature_ordering$original[[j]]) >= i) {
                name_values[(i - 1) * n_groups + j] <- feature_ordering$original[[j]][i]
            }
        }
        at_values <- c(at_values, seq_len(n_groups) * space_intra_groups + (n_groups * space_intra_groups + space_inter_groups) * (i - 1))

        if (i == 1) {
            next
        }

        abline_coords[i - 1] <- n_groups * space_intra_groups + (n_groups * space_intra_groups + space_inter_groups) * (i - 2) + (space_inter_groups + space_intra_groups) / 2
    }

    boundaries <- graphics::boxplot(
        formula = ecc ~ ftype + fsize,
        data = rhdf5::h5read("stability.h5", ifelse(resval == "overall",
            "feature_stability/overall/by_step",
            paste0("feature_stability/by_steps/", resval)
        )),
        col = col,
        at = at_values,
        xaxt = "n",
        xlab = ifelse(text_size > 1.5, "", "# features"),
        ylab = "ECC",
        boxwex = width * (n_groups + (space_intra_groups - 1) * (n_groups - 1)) / n_groups,
        outline = FALSE,
        cex.lab = text_size,
        cex.main = text_size * 1.2,
        cex.axis = text_size,
        main = plt_title
    )

    graphics::abline(v = abline_coords, lty = "dashed", col = "grey")

    y_text <- boundaries$stats[nrow(boundaries$stats), ]
    y_text[is.na(y_text)] <- 0

    graphics::axis(
        side = 1,
        at = at_values,
        labels = name_values,
        cex = text_size,
        cex.axis = text_size,
        las = 2,
        xpd = NA
    )
    
    graphics::legend(
        "bottomleft",
        legend = fgroups,
        col = col,
        pch = 15,
        cex = text_size,
        pt.cex = text_size * 2,
        bty = "n",
        ncol = number_columns,
        xpd = TRUE
    )
}

shiny_plot_feature_stability_incremental <- function(resval,
                                                     feature_ordering,
                                                     plt_height,
                                                     plt_width,
                                                     plt_title = "",
                                                     text_size = 1,
                                                     width = 0.5,
                                                     space_intra_groups = 1,
                                                     space_inter_groups = 1,
                                                     text_offset = 0.1) {
    # convert pixels to inches
    plt_height <- plt_height / 96
    plt_width <- plt_width / 96
    # calculate space needed for the legend
    fgroups <- names(feature_ordering$original_incremental)
    predicted_width <- graphics::strwidth(fgroups, units = "inches", cex = text_size)
    predicted_height <- graphics::strheight(fgroups[1], units = "inches", cex = text_size)
    number_columns <- plt_width %/% (1.2 * max(predicted_width))
    number_rows <- ceiling(length(fgroups) / number_columns)
    current_margins <- graphics::par("mai")
    old_margin <- current_margins[1]
    current_margins[1] <- current_margins[1] + (number_rows + 2) * predicted_height[1] * 1.2
    if (text_size > 1.5) {
        current_margins[2] <- current_margins[2] + 0.5 * predicted_height[1]
    }
    graphics::par(mai = current_margins)

    # col <- rhdf5::h5read("stability.h5", "feature_stability/colours")
    n_groups <- length(fgroups)
    col <- pkg_env$discrete_colors[[as.character(n_groups)]]
    n_fsizes <- max(sapply(feature_ordering$original_incremental, function(x) {
        length(x)
    }))

    at_values <- c()
    name_values <- rep("", n_fsizes * n_groups)
    abline_coords <- rep(0, n_groups - 1)

    for (i in seq_len(n_fsizes)) {
        for (j in seq_along(fgroups)) {
            if (length(feature_ordering$original_incremental[[j]]) >= i) {
                name_values[(i - 1) * n_groups + j] <- feature_ordering$original_incremental[[j]][i]
            }
        }
        at_values <- c(at_values, seq_len(n_groups) * space_intra_groups + (n_groups * space_intra_groups + space_inter_groups) * (i - 1))

        if (i == 1) {
            next
        }

        abline_coords[i - 1] <- n_groups * space_intra_groups + (n_groups * space_intra_groups + space_inter_groups) * (i - 2) + (space_inter_groups + space_intra_groups) / 2
    }

    boundaries <- graphics::boxplot(
        formula = ecc ~ ftype + fsize,
        data = rhdf5::h5read("stability.h5", ifelse(resval == "overall",
            "feature_stability/overall/incremental",
            paste0("feature_stability/incremental/", resval)
        )),
        col = col,
        at = at_values,
        xaxt = "n",
        xlab = "",
        ylab = "ECS",
        boxwex = width * (n_groups + (space_intra_groups - 1) * (n_groups - 1)) / n_groups,
        outline = FALSE,
        cex.lab = text_size,
        cex.main = text_size * 1.2,
        cex.axis = text_size,
        main = plt_title
    )

    graphics::abline(v = abline_coords, lty = "dashed", col = "grey")

    y_text <- boundaries$stats[nrow(boundaries$stats), ]
    y_text[is.na(y_text)] <- 0

    graphics::axis(
        side = 1,
        at = at_values,
        labels = name_values,
        cex = text_size,
        cex.axis = text_size,
        las = 2
    )

    graphics::legend(
        "bottomleft",
        legend = fgroups,
        col = col,
        pch = 15,
        cex = text_size,
        pt.cex = text_size * 2,
        bty = "n",
        ncol = number_columns,
        xpd = TRUE
    )
}
