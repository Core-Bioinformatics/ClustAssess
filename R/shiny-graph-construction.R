####### UI #######

#' UI - Graph construction module
#'
#' @description Creates the UI interface for the graph construction module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to identify the UI elements.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
ui_graph_construction <- function(id) {
    ns <- shiny::NS(id)
    shiny::absolutePanel(
        id = "selection_info",
        height = "20%",
        width = "100%",
        top = "100%"
    )
    shiny::tabPanel(
        "Graph Construction",
        shinyWidgets::circleButton(ns("info_title"),
            icon = shiny::icon("info"),
            size = "sm",
            status = "info",
            class = "page-info"
        ),
        shiny::h2("Relationship between the number of neighbours and the number of connected components", class = "first-element-tab"),
        shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info"), status = "success", icon = shiny::icon("info"), size = "sm"),
            shiny::h3(shiny::strong("Relationship between the number of neighbours and the number of connected components")),
            shiny::br(),
            shiny::h5("This plot describes the covariation between the number neighbours and the number of connected components obtained using both PCA and UMAP reductions as base for graph building. As the number of neighbours increases, the number of connected components decreases (this is an expected result, as increasing the number of neighbours result in a better connected graph). Please note that the number of connected components provides a lower bound on the number of clusters we can obtain by downstream community detection algorithms such as Louvain and Leiden."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href = "https://github.com/Core-Bioinformatics/ClustAssess", target = "_blank")),
            placement = "right",
            arrow = F,
            maxWidth = "700px"
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("download"),
            status = "success",
            size = "sm",
            shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
            shiny::textInput(ns("filename_con_comps"), "File name:", width = "80%"),
            shiny::numericInput(ns("width_con_comps"), "Width (in):", 7, 3, 100, 0.1),
            shiny::numericInput(ns("height_con_comps"), "Height (in):", 7, 3, 100, 0.1),
            shiny::selectInput(ns("filetype_con_comps"), "Filetype", choices = c("PDF", "PNG", "SVG"), width = "50%", selected = "PDF"),
            shiny::downloadButton(ns("download_con_comps"), label = "Download Plot")
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("cog"),
            status = "success",
            size = "sm",
            shiny::selectInput(ns("palette_plot_conn_comps"), "Colour scheme:",
                c(
                    "Colour scheme 1" = "col_1",
                    "Colour scheme 2" = "col_2",
                    "Colour scheme 3" = "col_3",
                    "Colour scheme 4" = "col_4"
                ),
                size = 4, selectize = F, selected = "col_4"
            ),
            shiny::sliderInput(
                inputId = ns("title_size1"),
                label = "Title size",
                min = 0,
                max = 40,
                value = 15,
                step = 0.5
            ),
            shiny::sliderInput(
                inputId = ns("axis_size1"),
                label = "Axis text size",
                min = 0,
                max = 40,
                value = 10,
                step = 0.5
            ),
            shiny::sliderInput(
                inputId = ns("legend_size1"),
                label = "Legend text size",
                min = 0,
                max = 40,
                value = 10,
                step = 0.5
            )
        ),
        shiny::splitLayout(
            cellWidths = c("15%", "85%"), shiny::div(style = "width:90%;", shiny::verticalLayout()),
            shiny::plotOutput(ns("plot_conn_comps"))
        ),
        shiny::h2("Relationship between the number of neighbours and then number of clusters"),
        shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info_2"), status = "success", icon = shiny::icon("info"), size = "sm"),
            shiny::h3(shiny::strong("Relationship between the number of neighbours and the number of connected components")),
            shiny::br(),
            shiny::h5("This plot summarizes the relationship between the number of nearest neighbours and k, the number of clusters. The plot also illustrates the difference between the two graph types: SNN has a tighter distribution of k (over multiple iterations) compared to NN. If the initial object contains graphs based on both UMAP and PCA embedding, this plot will showcase the impact of this choice, as well."),
            shiny::h5("This plot is interactive:"),
            shiny::h5("- You can select a desired set of configurations, and chose between three different colour schemes."),
            shiny::h5("- You can hover over the plot to show detailed information on the different boxes."),
            shiny::h5("- You can click anywhere in the plot to obtain a UMAP representation coloured with the stability (ECC) across different seeds for a specific number of nearest neighbours."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href = "https://github.com/Core-Bioinformatics/ClustAssess", target = "_blank")),
            placement = "right",
            arrow = F,
            maxWidth = "700px"
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("download"),
            status = "success",
            size = "sm",
            shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
            shiny::textInput(ns("filename_neigh_k"), "File name:", width = "80%"),
            shiny::numericInput(ns("width_neigh_k"), "Width (in):", 7, 3, 100, 0.1),
            shiny::numericInput(ns("height_neigh_k"), "Height (in):", 7, 3, 100, 0.1),
            shiny::selectInput(ns("filetype_neigh_k"), "Filetype", choices = c("PDF", "PNG", "SVG"), width = "50%", selected = "PDF"),
            shiny::downloadButton(ns("download_neigh_k"), label = "Download Plot")
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("cog"),
            status = "success",
            size = "sm",
            shiny::selectInput(ns("palette_plot_neigh_k"), "Colour scheme:",
                c(
                    "Colour scheme 1" = "col_1",
                    "Colour scheme 2" = "col_2",
                    "Colour scheme 3" = "col_3",
                    "Colour scheme 4" = "col_4"
                ),
                size = 4, selectize = F, selected = "col_4"
            ),
            shiny::uiOutput(ns("sel_conn_comps_render")),
            shiny::sliderInput(
                inputId = ns("title_size2"),
                label = "Title size",
                min = 0,
                max = 40,
                value = 15,
                step = 0.5
            ),
            shiny::sliderInput(
                inputId = ns("axis_size2"),
                label = "Axis text size",
                min = 0,
                max = 40,
                value = 10,
                step = 0.5
            ),
            shiny::sliderInput(
                inputId = ns("legend_size2"),
                label = "Legend text size",
                min = 0,
                max = 40,
                value = 10,
                step = 0.5
            )
        ),
        shiny::fluidRow(
            shiny::tags$head(shiny::tags$style(shiny::HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
            shiny::splitLayout(
                cellWidths = c("15%", "55%", "30%"),
                shiny::div(style = "width:90%;", shiny::verticalLayout(
                    shiny::htmlOutput(ns("hover_info")),
                    shiny::h1("\n")
                )),
                shiny::plotOutput(ns("plot_neigh_k"), hover = shiny::hoverOpts(id = ns("nn_hover"), delayType = "throttle"), click = shiny::clickOpts(id = ns("ecc_click")), width = "100%", height = "700px"),
                shiny::div(style = "width:90%;", shiny::verticalLayout(shiny::verbatimTextOutput(ns("click_info")), shiny::plotOutput(ns("umap_ecc"))))
            )
        ),
        shiny::h2("Stability"),
        shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info_3"), status = "success", icon = shiny::icon("info"), size = "sm"),
            shiny::h3(shiny::strong("Stability of different configurations across a varying number of nearest neighbours")),
            shiny::br(),
            shiny::h5("The stability of these parameters can be also evaluated using the Element-Centric Consistency applied on the partition list obtained over the multiple runs. The following summary plot underlines the evolution of the consistency as the number of neighbours increases."),
            shiny::h5("This plot is interactive:"),
            shiny::h5("- You can select a desired set of configurations, and chose between three different colour schemes."),
            shiny::h5("- You can hover over the plot to show detailed information on the different boxes."),
            shiny::h5("- You can click anywhere in the plot to obtain a UMAP representation coloured with the stability (ECC) across different seeds for a specific number of nearest neighbours."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href = "https://github.com/Core-Bioinformatics/ClustAssess", target = "_blank")),
            placement = "right",
            arrow = F,
            maxWidth = "700px"
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("download"),
            status = "success",
            size = "sm",
            shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
            shiny::textInput(ns("filename_neigh_stability"), "File name:", width = "80%", value = "neigh_stability"),
            shiny::numericInput(ns("width_neigh_stability"), "Width (in):", 7, 3, 100, 0.1),
            shiny::numericInput(ns("height_neigh_stability"), "Height (in):", 7, 3, 100, 0.1),
            shiny::selectInput(ns("filetype_neigh_stability"), "Filetype", choices = c("PDF", "PNG", "SVG"), width = "50%", selected = "PDF"),
            shiny::downloadButton(ns("download_neigh_stability"), label = "Download Plot")
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("cog"),
            status = "success",
            size = "sm",
            shiny::selectInput(ns("palette_neigh_stability"), "Colour scheme:",
                c(
                    "Colour scheme 1" = "col_1",
                    "Colour scheme 2" = "col_2",
                    "Colour scheme 3" = "col_3",
                    "Colour scheme 4" = "col_4"
                ),
                size = 4, selectize = F, selected = "col_4"
            ),
            shiny::uiOutput(ns("sel_stab_render")),
            shiny::sliderInput(
                inputId = ns("title_size3"),
                label = "Title size",
                min = 0,
                max = 40,
                value = 15,
                step = 0.5
            ),
            shiny::sliderInput(
                inputId = ns("axis_size3"),
                label = "Axis text size",
                min = 0,
                max = 40,
                value = 10,
                step = 0.5
            ),
            shiny::sliderInput(
                inputId = ns("legend_size3"),
                label = "Legend text size",
                min = 0,
                max = 40,
                value = 10,
                step = 0.5
            )
        ),
        shiny::fluidRow(
            shiny::splitLayout(
                cellWidths = c("15%", "55%", "30%"),
                shiny::verticalLayout(
                    shiny::div(style = "width:90%;"),
                    shiny::htmlOutput(ns("hover_info_nnstab")),
                    shiny::h1("\n")
                ),
                shiny::plotOutput(ns("plot_neigh_stability"), hover = shiny::hoverOpts(id = ns("nn_stability_hover"), delayType = "throttle"), click = shiny::clickOpts(id = ns("neigh_stability_click")), width = "100%", height = "700px"),
                shiny::div(style = "width:90%;", shiny::verticalLayout(shiny::verbatimTextOutput(ns("click_info_2")), shiny::plotOutput(ns("neigh_umap_ecc"))))
            )
        ),
        style = "margin-left:25px;margin-top:72px"
    )
}
####### SERVER #######

#' Server - Graph construction module
#'
#' @description Creates the backend interface for the graph construction
#' module inside the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to acess the UI elements.
#' @param chosen_config A reactive object that contains the chosen configuration
#' from the Dimensionality Reduction tab.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
server_graph_construction <- function(id, chosen_config) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            # Reactive for first plot

            chosen_config <- shiny::isolate(chosen_config())
            temp_list <- rhdf5::h5read("stability.h5", paste(unlist(chosen_config[1]), unlist(chosen_config[2]), sep = "/"))
            temp_list$pca <- NULL
            temp_list$stable_config <- NULL
            temp_list$clustering_stability <- NULL
            add_env_variable("stab_obj", temp_list)

            add_env_variable("current_tab", "Graph Construction")
            rm(temp_list)
            gc()

            options <- names(pkg_env$stab_obj$nn_stability$n_different_partitions)
            output$sel_conn_comps_render <- shiny::renderUI({
                ns <- session$ns
                shiny::checkboxGroupInput(
                    inputId = ns("sel_conn_comps"),
                    label = "Select configurations:",
                    choices = options,
                    selected = options,
                    width = "50%"
                )
            })
            shiny::outputOptions(output, "sel_conn_comps_render", suspendWhenHidden = FALSE)

            options <- names(pkg_env$stab_obj$nn_stability$n_different_partitions)
            output$sel_stab_render <- shiny::renderUI({
                ns <- session$ns
                shiny::checkboxGroupInput(
                    inputId = ns("sel_stab"),
                    label = "Select configurations:",
                    choices = options,
                    selected = options,
                    width = "50%"
                )
            })
            shiny::outputOptions(output, "sel_stab_render", suspendWhenHidden = FALSE)

            plot_conn_comps_Input <- shiny::reactive({
                title_size <- input$title_size1
                axis_size <- input$axis_size1
                legend_size <- input$legend_size1

                temp_gplot_obj <- plot_connected_comps_evolution(pkg_env$stab_obj$nn_conn_comps) +
                    ggplot2::theme(
                        axis.title = ggplot2::element_text(size = axis_size),
                        axis.text = ggplot2::element_text(size = axis_size),
                        legend.text = ggplot2::element_text(size = legend_size),
                        legend.title = ggplot2::element_text(size = legend_size),
                        plot.title = ggplot2::element_text(size = title_size)
                    )

                if (input$palette_plot_conn_comps == "col_4") {

                    return(
                        temp_gplot_obj +
                            ggplot2::scale_fill_manual(values = pkg_env$discrete_colors[[as.character(length(pkg_env$stab_obj$nn_conn_comps))]])
                    )
                }

                if (input$palette_plot_conn_comps == "col_1") {
                    option <- "RdBu"
                } else if (input$palette_plot_conn_comps == "col_2") {
                    option <- "RdYIGn"
                } else {
                    option <- "Greys"
                }
                temp_gplot_obj + ggplot2::scale_fill_brewer(palette = option)
            })
            # 1 Plot Relationship between the number of neighbours and then number of connected components
            output$plot_conn_comps <- shiny::renderPlot(
                {
                    plot_conn_comps_Input()
                },
                height = 400,
                width = 1000
            )
            # generate filenames
            con_comps_filetype <- shiny::reactive({
                if (input$filetype_con_comps == "PDF") {
                    filename <- paste0(input$filename_con_comps, ".pdf")
                    return(filename)
                } else if (input$filetype_con_comps == "PNG") {
                    filename <- paste0(input$filename_con_comps, ".png")
                    return(filename)
                } else {
                    filename <- paste0(input$filename_con_comps, ".svg")
                    return(filename)
                }
            })
            # Download Plots
            output$download_con_comps <- shiny::downloadHandler(
                filename = function() {
                    con_comps_filetype()
                },
                content = function(file) {
                    ggplot2::ggsave(file, plot_conn_comps_Input(),
                        width = input$width_con_comps,
                        height = input$height_con_comps,
                        units = "in",
                        limitsize = FALSE
                    )
                }
            )
            # Add hovering output for second plot, and add the text box
            output$hover_info <- shiny::renderUI({
                if (is.null(input$nn_hover)) {
                    shiny::HTML(paste(shiny::strong("Hover over the plot to see more")))
                } else {
                    col <- round(input$nn_hover$x)
                    obj <- pkg_env$stab_obj$nn_stability$n_neigh_k_corresp[input$sel_conn_comps]
                    obj <- reshape2::melt(obj)
                    val <- sort.int(as.integer(unique(obj$L2)))[col]
                    configs <- unique(obj$L1)

                    out <- paste0(shiny::strong("Selected number of neighbours: "), val, "<br/>", "<br/>")
                    for (config in configs) {
                        filt_obj <- obj[obj$L1 == config, ]

                        maxi <- max(filt_obj$value[filt_obj$L2 == val])
                        mini <- min(filt_obj$value[filt_obj$L2 == val])
                        iqr <- stats::IQR(filt_obj$value[filt_obj$L2 == val])

                        tmp_out <- paste0(shiny::strong(config), "<br/>", "Max: ", maxi, "<br/>", "Min: ", mini, "<br/>", "IQR: ", iqr, "<br/>", "<br/>")
                        out <- paste0(out, tmp_out)
                    }
                    shiny::HTML(paste(out, sep = ""))
                }
            })
            # Reactive for second plot
            plot_n_neigh_k_correspondence_Input <- shiny::reactive({
                title_size <- input$title_size2
                axis_size <- input$axis_size2
                legend_size <- input$legend_size2

                if (length(input$sel_conn_comps) == 0) {
                    return(ggplot2::ggplot() +
                        ggplot2::theme_void())
                }

                obj <- pkg_env$stab_obj
                obj$nn_stability$n_neigh_k_corresp <- obj$nn_stability$n_neigh_k_corresp[input$sel_conn_comps]

                temp_gplot_obj <- plot_n_neigh_k_correspondence(pkg_env$stab_obj$nn_stability) +
                    ggplot2::theme(
                        axis.title = ggplot2::element_text(size = axis_size),
                        axis.text = ggplot2::element_text(size = axis_size),
                        legend.text = ggplot2::element_text(size = legend_size),
                        legend.title = ggplot2::element_text(size = legend_size),
                        plot.title = ggplot2::element_text(size = title_size)
                    )


                if (input$palette_plot_neigh_k == "col_4") {
                    dis_colors <- pkg_env$discrete_colors[[as.character(length(input$sel_conn_comps))]]
                    names(dis_colors) <- names(obj$nn_stability$n_neigh_k_corresp)
                    return(
                        temp_gplot_obj + ggplot2::scale_fill_manual(values = dis_colors)
                    )
                }

                if (input$palette_plot_neigh_k == "col_1") {
                    option <- "RdBu"
                } else if (input$palette_plot_neigh_k == "col_2") {
                    option <- "RdYIGn"
                } else {
                    option <- "Greys"
                }
                
                temp_gplot_obj + ggplot2::scale_fill_brewer(palette = option)
            })
            # 2 Plot Relationship between the number of neighbours and then number of clusters
            output$plot_neigh_k <- shiny::renderPlot({
                plot_n_neigh_k_correspondence_Input()
            })
            neigh_k_filetype <- shiny::reactive({
                if (input$filetype_neigh_k == "PDF") {
                    filename <- paste0(input$filename_neigh_k, ".pdf")
                    return(filename)
                } else if (input$filetype_neigh_k == "PNG") {
                    filename <- paste0(input$filename_neigh_k, ".png")
                    return(filename)
                } else {
                    filename <- paste0(input$filename_neigh_k, ".svg")
                    return(filename)
                }
            })
            # Download Plots
            output$download_neigh_k <- shiny::downloadHandler(
                filename = function() {
                    neigh_k_filetype()
                },
                content = function(file) {
                    ggplot2::ggsave(file, plot_n_neigh_k_correspondence_Input(),
                        width = input$width_neigh_k,
                        height = input$height_neigh_k,
                        units = "in"
                    )
                }
            )
            # Plot the UMAP + add some text on the chosen selection
            output$click_info <- shiny::renderPrint({
                if (is.null(input$ecc_click)) {
                    return(cat("Click the plot to see more"))
                }
                col <- round(input$ecc_click$x)
                obj <- pkg_env$stab_obj$nn_stability$n_neigh_k_corresp[input$sel_conn_comps]
                obj <- reshape2::melt(obj)
                val <- sort.int(as.integer(unique(obj$L2)))[col]
                configs <- unique(obj$L1)
                if (input$ecc_click$x < 0.5) {
                    selection <- sort(names(pkg_env$stab_obj$nn_stability$n_different_partitions))[1]
                } else {
                    steps <- length(configs)
                    input_steps <- (input$ecc_click$x - (0.5 + (col - 1))) * steps
                    selected_config <- ceiling(input_steps)
                    selection <- sort(configs)[selected_config]
                }
                cat(paste0(selection, "\n", val, " NNs"))
            })
            output$umap_ecc <- shiny::renderPlot({
                if (is.null(input$ecc_click)) {
                    return(ggplot2::ggplot() +
                        ggplot2::theme_void())
                }
                # get click info for ecc
                col <- round(input$ecc_click$x)
                obj <- pkg_env$stab_obj$nn_stability$n_neigh_k_corresp[input$sel_conn_comps]
                obj <- reshape2::melt(obj)
                val <- sort.int(as.integer(unique(obj$L2)))[col]
                configs <- unique(obj$L1)
                if (input$ecc_click$x < 0.5) {
                    selection <- sort(names(pkg_env$stab_obj$nn_stability$n_different_partitions))[1]
                } else {
                    steps <- length(configs)
                    input_steps <- (input$ecc_click$x - (0.5 + (col - 1))) * steps
                    selected_config <- ceiling(input_steps)
                    selection <- sort(names(pkg_env$stab_obj$nn_stability$n_different_partitions))[selected_config]
                }
                ecc <- pkg_env$stab_obj$nn_stability$n_neigh_ec_consistency[input$sel_stab]
                ecc <- ecc[[selection]][[col]]

                ggplot2::ggplot(
                    data.frame(
                        pkg_env$stab_obj$umap
                    ),
                    ggplot2::aes(
                        x = .data$X1,
                        y = .data$X2,
                        color = ecc
                    )
                ) +
                    ggplot2::xlab("UMAP 1") +
                    ggplot2::ylab("UMAP 2") +
                    ggplot2::geom_point() +
                    ggplot2::scale_color_viridis_c() +
                    ggplot2::theme_bw() +
                    ggplot2::coord_fixed()
            })
            # Reactive for third plot
            plot_neigh_stability_Input <- shiny::reactive({
                if (length(input$sel_stab) == 0) {
                    return(ggplot2::ggplot() +
                        ggplot2::theme_void())
                }
                title_size <- input$title_size3
                axis_size <- input$axis_size3
                legend_size <- input$legend_size3

                obj <- pkg_env$stab_obj
                obj$nn_stability$n_neigh_ec_consistency <- obj$nn_stability$n_neigh_ec_consistency[input$sel_stab]
                temp_gplot_obj <- plot_n_neigh_ecs(obj$nn_stability, boxplot_width = 0.8) +
                    ggplot2::theme(
                        axis.title = ggplot2::element_text(size = axis_size),
                        axis.text = ggplot2::element_text(size = axis_size),
                        legend.text = ggplot2::element_text(size = legend_size),
                        legend.title = ggplot2::element_text(size = legend_size),
                        plot.title = ggplot2::element_text(size = title_size)
                    )

                if (input$palette_neigh_stability == "col_4") {
                    dis_colours <- pkg_env$discrete_colors[[as.character(length(input$sel_stab))]]
                    names(dis_colours) <- input$sel_stab
                    return(
                        temp_gplot_obj + ggplot2::scale_fill_manual(values = dis_colours)
                    )
                }

                if (input$palette_neigh_stability == "col_1") {
                    option <- "RdBu"
                } else if (input$palette_neigh_stability == "col_2") {
                    option <- "RdYIGn"
                } else {
                    option <- "Greys"
                }

                temp_gplot_obj + ggplot2::scale_fill_brewer(palette = option)
            })
            # 3. Plot Stability
            output$plot_neigh_stability <- shiny::renderPlot({
                plot_neigh_stability_Input()
            })
            # generate filenames
            neigh_stability_filetype <- shiny::reactive({
                if (input$filetype_neigh_stability == "PDF") {
                    filename <- paste0(input$filename_neigh_stability, ".pdf")
                    return(filename)
                } else if (input$filetype_neigh_stability == "PNG") {
                    filename <- filename <- paste0(input$filename_neigh_stability, ".png")
                    return(filename)
                } else {
                    filename <- filename <- paste0(input$filename_neigh_stability, ".svg")
                    return(filename)
                }
            })
            # Download Plots
            output$download_neigh_stability <- shiny::downloadHandler(
                filename = function() {
                    neigh_stability_filetype()
                },
                content = function(file) {
                    ggplot2::ggsave(file, plot_neigh_stability_Input(),
                        width = input$width_neigh_stability,
                        height = input$height_neigh_stability,
                        units = "in"
                    )
                }
            )
            # Hover info
            hoverdf <- pkg_env$stab_obj$nn_stability$n_neigh_ec_consistency
            hoverdf <- reshape2::melt(hoverdf)
            hoverdf$value <- round(hoverdf$value, digits = 2)

            output$hover_info_nnstab <- shiny::renderUI({
                if (is.null(input$nn_stability_hover)) {
                    shiny::HTML(paste(shiny::strong("Hover over the plot to see more")))
                } else {
                    col <- round(input$nn_stability_hover$x)
                    obj <- subset(hoverdf, .data$L1 == input$sel_stab)
                    val <- sort.int(as.integer(unique(obj$L2)))[col]
                    configs <- unique(obj$L1)

                    out <- paste0(shiny::strong("Selected number of neighbours: "), val, "<br/>", "<br/>")
                    for (config in configs) {
                        filt_obj <- obj[obj$L1 == config, ]

                        maxi <- max(filt_obj$value[filt_obj$L2 == val])
                        mini <- min(filt_obj$value[filt_obj$L2 == val])
                        iqr <- stats::IQR(filt_obj$value[filt_obj$L2 == val])

                        tmp_out <- paste0(shiny::strong(config), "<br/>", "Max: ", maxi, "<br/>", "Min: ", mini, "<br/>", "IQR: ", iqr, "<br/>", "<br/>")
                        out <- paste0(out, tmp_out)
                    }
                    shiny::HTML(paste(out, sep = ""))
                }
            })
            # Plot UMAP + add some text on the chosen selection
            output$click_info_2 <- shiny::renderPrint({
                if (is.null(input$neigh_stability_click)) {
                    return(cat("Click the plot to see more"))
                }
                col <- round(input$neigh_stability_click$x)
                obj <- pkg_env$stab_obj$nn_stability$n_neigh_ec_consistency[input$sel_stab]
                obj <- reshape2::melt(obj)
                val <- sort.int(as.integer(unique(obj$L2)))[col]
                configs <- unique(obj$L1)
                if (input$neigh_stability_click$x < 0.5) {
                    selection <- sort(names(pkg_env$stab_obj$nn_stability$n_different_partitions))[1]
                } else {
                    steps <- length(configs)
                    input_steps <- (input$neigh_stability_click$x - (0.5 + (col - 1))) * steps
                    selected_config <- ceiling(input_steps)
                    selection <- sort(configs)[selected_config]
                }
                cat(paste0(selection, "\n", val, " NNs"))
            })
            output$neigh_umap_ecc <- shiny::renderPlot({
                if (is.null(input$neigh_stability_click)) {
                    return(ggplot2::ggplot() +
                        ggplot2::theme_void())
                }
                # get click info for ecc
                col <- round(input$neigh_stability_click$x)
                obj <- pkg_env$stab_obj$nn_stability$n_neigh_ec_consistency[input$sel_stab]
                obj <- reshape2::melt(obj)
                val <- sort.int(as.integer(unique(obj$L2)))[col]
                configs <- unique(obj$L1)
                if (input$neigh_stability_click$x < 0.5) {
                    selection <- sort(names(pkg_env$stab_obj$nn_stability$n_different_partitions))[1]
                } else {
                    steps <- length(configs)
                    input_steps <- (input$neigh_stability_click$x - (0.5 + (col - 1))) * steps
                    selected_config <- ceiling(input_steps)
                    selection <- sort(configs)[selected_config]
                }
                ecc <- pkg_env$stab_obj$nn_stability$n_neigh_ec_consistency[input$sel_stab]
                ecc <- ecc[[selection]][[col]]

                ggplot2::ggplot(
                    data.frame(
                        pkg_env$stab_obj$umap
                    ),
                    ggplot2::aes(
                        x = .data$X1,
                        y = .data$X2,
                        color = ecc
                    )
                ) +
                    ggplot2::xlab("UMAP 1") +
                    ggplot2::ylab("UMAP 2") +
                    ggplot2::geom_point() +
                    ggplot2::scale_color_viridis_c() +
                    ggplot2::theme_bw() +
                    ggplot2::coord_fixed()
            })
        }
    )
}
