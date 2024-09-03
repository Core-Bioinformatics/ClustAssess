####### UI #######

ui_comparison_markers <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Identification of markers"),
        shiny::splitLayout(
            cellWidths = c("90%", "10%"),
            shiny::plotOutput(ns("avg_expression_violin"), height = "auto"),
            shiny::tableOutput(ns("avg_expression_table"))
        ),
        shinyWidgets::dropdownButton(
            shiny::tagList(
                shiny::sliderInput(
                    inputId = ns("logfc"),
                    label = "logFC threshold",
                    min = 0.00, max = 10.00, value = 0.50, step = 0.01
                ),
                shiny::sliderInput(
                    inputId = ns("avg_expr_thresh"),
                    label = "Average expression threshold",
                    min = 0.00, max = 0.01, value = 0.00
                ),
                shiny::sliderInput(
                    inputId = ns("avg_expr_thresh_gr1"),
                    label = "Average expression threshold - group 1",
                    min = 0.00, max = 0.01, value = 0.00
                ),
                shiny::sliderInput(
                    inputId = ns("min_pct"),
                    label = "Minimum gene frequency",
                    min = 0.01, max = 1.00, value = 0.10, step = 0.01
                ),
                shiny::sliderInput(
                    inputId = ns("pval"),
                    label = "Maximum adj-pval",
                    min = 0.001, max = 1.00, value = 0.01, step = 0.001
                ),
                shinyWidgets::prettySwitch(
                    inputId = ns("norm_type"),
                    label = "Data is normalised",
                    value = TRUE,
                    status = "success",
                    fill = TRUE
                )
            ),
            circle = TRUE,
            status = "success",
            size = "sm",
            icon = shiny::icon("cog")
        ),
        shiny::htmlOutput(ns("marker_text")),
        shiny::actionButton(ns("enable_markers"),
            "Enable DEG analysis",
            style = "font-size:20px;",
            class = "btn-danger"
        ),
        shiny::fluidRow(
            shiny::column(
                6,
                ui_comparison_markers_panel(ns("group_left"))
            ),
            shiny::column(
                6,
                ui_comparison_markers_panel(ns("group_right"))
            )
        ),
        shiny::actionButton(ns("markers_button"), "Find markers!", class = "btn-danger"),
        DT::dataTableOutput(ns("markers_dt")),
        shiny::downloadButton(ns("markers_download_button"), "Download markers!")
    )
}

ui_comparison_markers_panel <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::selectInput(
            inputId = ns("select_k_markers"),
            label = "Select the number of clusters (k) or metadata",
            choices = ""
        ),
        shinyWidgets::pickerInput(
            inputId = ns("select_clusters_markers"),
            label = "Select the groups of cells",
            choices = "",
            inline = FALSE,
            options = list(
                `actions-box` = TRUE,
                title = "Select/deselect subgroups",
                size = 10,
                width = "90%",
                `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
        )
    )
}

ui_comparison_metadata_panel <- function(id, draw_line) {
    ns <- shiny::NS(id)
    style <- ifelse(draw_line, "border-right:5px solid", "")

    shinyWidgets::panel(
        style = style,
        shiny::selectizeInput(
            inputId = ns("metadata"),
            label = "Metadata",
            choices = NULL
        ),
        shiny::splitLayout(
            shiny::selectInput(
                inputId = ns("metadata_subset"),
                choices = NULL,
                label = "Subset by metadata"
            ),
            shiny::verticalLayout(
                shiny::tags$b("Select groups"),
                shinyWidgets::pickerInput(
                    inputId = ns("metadata_groups_subset"),
                    choices = NULL,
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
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            gear_umaps(ns, "metadata"),
            gear_download(ns, "metadata", "metadata")
        ),
        shiny::plotOutput(ns("umap_metadata"), height = "auto"),
        shiny::plotOutput(ns("umap_metadata_legend"), height = "auto")
    )
}

ui_comparison_gene_panel <- function(id, draw_line) {
    ns <- shiny::NS(id)
    style <- ifelse(draw_line, "border-right:5px solid", "")

    shinyWidgets::panel(
        style = style,
        shiny::selectizeInput(
            inputId = ns("gene_expr"),
            choices = NULL,
            label = "Gene name(s)",
            width = "95%",
            multiple = TRUE
        ),
        shiny::verticalLayout(
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
                shiny::selectInput(
                    inputId = ns("metadata_subset"),
                    choices = NULL,
                    label = "Subset by metadata"
                ),
                shiny::verticalLayout(
                    shiny::tags$b("Select groups"),
                    shinyWidgets::pickerInput(
                        inputId = ns("metadata_groups_subset"),
                        choices = NULL,
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
            )
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            gear_umaps(ns, "gene", FALSE, "highest"),
            gear_download(ns, "gene", "gene"),
        ),
        shiny::plotOutput(ns("umap_gene"), height = "auto"),
        shiny::plotOutput(ns("umap_gene_legend"), height = "auto")
    )
}

ui_comparison_jsi_panel <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h1("Jaccard Simmilarity Index (JSI)/Cells per cluster"),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            shinyWidgets::dropdownButton(
                label = "",
                icon = shiny::icon("cog"),
                status = "success",
                size = "sm",
                shiny::radioButtons(ns("heatmap_type"), "Calculate similarity", choices = c("JSI", "Cells per cluster"), width = "100%")
            ),
            shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info"), status = "success", icon = shiny::icon("info"), size = "sm"),
                shiny::h3(shiny::strong("Jaccard Simmilarity Index (JSI) between clusters")),
                shiny::br(),
                shiny::h5("This plot aims to showcase the behaviour of the individual clusters on the different partitions. JSI is calculated for the cell barcodes for every cluster, in both configurations, in a pair-wise manner."),
                shiny::h1("\n"),
                shiny::h5("For more information please go to:"),
                shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href = "https://github.com/Core-Bioinformatics/ClustAssess", target = "_blank")),
                placement = "right",
                arrow = FALSE
            ),
            #  maxWidth = '700px'),
            shinyWidgets::dropdownButton(
                label = "",
                icon = shiny::icon("download"),
                status = "success",
                size = "sm",
                shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                shiny::textInput(ns("filename_heatmap"), "File name:", width = "80%"),
                shiny::numericInput(ns("width_heatmap"), "Width (in):", 7, 3, 100, 0.1),
                shiny::numericInput(ns("height_heatmap"), "Height (in):", 7, 3, 100, 0.1),
                shiny::selectInput(ns("heatmap_filetype"), "Filetype", choices = c("PDF", "PNG", "SVG"), selected = "PDF", width = "100%"),
                shiny::downloadButton(ns("download_heatmap"), label = "Download Plot")
            )
        ),
        shiny::selectInput(
            inputId = ns("jsi_k_1"),
            label = "Select the number of clusters (k) or metadata, for the first comparison",
            choices = ""
        ),
        shiny::selectInput(
            inputId = ns("jsi_k_2"),
            label = "Select the number of clusters (k)or metadata, for the second comparison",
            choices = ""
        ),
        shiny::plotOutput(ns("barcode_heatmap"), height = "auto", width = "98%")
    )
}

ui_comparison_gene_heatmap <- function(id) {
    ns <- shiny::NS(id) 

    shiny::tagList(
        shiny::h2("Gene expression heatmap"),
        shiny::splitLayout(
            cellWidths = "40px",
            shinyWidgets::dropdownButton(
                shiny::sliderInput(
                    inputId = ns("text_size"),
                    label = "Text size",
                    min = 5, max = 50, value = 15, step = 0.5
                ),
                shiny::sliderInput(
                    inputId = ns("point_size"),
                    label = "Point max size",
                    min = 0.1, max = 30, value = c(0.5, 3), step = 0.1
                ),
                shinyWidgets::prettySwitch(
                    inputId = ns("scale"),
                    label = "Apply scaling",
                    status = "success",
                    fill = TRUE,
                    value = FALSE
                ),
                shiny::sliderInput(
                    inputId = ns("clipping_value"),
                    label = "Clipping value",
                    min = 0.01, max = 20, value = 6, step = 0.1
                ),
                shinyWidgets::prettySwitch(
                    inputId = ns("show_numbers"),
                    label = "Show the values on the heatmap",
                    status = "success",
                    fill = TRUE,
                    value = FALSE 
                ),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "heatmap", "heatmap")
        ),
        shiny::splitLayout(
            shiny::selectizeInput(
                inputId = ns("gene_expr"),
                choices = NULL,
                label = "Select a metadata or a gene",
                width = "95%",
                multiple = TRUE
            ),
            shiny::selectInput(
                inputId = ns("metadata"),
                choices = NULL,
                label = "Split by metadata"
            ),
            shinyWidgets::radioGroupButtons(
                inputId = ns("plot_type"),
                label = "Plot type",
                choices = c("Heatmap", "Bubbleplot"),
                selected = "Heatmap"
            )
        ),
        shiny::plotOutput(ns("gene_heatmap"), height = "auto"),
    )
}

ui_comparison_violin_gene <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Violin Plots - Gene expression"),
        shiny::splitLayout(
            cellWidths = "40px",
            shinyWidgets::dropdownButton(
                shiny::sliderInput(
                    inputId = ns("text_size"),
                    label = "Text size",
                    min = 5, max = 30, value = 10, step = 0.5
                ),
                shinyWidgets::prettySwitch(
                    inputId = ns("log_scale"),
                    label = "Use log transform",
                    status = "success",
                    fill = TRUE,
                    value = FALSE
                ),
                shiny::sliderInput(
                    inputId = ns("boxplot_width"),
                    label = "Boxplot width",
                    min = 0.01, max = 1, value = 0.1, step = 0.01
                ),
                shiny::sliderInput(
                    inputId = ns("boxplot_dodge"),
                    label = "Distance between boxplots",
                    min = 0.01, max = 1, value = 0.1, step = 0.01
                ),
                shinyWidgets::checkboxGroupButtons(
                    inputId = ns("graph_type"),
                    label = "Graph type (multiple)",
                    choices = c("Violin", "Boxplot"),
                    selected = c("Violin", "Boxplot")
                ),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "violin", "violin")
        ),
        shiny::splitLayout(
            shiny::selectizeInput(
                inputId = ns("gene_expr"),
                choices = NULL,
                label = "Select a metadata or a gene",
                # width = "95%",
                multiple = FALSE
            ),
            shiny::selectInput(
                inputId = ns("metadata_split"),
                choices = NULL,
                label = "Split by metadata"
            ),
            shiny::selectInput(
                inputId = ns("metadata_group"),
                choices = NULL,
                label = "Group by metadata"
            ),
            shiny::selectInput(
                inputId = ns("metadata_subset"),
                label = "Subset by metadata",
                choices = NULL
            ),
            shiny::verticalLayout(
                shiny::tags$b("Select subset groups"),
                shinyWidgets::pickerInput(
                    inputId = ns("metadata_groups_subset"),
                    choices = NULL,
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
        shiny::plotOutput(ns("violin_gene"), height = "auto"),
        shiny::selectInput(
            inputId = ns("stat_mtd_group"),
            label = "Select the group for the stats table",
            choices = NULL
        ),
        shiny::tableOutput(ns("stats"))
    )
}

ui_comparison_enrichment <- function(id) {
    ns <- shiny::NS(id)

    shiny::div(
        shiny::h2("Enrichment analysis"),
        shiny::splitLayout(
            shinyWidgets::pickerInput(
                inputId = ns("gprofilerSources"),
                label = "Select data sources",
                choices = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP"),
                selected = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA"),
                multiple = TRUE,
                options = list(
                    title = "Select data sources"
                )
            ),
            shinyWidgets::radioGroupButtons(
                inputId = ns("group"),
                label = "Select the group containing the target markers",
                choices = c("group 1", "group 2")
            )
        ),
        shiny::actionButton(
            inputId = ns("enrichment_button"),
            label = "Perform enrichment analysis!",
            style = "font-size:20px;",
            class = "btn-danger"
        ),
        plotly::plotlyOutput(
            outputId = ns("gost_plot"),
            height = "auto"
        ),
        DT::DTOutput(outputId = ns("gost_table")),
        shiny::downloadButton(
            outputId = ns("download_gost"),
            label = "Download enriched terms as CSV",
            class = "btn-info"
        ),
        id = ns("enrichment_id")
    )
}

#' UI - Comparison module
#'
#' @description Creates the UI interface for the comparison module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to identify the UI elements.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
ui_comparisons <- function(id) {
    ns <- shiny::NS(id)
    shiny::tabPanel(
        "Comparison",
        shinyWidgets::circleButton(ns("info_title"),
            icon = shiny::icon("info"),
            size = "sm",
            status = "info",
            class = "page-info"
        ),
        shiny::actionButton(ns("show_config"), "Show config", type = "info", class = "btn-info show_config"),
        shiny::h2("Compare your current configuration", class = "first-element-tab"), # style = "margin-bottom:10px ")
        shiny::splitLayout(
            cellWidths = c("48%", "48%"),
            ui_comparison_metadata_panel(ns("metadata_panel_left"), TRUE),
            ui_comparison_metadata_panel(ns("metadata_panel_right"), FALSE)
        ),
        shiny::splitLayout(
            cellWidths = c("48%", "48%"),
            ui_comparison_gene_panel(ns("gene_panel_left"), TRUE),
            ui_comparison_gene_panel(ns("gene_panel_right"), FALSE)
        ),
        ui_comparison_jsi_panel(ns("jsi_plot")),
        ui_comparison_violin_gene(ns("violin_gene")),
        ui_comparison_gene_heatmap(ns("gene_heatmap")),
        ui_comparison_markers(ns("markers")),
        ui_comparison_enrichment(ns("enrichment")),
        style = "margin-bottom:30px;"
    )
}

####### SERVER #######
server_comparison_markers <- function(id, k_choices) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            if ("genes" %in% names(pkg_env)) { # for backward-compatibility purposes
                output$avg_expression_violin <- shiny::renderPlot(
                    {
                        vioplot::vioplot(
                            x = rhdf5::h5read("expression.h5", "average_expression"),
                            horizontal = TRUE,
                            xlab = "Average expression",
                            main = "Average gene expression",
                            ylab = "",
                            xaxt = "n"
                        )
                    },
                    height = function() {
                        400
                    }
                )

                avg_stats <- stats::fivenum(rhdf5::h5read("expression.h5", "average_expression"))

                output$avg_expression_table <- shiny::renderTable(
                    {
                        data.frame(
                            row.names = c("min", "Q1", "median", "Q3", "max"),
                            b = as.character(trunc(avg_stats * 1e5) / 1e5)
                        )
                    },
                    colnames = FALSE,
                    rownames = TRUE
                )

                shiny::updateSliderInput(
                    session = session,
                    inputId = "avg_expr_thresh",
                    min = round(avg_stats[1], digits = 3),
                    max = round(avg_stats[5], digits = 3),
                    step = 0.01
                )

                shiny::updateSliderInput(
                    session = session,
                    inputId = "avg_expr_thresh_gr1",
                    min = round(avg_stats[1], digits = 3),
                    max = round(avg_stats[5], digits = 3),
                    step = 0.01
                )
            }

            marker_genes <- shiny::reactiveVal(NULL)

            # it would be nice to have gene umaps
            shinyjs::html("marker_text", "Warning: Enabling DEG analysis will results into loading the memory. This process might take some time.")
            shinyjs::hide("group_left-select_k_markers")
            shinyjs::hide("group_left-select_clusters_markers")
            shinyjs::hide("group_right-select_k_markers")
            shinyjs::hide("group_right-select_clusters_markers")
            shinyjs::hide("markers_download_button")
            shinyjs::hide("markers_button")
            shinyjs::hide("markers_dt")
            shinyjs::show("enable_markers")

            server_comparison_markers_panels(session, k_choices)
            # server_comparison_markers_panels("group_left", k_choices)

            shiny::observe({
                current_button_value <- as.integer(shiny::isolate(input$enable_markers))
                shiny::req(pkg_env$enable_markers_button() != current_button_value)
                pkg_env$enable_markers_button(current_button_value)
                shinyjs::html("marker_text", "Preparing the objects for the analysis...")

                if (!("genes" %in% names(pkg_env))) {
                    expr_matrix <- rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(pkg_env$genes_of_interest[pkg_env$used_genes], NULL))
                    rownames(expr_matrix) <- pkg_env$used_genes

                    add_env_variable("rank_matrix", rhdf5::h5read("expression.h5", "rank_of_interest", index = list(pkg_env$genes_of_interest[pkg_env$used_genes], NULL)))
                    add_env_variable("expr_matrix", expr_matrix)
                }
                shinyjs::hide("enable_markers")
                shinyjs::show("group_left-select_k_markers")
                shinyjs::show("group_left-select_clusters_markers")
                shinyjs::show("group_right-select_k_markers")
                shinyjs::show("group_right-select_clusters_markers")
                shinyjs::show("markers_button")
                shinyjs::html("marker_text", "")
            }) %>% shiny::bindEvent(input$enable_markers)

            markers_val <- shiny::reactive({
                current_button_value <- as.integer(shiny::isolate(input$markers_button))
                shiny::req(
                    input$"group_left-select_clusters_markers",
                    input$"group_right-select_clusters_markers",
                    pkg_env$find_markers_button() != current_button_value
                )
                pkg_env$find_markers_button(current_button_value)

                shinyjs::disable("markers_button")
                shinyjs::html("marker_text", "Calculating the markers...")
                subgroup_left <- input$"group_left-select_k_markers"
                subgroup_right <- input$"group_right-select_k_markers"

                if (is.na(as.numeric(subgroup_left))) {
                    mb1 <- pkg_env$metadata[[subgroup_left]]
                } else {
                    mb1 <- factor(pkg_env$stab_obj$mbs[[subgroup_left]])
                }

                if (is.na(as.numeric(subgroup_right))) {
                    mb2 <- pkg_env$metadata[[subgroup_right]]
                } else {
                    mb2 <- factor(pkg_env$stab_obj$mbs[[subgroup_right]])
                }

                cells_index_left <- which(mb1 %in% input$"group_left-select_clusters_markers")
                cells_index_right <- which(mb2 %in% input$"group_right-select_clusters_markers")

                if ("genes" %in% names(pkg_env)) {
                    markers_result <- calculate_markers_shiny(
                        cells1 = cells_index_left,
                        cells2 = cells_index_right,
                        norm_method = ifelse(input$norm_type, "LogNormalize", ""),
                        used_slot = "data",
                        min_pct_threshold = input$min_pct,
                        logfc_threshold = input$logfc,
                        average_expression_threshold = input$avg_expr_thresh,
                        average_expression_group1_threshold = input$avg_expr_thresh_gr1
                    )
                } else { # for backward-compatibility reasons
                    markers_result <- calculate_markers(
                        expression_matrix = pkg_env$expr_matrix, # expression matrix
                        cells1 = cells_index_left,
                        cells2 = cells_index_right,
                        rank_matrix = pkg_env$rank_matrix, # rank matrix
                        norm_method = ifelse(input$norm_type, "LogNormalize", ""),
                        min_pct_threshold = input$min_pct,
                        logfc_threshold = input$logfc
                    )
                }

                all_genes <- as.vector(markers_result$gene)

                markers_result <- markers_result %>%
                    dplyr::filter(.data$p_val_adj <= input$pval) %>%
                    dplyr::arrange(dplyr::desc(.data$avg_log2FC), .data$p_val_adj)
                genes_group1 <- (markers_result %>% dplyr::filter(.data$avg_log2FC >= 0))$gene

                marker_genes(list(
                    all_genes = all_genes,
                    group_1 = as.vector(genes_group1),
                    group_2 = as.vector(markers_result$gene[seq(from = length(genes_group1) + 1, to = nrow(markers_result))])
                ))

                shinyjs::show("markers_dt")
                shinyjs::show("markers_download_button")
                shinyjs::enable("markers_button")
                shinyjs::html("marker_text", "")

                return(markers_result)
            }) %>% shiny::bindEvent(input$markers_button)


            shiny::observe(
                output$markers_dt <- DT::renderDataTable(
                    {
                        shiny::req(markers_val())
                        markers_val()
                    },
                    rownames = FALSE
                )
            ) %>% shiny::bindEvent(markers_val())

            output$markers_download_button <- shiny::downloadHandler(
                filename = function() {
                    "markers.csv"
                },
                content = function(file) {
                    utils::write.csv(markers_val(), file)
                }
            )

            shiny::observe({
                shiny::req(markers_val())
                shinyjs::show("markers_download_button")
            }) %>% shiny::bindEvent(markers_val())

            return(shiny::reactive(marker_genes()))
        }
    )
}

server_comparison_markers_panels <- function(session, k_choices) {
    available_choices <- c(names(pkg_env$metadata_unique), k_choices)
    input <- session$input

    shiny::updateSelectInput(
        session = session,
        inputId = "group_left-select_k_markers",
        choices = available_choices,
        selected = available_choices[1]
    )

    shiny::observe({
        shiny::req(input$"group_left-select_k_markers" %in% available_choices)

        if (is.na(as.numeric(input$"group_left-select_k_markers"))) {
            available_subgroups <- pkg_env$metadata_unique[[input$"group_left-select_k_markers"]]
        } else {
            available_subgroups <- seq_len(as.numeric(input$"group_left-select_k_markers"))
        }

        shinyWidgets::updatePickerInput(
            session = session,
            inputId = "group_left-select_clusters_markers",
            choices = available_subgroups,
            selected = available_subgroups[1]
        )

        shiny::updateSelectInput(
            session = session,
            inputId = "group_right-select_k_markers",
            choices = available_choices,
            selected = input$"group_left-select_k_markers"
        )
    }) %>% shiny::bindEvent(input$"group_left-select_k_markers")


    shiny::observe({
        shiny::req(input$"group_right-select_k_markers" %in% available_choices)

        if (is.na(as.numeric(input$"group_right-select_k_markers"))) {
            available_subgroups <- pkg_env$metadata_unique[[input$"group_right-select_k_markers"]]
        } else {
            available_subgroups <- seq_len(as.numeric(input$"group_right-select_k_markers"))
        }

        if (input$"group_left-select_k_markers" == input$"group_right-select_k_markers") {
            selected_groups <- available_subgroups[!(available_subgroups %in% input$"group_left-select_clusters_markers")]
        } else {
            selected_groups <- available_subgroups[1]
        }

        shinyWidgets::updatePickerInput(
            session = session,
            inputId = "group_right-select_clusters_markers",
            choices = available_subgroups,
            selected = selected_groups
        )
    })

    shiny::observe({
        shiny::updateSelectInput(
            session = session,
            inputId = "group_right-select_k_markers",
            choices = available_choices,
            selected = input$"group_left-select_k_markers"
        )
    }) %>% shiny::bindEvent(input$"group_left-select_clusters_markers")
}

server_comparison_metadata_panel <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            plt_height <- shiny::reactive(
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.43))
            )
            changed_metadata <- shiny::reactiveVal(FALSE)

            shiny::observe({
                shiny::req(input$metadata_subset)
                is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")

                if (is_cluster) {
                    k_value <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = "metadata_groups_subset",
                        choices = seq_len(k_value),
                        selected = seq_len(k_value),
                        clearOptions = TRUE
                    )
                    changed_metadata(TRUE)

                    return()
                }

                mtd_names <- pkg_env$metadata_unique[[input$metadata_subset]]

                shinyWidgets::updatePickerInput(
                    session,
                    inputId = "metadata_groups_subset",
                    choices = mtd_names,
                    selected = mtd_names
                )

                changed_metadata(TRUE)
            }) %>% shiny::bindEvent(input$metadata_subset)

            metadata_mask <- shiny::reactive({
                shiny::req(input$metadata_groups_subset, input$metadata_subset, cancelOutput = TRUE)
                shiny::isolate({
                    is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")
                    if (is_cluster) {
                        k <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                        all_unique_values <- as.character(seq_len(k))
                    } else {
                        all_unique_values <- pkg_env$metadata_unique[[input$metadata_subset]]
                    }

                    if (changed_metadata()) {
                        shiny::req(
                            all(input$metadata_groups_subset %in% all_unique_values),
                            cancelOutput = TRUE
                        )

                        changed_metadata(FALSE)
                    }

                    if (is_cluster) {
                        return(pkg_env$stab_obj$mbs[[as.character(k)]] %in% as.integer(input$metadata_groups_subset))
                    }
                    return(pkg_env$metadata[[input$metadata_subset]] %in% input$metadata_groups_subset)
                })
            })

            metadata_legend_height <- shiny::reactiveVal(0)

            plot_data <- shiny::reactive({
                shiny::req(input$metadata)

                is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
                is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

                if (is_cluster || is_ecc) {
                    k <- strsplit(input$metadata, "_")[[1]][2]
                    if (is_ecc) {
                        cl_method <- strsplit(names(pkg_env$stab_obj$ecc)[1], ";")[[1]][2]
                        unique_values <- NULL
                        color_values <- NULL
                        color_info <- pkg_env$stab_obj$ecc[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]
                        color_info <- color_info[pkg_env$stab_obj$ecc_order[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]]
                    } else {
                        color_values <- rhdf5::h5read("stability.h5", paste0("colors/", k))
                        unique_values <- seq_len(as.integer(k))
                        color_info <- factor(pkg_env$stab_obj$mbs[[k]])
                    }
                } else {
                    unique_values <- pkg_env$metadata_unique[[input$metadata]]
                    color_values <- pkg_env$metadata_colors[[input$metadata]]
                    color_info <- pkg_env$metadata[[input$metadata]]
                }

                list(
                    unique_values = unique_values,
                    color_values = color_values,
                    color_info = color_info
                )
            }) %>% shiny::bindEvent(input$metadata)

            output$umap_metadata <- shiny::renderPlot(
                height = function() {
                    plt_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    shiny::req(metadata_mask(), input$metadata_groups_subset, cancelOutput = TRUE)
                    plot_data()
                    input$metadata_subset
                    input$metadata_pt_size
                    input$metadata_axis_size
                    input$metadata_text_size
                    input$metadata_labels
                    input$metadata_pt_type
                    input$metadata_pt_order
                    input$metadata_legend_size
                    plt_height()

                    shiny::isolate({
                        if (is.null(plot_data()$unique_values)) {
                            old_par <- graphics::par(mai = c(0.1, 0, 0.1, 0))
                            text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$metadata_legend_size)
                        } else {
                            old_par <- graphics::par(mar = c(0, 0, 0, 0))
                            predicted_width <- graphics::strwidth(c(" ", plot_data()$unique_values), units = "inches", cex = input$metadata_legend_size) * ppi
                            space_width <- predicted_width[1]
                            predicted_width <- predicted_width[2:length(predicted_width)]

                            number_columns <- min(
                                max(
                                    plt_height() %/% (6 * space_width + max(predicted_width)),
                                    1
                                ),
                                length(plot_data()$unique_values)
                            )
                            number_rows <- ceiling(length(plot_data()$unique_values) / number_columns)

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
                        color_plot2(
                            embedding = pkg_env$stab_obj$umap,
                            color_info = plot_data()$color_info,
                            color_values = plot_data()$color_values,
                            unique_values = plot_data()$unique_values,
                            plt_height = plt_height(),
                            plt_width = plt_height(),
                            display_legend = FALSE,
                            pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
                            pt_size = input$metadata_pt_size,
                            sort_cells = input$metadata_pt_order,
                            text_size = input$metadata_text_size,
                            axis_size = input$metadata_axis_size,
                            labels = input$metadata_labels,
                            cell_mask = metadata_mask()
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
                        shiny::req(metadata_mask(), input$metadata_groups_subset, cancelOutput = TRUE)
                        plot_data()
                        plt_height()
                        input$select_groups
                        input$metadata_legend_size

                        shiny::isolate({
                            if (!is.null(plot_data()$unique_values)) {
                                unique_values <- unique(plot_data()$color_info[metadata_mask()])
                                unique_values <- plot_data()$unique_values[plot_data()$unique_values %in% unique_values]
                                color_values <- plot_data()$color_values[which(plot_data()$unique_values %in% unique_values)]
                            } else {
                                color_values <- NULL
                                unique_values <- NULL
                            }


                            only_legend_plot(
                                unique_values = unique_values,
                                color_values = color_values,
                                color_info = plot_data()$color_info,
                                plt_width = plt_height(),
                                text_size = input$metadata_legend_size
                            )
                        })
                    }
                )
            })

            output$download_metadata <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_metadata, ".", tolower(input$filetype_metadata))
                },
                content = function(file) {
                    shiny::req(input$metadata, input$width_metadata, input$height_metadata)

                    ggplot_obj <- color_ggplot(
                        embedding = pkg_env$stab_obj$umap,
                        color_info = plot_data()$color_info,
                        sort_cells = input$metadata_pt_order,
                        cell_mask = metadata_mask(),
                        legend_text_size = input$metadata_legend_size * 10,
                        axis_text_size = input$metadata_axis_size * 10,
                        text_size = input$metadata_text_size * 3,
                        labels = input$metadata_labels,
                        pt_size = input$metadata_pt_size
                    ) + ggplot2::ggtitle(input$metadata)

                    if (!is.null(plot_data()$unique_values)) {
                        color_vector <- plot_data()$color_values
                        names(color_vector) <- plot_data()$unique_values
                        ggplot_obj <- ggplot_obj +
                            ggplot2::scale_colour_manual(values = color_vector) +
                            ggplot2::guides(color = ggplot2::guide_legend(
                                override.aes = list(
                                    size = input$metadata_pt_size * 10,
                                    shape = 15
                                )
                            ))
                    } else {
                        ggplot_obj <- ggplot_obj +
                            ggplot2::scale_colour_gradientn(colours = viridis::viridis(50)) +
                            ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = grid::unit(input$width_metadata * 3 / 4, "inches")))
                    }

                    ggplot2::ggsave(
                        filename = file,
                        plot = ggplot_obj,
                        height = input$height_metadata,
                        width = input$width_metadata
                    )
                }
            )
        }
    )
}

server_comparison_gene_panel <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            gene_legend_height <- shiny::reactiveVal(0)
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

            changed_metadata <- shiny::reactiveVal(FALSE)

            shiny::observe({
                shiny::req(input$metadata_subset)
                is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")

                if (is_cluster) {
                    k_value <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = "metadata_groups_subset",
                        choices = seq_len(k_value),
                        selected = seq_len(k_value),
                        clearOptions = TRUE
                    )
                    changed_metadata(TRUE)

                    return()
                }

                mtd_names <- pkg_env$metadata_unique[[input$metadata_subset]]

                shinyWidgets::updatePickerInput(
                    session,
                    inputId = "metadata_groups_subset",
                    choices = mtd_names,
                    selected = mtd_names
                )

                changed_metadata(TRUE)
            }) %>% shiny::bindEvent(input$metadata_subset)

            shiny::observe({
                shiny::updateNumericInput(session,
                    inputId = "expr_threshold",
                    max = round(max_level_expr(), 3)
                    # step = round(max_level_expr() / 10, 3)
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
                    shiny::req(input$gene_expr, input$metadata_groups_subset, cancelOutput = TRUE)
                    plt_height()
                    relaxation <- input$relaxation
                    expr_threshold <- input$expr_threshold
                    input$gene_pt_type
                    input$gene_legend_size
                    input$gene_axis_size
                    input$gene_pt_size
                    input$gene_pt_order
                    input$metadata_groups_subset

                    shiny::isolate({
                        is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")
                        if (is_cluster) {
                            k <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                            all_unique_values <- as.character(seq_len(k))
                        } else {
                            all_unique_values <- pkg_env$metadata_unique[[input$metadata_subset]]
                        }

                        if (changed_metadata()) {
                            shiny::req(
                                all(input$metadata_groups_subset %in% all_unique_values),
                                cancelOutput = TRUE
                            )

                            changed_metadata(FALSE)
                        }

                        if (is_cluster) {
                            metadata_mask <- (pkg_env$stab_obj$mbs[[as.character(k)]] %in% as.integer(input$metadata_groups_subset))
                        } else {
                            metadata_mask <- (pkg_env$metadata[[input$metadata_subset]] %in% input$metadata_groups_subset)
                        }


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
                            embedding = pkg_env$stab_obj$umap,
                            color_info = used_matrix,
                            plt_height = plt_height(),
                            plt_width = plt_height(),
                            display_legend = FALSE,
                            cell_mask = metadata_mask,
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
                        shiny::req(input$metadata_groups_subset, cancelOutput = TRUE)
                        plt_height()
                        input$select_groups
                        expr_threshold <- input$expr_threshold
                        relaxation <- input$relaxation
                        input$gene_expr
                        input$gene_legend_size
                        input$metadata_groups_subset

                        shiny::isolate({
                            is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")
                            if (is_cluster) {
                                k <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                                all_unique_values <- as.character(seq_len(k))
                            } else {
                                all_unique_values <- pkg_env$metadata_unique[[input$metadata_subset]]
                            }

                            if (changed_metadata()) {
                                shiny::req(
                                    all(input$metadata_groups_subset %in% all_unique_values),
                                    cancelOutput = TRUE
                                )
                            }

                            if (is_cluster) {
                                metadata_mask <- (pkg_env$stab_obj$mbs[[as.character(k)]] %in% as.integer(input$metadata_groups_subset))
                            } else {
                                metadata_mask <- (pkg_env$metadata[[input$metadata_subset]] %in% input$metadata_groups_subset)
                            }

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
                                color_info = used_matrix[metadata_mask],
                                plt_width = plt_height(),
                                text_size = input$gene_legend_size
                            )
                        })
                    }
                )
            })

            output$download_gene <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_gene, ".", tolower(input$filetype_gene))
                },
                content = function(file) {
                    shiny::req(input$expr_threshold, input$gene_expr, input$width_gene, input$height_gene, expr_matrix())
                    is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")
                    if (is_cluster) {
                        k <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                        all_unique_values <- as.character(seq_len(k))
                    } else {
                        all_unique_values <- pkg_env$metadata_unique[[input$metadata_subset]]
                    }

                    if (changed_metadata()) {
                        shiny::req(
                            all(input$metadata_groups_subset %in% all_unique_values),
                            cancelOutput = TRUE
                        )
                    }

                    if (is_cluster) {
                        metadata_mask <- (pkg_env$stab_obj$mbs[[as.character(k)]] %in% as.integer(input$metadata_groups_subset))
                    } else {
                        metadata_mask <- (pkg_env$metadata[[input$metadata_subset]] %in% input$metadata_groups_subset)
                    }

                    # filetypes[[input$filetype_gene]](file, width = input$width_gene, height = input$height_gene)
                    unique_values <- NULL
                    used_matrix <- expr_matrix()
                    color_values <- function(n) {
                        grDevices::colorRampPalette(c("grey85", RColorBrewer::brewer.pal(9, "OrRd")))(n)
                    }
                    if (length(input$gene_expr) > 1) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("other" = "lightgray", "cells above threshold" = "red")
                        used_matrix <- factor(ifelse(matrixStats::colSums2(used_matrix > input$expr_threshold) >= (length(input$gene_expr) - input$relaxation), "cells above threshold", "other"))
                    } else if (input$expr_threshold > 0) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("other" = "lightgray", "cells above threshold" = "red")
                        used_matrix <- factor(ifelse(used_matrix > input$expr_threshold, "cells above threshold", "other"))
                    } else {
                        used_matrix <- as.numeric(used_matrix)
                    }


                    ggplot_obj <- color_ggplot(
                        embedding = pkg_env$stab_obj$umap,
                        color_info = used_matrix,
                        sort_cells = input$gene_pt_order,
                        cell_mask = metadata_mask,
                        pt_size = input$gene_pt_size
                    ) + ggplot2::ggtitle(paste(input$gene_expr, collapse = " ")) +
                        ggplot2::theme(
                            legend.position = "bottom",
                            legend.title = ggplot2::element_blank(),
                            legend.text = ggplot2::element_text(size = input$gene_legend_size * 10),
                            axis.text = ggplot2::element_text(size = input$gene_axis_size * 10),
                            axis.title = ggplot2::element_text(size = input$gene_axis_size * 10),
                            plot.title = ggtext::element_textbox_simple(hjust = 0.5, size = input$gene_axis_size * 10 * 1.5),
                            aspect.ratio = 1
                        )


                    if (length(input$gene_expr) > 1 || input$expr_threshold > 0) {
                        ggplot_obj <- ggplot_obj + ggplot2::scale_colour_manual(values = color_values)
                    } else {
                        ggplot_obj <- ggplot_obj +
                            ggplot2::scale_colour_gradientn(colours = color_values(50)) +
                            ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = grid::unit(input$width_gene * 3 / 4, "inches")))
                    }

                    ggplot2::ggsave(
                        filename = file,
                        plot = ggplot_obj,
                        height = input$height_gene,
                        width = input$width_gene
                    )
                }
            )
        }
    )
}

server_comparison_jsi <- function(id, k_choices) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::updateSelectInput(
                session = session,
                inputId = "jsi_k_1",
                choices = k_choices,
                selected = k_choices[1]
            )

            shiny::updateSelectInput(
                session = session,
                inputId = "jsi_k_2",
                choices = k_choices,
                selected = k_choices[1]
            )



            plt_height <- shiny::reactive(
                floor(pkg_env$height_ratio * pkg_env$dimension()[2])
            )

            barcode_heatmap <- shiny::reactive({
                shiny::req(input$jsi_k_1, input$jsi_k_2)
                if (!is.na(as.numeric(input$jsi_k_1))) {
                    clustering_1 <- as.matrix(pkg_env$stab_obj$mbs[[as.character(input$jsi_k_1)]])
                    df_1 <- data.frame(clustering_1)
                } else {
                    meta_category <- pkg_env$metadata[, input$jsi_k_1]
                    df_1 <- data.frame(meta_category)
                }
                df_1$cell <- rownames(df_1)
                if (!is.na(as.numeric(input$jsi_k_2))) {
                    clustering_2 <- as.matrix(pkg_env$stab_obj$mbs[[as.character(input$jsi_k_2)]])
                    df_2 <- data.frame(clustering_2)
                } else {
                    meta_category <- pkg_env$metadata[, input$jsi_k_2]
                    df_2 <- data.frame(meta_category)
                }
                df_2$cell <- rownames(df_2)
                all_clusters_1 <- unique(df_1[, 1])
                all_clusters_2 <- unique(df_2[, 1])

                mat <- matrix(,
                    nrow = length(all_clusters_2),
                    ncol = length(all_clusters_1)
                )
                colnames(mat) <- sort(all_clusters_1)
                rownames(mat) <- sort(all_clusters_2)
                if (input$heatmap_type == "JSI") {
                    for (m in all_clusters_1) {
                        cluster_1 <- rownames(df_1[df_1[, 1] == m, ])
                        for (n in all_clusters_2) {
                            cluster_2 <- rownames(df_2[df_2[, 1] == n, ])
                            mat[as.character(n), as.character(m)] <- jaccard_index(cluster_1, cluster_2)
                            label <- "JSI"
                        }
                    }
                } else {
                    for (m in all_clusters_1) {
                        cluster_1 <- rownames(df_1[df_1[, 1] == m, ])
                        for (n in all_clusters_2) {
                            cluster_2 <- rownames(df_2[df_2[, 1] == n, ])
                            mat[as.character(n), as.character(m)] <- length(intersect(cluster_1, cluster_2))
                            label <- "Shared cells"
                        }
                    }
                }
                df_mat <- reshape2::melt(mat)

                ggplot2::ggplot(df_mat, ggplot2::aes(as.factor(.data$Var1), as.factor(.data$Var2))) +
                    ggplot2::geom_tile(ggplot2::aes(fill = .data$value)) +
                    ggplot2::geom_text(ggplot2::aes(label = round(.data$value, 2))) +
                    ggplot2::scale_fill_gradient2(
                        low = scales::muted("darkred"),
                        mid = "white",
                        high = scales::muted("green"),
                        midpoint = 0
                    ) +
                    ggplot2::theme(
                        panel.background = ggplot2::element_rect(fill = "white"),
                        axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, size = 10, face = "bold"),
                        axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
                        axis.title = ggplot2::element_text(size = 14, face = "bold"),
                        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
                        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))
                    ) +
                    ggplot2::xlab("Configuration 2") +
                    ggplot2::ylab("Configuration 1") +
                    ggplot2::labs(fill = label)
            })

            output$barcode_heatmap <- shiny::renderPlot(
                {
                    if (is.null(input$jsi_k_1) | is.null(input$jsi_k_2)) {
                        return(ggplot2::ggplot() +
                            ggplot2::theme_void())
                    } else {
                        barcode_heatmap()
                    }
                },
                height = plt_height()
            )

            heatmap_filetype <- shiny::reactive({
                if (input$heatmap_filetype == "PDF") {
                    filename <- paste0(input$filename_heatmap, ".pdf")
                    return(filename)
                } else if (input$heatmap_filetype == "PNG") {
                    filename <- paste0(input$filename_heatmap, ".png")
                    return(filename)
                } else {
                    filename <- paste0(input$filename_heatmap, ".svg")
                    return(filename)
                }
            })

            output$download_heatmap <- shiny::downloadHandler(
                filename = function() {
                    heatmap_filetype()
                },
                content = function(file) {
                    ggplot2::ggsave(file, barcode_heatmap(),
                        width = input$width_heatmap,
                        height = input$height_heatmap,
                        units = "in",
                        limitsize = FALSE
                    )
                }
            )
        }
    )
}

server_comparison_violin_gene <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            changed_metadata <- shiny::reactiveVal(FALSE)

            shiny::observe({
                shiny::req(input$metadata_subset)
                is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")

                if (is_cluster) {
                    k_value <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = "metadata_groups_subset",
                        choices = seq_len(k_value),
                        selected = seq_len(k_value),
                        clearOptions = TRUE
                    )
                    changed_metadata(TRUE)

                    return()
                }

                mtd_names <- pkg_env$metadata_unique[[input$metadata_subset]]

                shinyWidgets::updatePickerInput(
                    session,
                    inputId = "metadata_groups_subset",
                    choices = mtd_names,
                    selected = mtd_names
                )

                changed_metadata(TRUE)
            }) %>% shiny::bindEvent(input$metadata_subset)

            metadata_mask <- shiny::reactive({
                shiny::req(input$metadata_groups_subset, input$metadata_subset, cancelOutput = TRUE)
                shiny::isolate({
                    is_cluster <- stringr::str_detect(input$metadata_subset, "stable_[0-9]+_clusters")
                    if (is_cluster) {
                        k <- as.numeric(strsplit(input$metadata_subset, "_")[[1]][2])
                        all_unique_values <- as.character(seq_len(k))
                    } else {
                        all_unique_values <- pkg_env$metadata_unique[[input$metadata_subset]]
                    }

                    if (changed_metadata()) {
                        shiny::req(
                            all(input$metadata_groups_subset %in% all_unique_values),
                            cancelOutput = TRUE
                        )

                        changed_metadata(FALSE)
                    }

                    if (is_cluster) {
                        return(pkg_env$stab_obj$mbs[[as.character(k)]] %in% as.integer(input$metadata_groups_subset))
                    }
                    return(pkg_env$metadata[[input$metadata_subset]] %in% input$metadata_groups_subset)
                })
            })


            distr_val <- shiny::reactive({
                shiny::req(input$gene_expr, length(input$gene_expr) == 1, cancelOutput = TRUE)

                is_ecc <- stringr::str_detect(input$gene_expr, "ecc_[0-9]+")
                is_continuous <- (!is_ecc && !(input$gene_expr %in% names(pkg_env$metadata_unique)) && (input$gene_expr %in% colnames(pkg_env$metadata)))

                if (is_ecc) {
                    k <- strsplit(input$gene_expr, "_")[[1]][2]
                    cl_method <- strsplit(names(pkg_env$stab_obj$ecc)[1], ";")[[1]][2]
                    temp_val <- pkg_env$stab_obj$ecc[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]
                    return(temp_val[pkg_env$stab_obj$ecc_order[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]])
                }

                if (is_continuous) {
                    return(pkg_env$metadata[[input$gene_expr]])
                }

                if ("genes" %in% names(pkg_env)) {
                    index_gene <- pkg_env$genes[input$gene_expr]
                    return(rhdf5::h5read("expression.h5", "expression_matrix", index = list(index_gene, NULL))[1, ])
                }

                # for backward-compatibility purposes
                index_interest <- pkg_env$genes_of_interest[input$gene_expr]
                index_others <- pkg_env$genes_others[input$gene_expr]

                if (is.na(index_interest)) {
                    return(rhdf5::h5read("expression.h5", "matrix_others", index = list(index_others, NULL))[1, ])
                }
                return(rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(index_interest, NULL))[1, ])
            })

            metadata_split_info <- shiny::reactive({
                metadata_split <- input$metadata_split
                shiny::req(metadata_split, cancelOutput = TRUE)
                is_cluster <- stringr::str_detect(metadata_split, "stable_[0-9]+_clusters")

                if (is_cluster) {
                    k <- strsplit(metadata_split, "_")[[1]][2]
                    return(list(
                        color_values = rhdf5::h5read("stability.h5", paste0("colors/", k)),
                        color_info = factor(pkg_env$stab_obj$mbs[[k]]),
                        unique_values = as.character(seq_len(as.integer(k)))
                    ))
                }

                return(list(
                    color_values = pkg_env$metadata_colors[[metadata_split]],
                    color_info = pkg_env$metadata[[metadata_split]],
                    unique_values = pkg_env$metadata_unique[[metadata_split]]
                ))
            })

            metadata_group_info <- shiny::reactive({
                metadata_group <- input$metadata_group
                shiny::req(metadata_group, cancelOutput = TRUE)
                is_cluster <- stringr::str_detect(metadata_group, "stable_[0-9]+_clusters")

                if (is_cluster) {
                    k <- strsplit(metadata_group, "_")[[1]][2]
                    return(list(
                        color_values = rhdf5::h5read("stability.h5", paste0("colors/", k)),
                        color_info = factor(pkg_env$stab_obj$mbs[[k]]),
                        unique_values = as.character(seq_len(as.integer(k)))
                    ))
                }

                return(list(
                    color_values = pkg_env$metadata_colors[[metadata_group]],
                    color_info = pkg_env$metadata[[metadata_group]],
                    unique_values = pkg_env$metadata_unique[[metadata_group]]
                ))
            })

            shiny::observe({
                mtd_grp_info <- metadata_group_info()
                shiny::req(mtd_grp_info)
                unique_values <- mtd_grp_info$unique_values

                shiny::updateSelectInput(
                    session = session,
                    inputId = "stat_mtd_group",
                    choices = unique_values,
                    selected = unique_values[1]
                )
            })

            ggplot_object <- shiny::reactive({
                shiny::req(distr_val(), metadata_split_info(), metadata_group_info(), cancelOutput = TRUE)
                mtd_mask <- metadata_mask()
                graph_types <- input$graph_type
                if (is.null(graph_types)) {
                    graph_types <- "Violin"
                }
                boxplot_dodge <- input$boxplot_dodge

                is_ecc <- stringr::str_detect(input$gene_expr, "ecc_[0-9]+")
                is_continuous <- (!is_ecc && !(input$gene_expr %in% names(pkg_env$metadata_unique)) && (input$gene_expr %in% colnames(pkg_env$metadata)))
                # is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")

                function_applied <- ifelse(!input$log_scale,
                    function(x) {
                        x
                    },
                    function(x) {
                        log10(x)
                    }
                )

                df <- data.frame(
                    gene_expr = function_applied(distr_val()),
                    metadata_split = metadata_split_info()$color_info,
                    metadata_group = metadata_group_info()$color_info
                )
                df$metadata_group <- factor(df$metadata_group, levels = unique(df$metadata_group))
                df$metadata_split <- factor(df$metadata_split, levels = unique(df$metadata_split))

                df <- df[mtd_mask, ]

                gplot_object <- ggplot2::ggplot(
                    df,
                    ggplot2::aes(x = .data$metadata_split, y = .data$gene_expr, fill = .data$metadata_group)
                ) +
                    ggplot2::scale_fill_manual(values = metadata_group_info()$color_values, name = input$metadata_group) +
                    # ggplot2::geom_violin(width = input$boxplot_width) +
                    # ggplot2::geom_boxplot(width = input$boxplot_width, outlier.shape = NA) +
                    ggplot2::theme(
                        axis.text = ggplot2::element_text(size = input$text_size),
                        axis.title = ggplot2::element_text(size = input$text_size)
                    ) +
                    ggplot2::ylab(paste0(ifelse(
                        is_ecc,
                        "ECC",
                        ifelse(
                            is_continuous,
                            input$gene_expr,
                            "Gene Expression"
                        )
                    ), ifelse(input$log_scale, " (log10 scale)", ""))) +
                    ggplot2::xlab(input$metadata)

                # FIXME empty distribution cause shifts in violin plots and mismatch with the boxplots
                # TODO add violin plots for multiple genes
                if ("Violin" %in% graph_types) {
                    gplot_object <- gplot_object + ggplot2::geom_violin(
                        width = input$boxplot_width,
                        trim = FALSE,
                        position = ggplot2::position_dodge(width = boxplot_dodge)
                    )
                }

                if ("Boxplot" %in% graph_types) {
                    gplot_object <- gplot_object + ggplot2::geom_boxplot(
                        width = input$boxplot_width,
                        position = ggplot2::position_dodge(width = boxplot_dodge),
                        outlier.shape = NA
                    )
                }

                gplot_object
            })

            output$violin_gene <- shiny::renderPlot(
                {
                    shiny::req(ggplot_object(), cancelOutput = TRUE)
                    ggplot_object()
                },
                height = function() {
                    pkg_env$plt_height()
                }
            )

            output$download_violin <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_violin, ".", tolower(input$filetype_violin))
                },
                content = function(file) {
                    shiny::req(ggplot_object())

                    ggplot2::ggsave(
                        filename = file,
                        plot = ggplot_object() + ggplot2::ggtitle(
                            glue::glue("Distribution of {input$gene_expr} - Split by {input$metadata}")
                        ),
                        height = input$height_violin,
                        width = input$width_violin
                    )
                }
            )

            output$stats <- shiny::renderTable(
                {
                    mtd_mask <- metadata_mask()
                    shiny::req(mtd_mask)
                    distr_vector <- distr_val()[mtd_mask]
                    mtd_split_info <- metadata_split_info()
                    mtd_split_info$color_info <- mtd_split_info$color_info[mtd_mask]
                    mtd_group_info <- metadata_group_info()
                    mtd_group_info$color_info <- mtd_group_info$color_info[mtd_mask]
                    selected_group <- input$stat_mtd_group
                    shiny::req(
                        distr_vector,
                        mtd_split_info,
                        selected_group,
                        mtd_group_info,
                        selected_group %in% mtd_group_info$unique_values,
                        cancelOutput = TRUE
                    )

                    is_ecc <- stringr::str_detect(input$gene_expr, "ecc_[0-9]+")
                    is_continuous <- (!is_ecc && !(input$gene_expr %in% names(pkg_env$metadata_unique)) && (input$gene_expr %in% colnames(pkg_env$metadata)))

                    if (length(mtd_group_info$unique_values) > 1) {
                        mask <- mtd_group_info$color_info == selected_group
                        distr_vector <- distr_vector[mask]
                        mtd_split_info$color_info <- mtd_split_info$color_info[mask]
                    }

                    distr_stats <- stats::fivenum(distr_vector)
                    distance_breaks <- (distr_stats[5] - distr_stats[1]) / 4
                    break_points <- c(
                        distr_stats[1],
                        distr_stats[1] + distance_breaks,
                        distr_stats[1] + 2 * distance_breaks,
                        distr_stats[1] + 3 * distance_breaks,
                        distr_stats[5]
                    )

                    split_vals <- split(
                        distr_vector,
                        mtd_split_info$color_info
                    )
                    for (i in names(split_vals)) {
                        if (is.null(split_vals[[i]]) || length(split_vals[[i]]) == 0) {
                            split_vals[[i]] <- NULL
                        }
                    }

                    stats_df <- rbind(
                        data.frame(sapply(seq_along(split_vals), function(i) {
                            stats::fivenum(split_vals[[i]])
                        })),
                        sapply(split_vals, length)
                    )

                    stats_df <- rbind(
                        stats_df,
                        sapply(seq_along(split_vals), function(x) { sum(split_vals[[x]] > stats_df[1, x])}),
                        sapply(seq_along(split_vals), function(x) { sum(split_vals[[x]] < stats_df[5, x])})
                    )

                    colnames(stats_df) <- names(split_vals)
                    rownames(stats_df) <- c("Min", "Q1", "Median", "Q3", "Max", "# cells", "# cells above min", "# cells under max")

                    breaks_df <- sapply(split_vals, function(x) {
                        table(cut(x, breaks = break_points))
                    })

                    rbind(stats_df, breaks_df)
                },
                rownames = TRUE
            )
        }
    )
}

server_comparison_gene_heatmap <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            heatmap_plot <- shiny::reactive({
                shiny::req(input$gene_expr, length(input$gene_expr) > 0, input$metadata, input$text_size, !is.na(input$scale), input$clipping_value, !is.na(input$show_numbers), input$plot_type, input$point_size, cancelOutput = TRUE)

                shiny::isolate({
                    is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")

                    if (is_cluster) {
                        k <- strsplit(input$metadata, "_")[[1]][2]
                        mtd_val <- pkg_env$stab_obj$mbs[[k]]
                        unique_vals <- as.character(seq_len(as.numeric(k)))
                    } else {
                        mtd_val <- pkg_env$metadata[[input$metadata]]
                        unique_vals <- levels(mtd_val)
                    }

                    htmp_matrix <- matrix(0, nrow = length(input$gene_expr), ncol = length(unique_vals))
                    perc_expressed <- matrix(0, nrow = length(input$gene_expr), ncol = length(unique_vals))
                    rownames(htmp_matrix) <- input$gene_expr
                    rownames(perc_expressed) <- input$gene_expr
                    colnames(htmp_matrix) <- unique_vals
                    colnames(perc_expressed) <- unique_vals

                    if ("genes" %in% names(pkg_env)) {
                        index_gene <- pkg_env$genes[input$gene_expr]

                        for (i in seq_along(input$gene_expr)) {
                            expr_profile <- rhdf5::h5read("expression.h5", "expression_matrix", index = list(index_gene[i], NULL))
                            for (j in seq_along(unique_vals)) {
                                filtered_expr <- expr_profile[mtd_val == unique_vals[j]]
                                htmp_matrix[i, j] <- mean(filtered_expr, na.rm = TRUE)
                                perc_expressed[i, j] <- sum(filtered_expr > 0) / length(filtered_expr)
                            }
                        }
                    } else { # for backward-compatibility purposes
                        index_interest <- pkg_env$genes_of_interest[input$gene_expr]
                        index_others <- pkg_env$genes_others[input$gene_expr]

                        for (i in seq_along(index_interest)) {
                            expr_profile <- rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(index_interest[i], NULL))
                            for (j in seq_along(unique_vals)) {
                                filtered_expr <- expr_profile[mtd_val == unique_vals[j]]
                                htmp_matrix[i, j] <- mean(filtered_expr, na.rm = TRUE)
                                perc_expressed[i, j] <- sum(filtered_expr > 0) / length(filtered_expr)
                            }
                        }

                        offset <- length(index_interest)
                        for (i in seq_along(index_others)) {
                            expr_profile <- rhdf5::h5read("expression.h5", "matrix_others", index = list(index_others[i], NULL))
                            for (j in seq_along(unique_vals)) {
                                filtered_expr <- expr_profile[mtd_val == unique_vals[j]]
                                htmp_matrix[i + offset, j] <- mean(filtered_expr, na.rm = TRUE)
                                perc_expressed[i, j] <- sum(filtered_expr > 0) / length(filtered_expr)
                            }
                        }
                    }

                    nelems <- nrow(htmp_matrix) * ncol(htmp_matrix)

                    if (input$scale && ncol(htmp_matrix) > 1) {
                        htmp_matrix <- t(scale(t(htmp_matrix)))
                        original_htmp_matrix <- htmp_matrix
                        htmp_matrix[htmp_matrix > input$clipping_value] <- input$clipping_value
                        htmp_matrix[htmp_matrix < -input$clipping_value] <- -input$clipping_value

                        colour_scheme <- colorRampPalette(c("blue", "white", "red"))(nelems*2)
                    } else {
                        original_htmp_matrix <- htmp_matrix
                        htmp_matrix[htmp_matrix > input$clipping_value] <- input$clipping_value
                        colour_scheme <- colorRampPalette(c("grey85", "#004c00"))(nelems*2)
                    }

                    if (min(htmp_matrix) == max(htmp_matrix)) {
                        colour_scheme <- "grey85"
                    }

                    if (input$plot_type == "Heatmap") {
                        return(ComplexHeatmap::Heatmap(
                            htmp_matrix,
                            row_order = seq_len(nrow(htmp_matrix)),
                            column_order = seq_len(ncol(htmp_matrix)),
                            row_names_side = "left",
                            heatmap_legend_param = list(direction = "horizontal", legend_width = grid::unit(5,  "cm")),
                            name = paste0(ifelse(input$scale, "scaled ", ""), "expression level"),
                            col = colour_scheme,
                            cell_fun = function(j, i, x, y, width, height, fill) {
                                if (input$show_numbers) {
                                    grid::grid.text(sprintf("%.2f", original_htmp_matrix[i, j]), x, y, just = "center", gp = grid::gpar(fontsize = input$text_size))
                                }
                            },
                            column_title = paste0("Gene expression heatmap split by ", input$metadata)
                        ))
                    }

                    bubbleplot_df <- reshape2::melt(htmp_matrix)
                    colnames(bubbleplot_df) <- c("gene", "metadata", "expr_level")
                    bubbleplot_df$perc <- reshape2::melt(perc_expressed)$value
                    bubbleplot_df$gene <- factor(bubbleplot_df$gene, levels = input$gene_expr)
                    bubbleplot_df$metadata <- factor(bubbleplot_df$metadata, levels = unique_vals)

                    ggplot2::ggplot(bubbleplot_df, ggplot2::aes(x = .data$metadata, y = .data$gene, size = .data$perc, fill = .data$expr_level)) +
                        ggplot2::geom_point(shape = 21, alpha = 0.7) +
                        ggplot2::scale_fill_gradientn(colours = colour_scheme) +
                        ggplot2::scale_size_continuous(range = input$point_size) +
                        ggplot2::theme(
                            axis.text.x = ggplot2::element_text(size = input$text_size),
                            axis.text.y = ggplot2::element_text(size = input$text_size),
                            axis.title = ggplot2::element_text(size = input$text_size),
                            legend.position = "bottom"
                        ) +
                        ggplot2::xlab("") +
                        ggplot2::ylab("")




                })
            })

            shiny::observe({
                shiny::req(input$gene_expr, heatmap_plot())
                shiny::req(pkg_env$dimension())
                
                shiny::isolate({
                    output$gene_heatmap <- shiny::renderPlot(
                        width = pkg_env$dimension()[1],
                        height = 200 + length(input$gene_expr) * 50,
                        {
                            if (input$plot_type == "Bubbleplot") {
                                heatmap_plot()
                            } else {
                                ComplexHeatmap::draw(heatmap_plot(), heatmap_legend_side = "bottom")
                            }
                        }
                    )

                    output$download_heatmap <- shiny::downloadHandler(
                        filename = function() {
                            paste0(input$filename_heatmap, ".", tolower(input$filetype_heatmap))
                        },
                        content = function(file) {
                            shiny::req(heatmap_plot())

                            if (input$plot_type == "Bubbleplot") {
                                ggplot2::ggsave(
                                    filename = file,
                                    plot = heatmap_plot(),
                                    height = input$height_heatmap,
                                    width = input$width_heatmap
                                )
                                return()
                            }

                            pdf(file, width = input$width_heatmap, height = input$height_heatmap)
                            ComplexHeatmap::draw(heatmap_plot(), heatmap_legend_side = "bottom")
                            dev.off()
                        }
                    )
                })
            })
            

        }
    )
}

server_comparison_enrichment <- function(id, marker_genes) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            # isolated_marker_genes <- shiny::isolate(marker_genes())
            shiny::observe({
                shiny::req(marker_genes())
                shinyjs::show("enrichment_id")
                shinyjs::hide("download_gost")
            }) %>% shiny::bindEvent(marker_genes(), once = TRUE)

            gprof_result <- shiny::reactive({
                shiny::req(marker_genes())
                shinyjs::disable("enrichment_button")
                selected_group <- stringr::str_replace(input$group, " ", "_")

                # TODO use as background only the genes with avg expression over a specific
                if ("genes" %in% names(pkg_env)) {
                    custom_bg <- names(pkg_env$genes)
                } else { # backward-compatibility
                    custom_bg <- c(names(pkg_env$genes_others), names(pkg_env$genes_of_interest))
                }

                gprf_res <- gprofiler2::gost(
                    query = marker_genes()[[selected_group]],
                    sources = input$gprofilerSources,
                    organism = pkg_env$organism,
                    evcodes = TRUE,
                    domain_scope = "custom",
                    custom_bg = custom_bg
                )

                if (!is.null(gprf_res)) {
                    gprf_res$result$parents <- sapply(gprf_res$result$parents, toString)
                }

                shinyjs::enable("enrichment_button")
                shinyjs::show("download_gost")

                gprf_res
            }) %>% shiny::bindEvent(input$enrichment_button)

            shiny::observe({
                shiny::req(gprof_result())

                output$gost_table <- DT::renderDT({
                    gprof_result()$result[, seq_len(ncol(gprof_result()$result) - 2)]
                })

                output$gost_plot <- plotly::renderPlotly(
                    gprofiler2::gostplot(gprof_result())
                )

                output$download_gost <- shiny::downloadHandler(
                    filename = function() {
                        "enrichment_results.csv"
                    },
                    content = function(file) {
                        utils::write.csv(gprof_result()$result, file)
                    }
                )
            })
        }
    )
}

#' Server - Comparison module
#'
#' @description Creates the backend interface for the comparison module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to acess the UI elements.
#' @param chosen_config A reactive object that contains the chosen configuration
#' from the Dimensionality Reduction tab.
#' @param chosen_method A reactive object that contains the chosen method from
#' the Clustering tab.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
server_comparisons <- function(id, chosen_config, chosen_method) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shinyjs::hide("enrichment-enrichment_id")
            isolated_chosen_config <- shiny::isolate(chosen_config())
            ftype <- isolated_chosen_config$chosen_feature_type
            fsize <- isolated_chosen_config$chosen_set_size

            isolated_chosen_method <- shiny::isolate(chosen_method())
            cl_method <- isolated_chosen_method$method_name

            discrete <- c()
            for (category in colnames(pkg_env$metadata)) {
                if (length(unique(pkg_env$metadata[, category])) < 20) {
                    discrete <- append(discrete, category)
                }
            }
            k_values <- isolated_chosen_method$n_clusters
            stable_config <- rhdf5::h5read("stability.h5", paste(ftype, fsize, "stable_config", sep = "/"))

            add_env_variable("stab_obj", list(
                mbs = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "mbs", cl_method, sep = "/"))[k_values],
                ecc = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "ecc", sep = "/"))[paste(sprintf("%06d", as.integer(k_values)), cl_method, sep = ";")],
                ecc_order = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "ecc_order", sep = "/"))[paste(sprintf("%06d", as.integer(k_values)), cl_method, sep = ";")],
                umap = rhdf5::h5read("stability.h5", paste(ftype, fsize, "umap", sep = "/"))
            ))

            add_env_variable("used_genes", rhdf5::h5read("stability.h5", paste(ftype, "feature_list", sep = "/"))[seq_len(as.numeric(fsize))])
            add_env_variable("current_tab", "Comparison")

            for (panels in c("left", "right")) {
                shiny::updateSelectizeInput(
                    session,
                    inputId = glue::glue("metadata_panel_{panels}-metadata"),
                    server = FALSE,
                    choices = c(colnames(pkg_env$metadata), paste0("stable_", k_values, "_clusters"), paste0("ecc_", k_values)),
                    selected = paste0("stable_", k_values[1], "_clusters")
                )

                shiny::updateSelectizeInput(
                    session,
                    inputId = glue::glue("gene_panel_{panels}-metadata_subset"),
                    server = FALSE,
                    choices = c(names(pkg_env$metadata_unique), paste0("stable_", k_values, "_clusters")),
                    selected = paste0("stable_", k_values[1], "_clusters")
                )

                shiny::updateSelectizeInput(
                    session,
                    inputId = glue::glue("metadata_panel_{panels}-metadata_subset"),
                    server = FALSE,
                    choices = c(names(pkg_env$metadata_unique), paste0("stable_", k_values, "_clusters")),
                    selected = paste0("stable_", k_values[1], "_clusters")
                )

                if ("genes" %in% names(pkg_env)) {
                    gene_choices <- names(pkg_env$genes)
                } else { # for backward-compatibility purposes
                    gene_choices <- c(names(pkg_env$genes_of_interest), names(pkg_env$genes_others))
                }

                shiny::updateSelectizeInput(
                    session,
                    inputId = glue::glue("gene_panel_{panels}-gene_expr"),
                    choices = gene_choices,
                    selected = gene_choices[1],
                    server = TRUE,
                    options = list(
                        maxOptions = 7,
                        create = TRUE
                    )
                )
            }

            shiny::updateSelectizeInput(
                session,
                inputId = glue::glue("violin_gene-metadata_split"),
                server = FALSE,
                choices = c(names(pkg_env$metadata_unique), paste0("stable_", k_values, "_clusters")),
                selected = paste0("stable_", k_values[1], "_clusters")
            )

            shiny::updateSelectizeInput(
                session,
                inputId = glue::glue("violin_gene-metadata_group"),
                server = FALSE,
                choices = c(names(pkg_env$metadata_unique), paste0("stable_", k_values, "_clusters")),
                selected = "one_level"
            )

            shiny::updateSelectizeInput(
                session,
                inputId = glue::glue("violin_gene-metadata_subset"),
                server = FALSE,
                choices = c(names(pkg_env$metadata_unique), paste0("stable_", k_values, "_clusters")),
                selected = "one_level"
            )

            shiny::updateSelectizeInput(
                session,
                inputId = glue::glue("gene_heatmap-metadata"),
                server = FALSE,
                choices = c(names(pkg_env$metadata_unique), paste0("stable_", k_values, "_clusters")),
                selected = paste0("stable_", k_values[1], "_clusters")
            )

            continuous_metadata <- c(
                setdiff(colnames(pkg_env$metadata), names(pkg_env$metadata_unique)),
                paste0("ecc_", k_values)
            )

            shiny::updateSelectizeInput(
                session,
                inputId = glue::glue("violin_gene-gene_expr"),
                choices = c(continuous_metadata, gene_choices),
                selected = continuous_metadata[1],
                server = TRUE,
                options = list(
                    maxOptions = length(continuous_metadata),
                    create = TRUE
                )
            )

            shiny::updateSelectizeInput(
                session,
                inputId = glue::glue("gene_heatmap-gene_expr"),
                choices = gene_choices,
                selected = gene_choices[1],
                server = TRUE,
                options = list(
                    maxOptions = 7,
                    create = TRUE
                )
            )

            server_comparison_metadata_panel("metadata_panel_left")
            server_comparison_metadata_panel("metadata_panel_right")
            server_comparison_gene_panel("gene_panel_left")
            server_comparison_gene_panel("gene_panel_right")
            server_comparison_jsi("jsi_plot", append(k_values, discrete))
            server_comparison_violin_gene("violin_gene")
            server_comparison_gene_heatmap("gene_heatmap")
            marker_genes <- server_comparison_markers("markers", k_values)
            server_comparison_enrichment("enrichment", marker_genes)

            shiny::observe({
                shiny::showModal(
                    stable_config_info(stable_config),
                    session = session
                )
            }) %>% shiny::bindEvent(input$show_config, ignoreInit = TRUE)

            shiny::observe(compar_info(session)) %>% shiny::bindEvent(input$info_title, ignoreInit = TRUE)

            shiny::observe({
                gc()
            }) %>% shiny::bindEvent(input$"metadata_panel_right-metadata_groups_subset", once = TRUE)
        }
    )
}
