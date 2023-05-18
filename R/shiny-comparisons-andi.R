####### UI #######

ui_comparison_markers <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::h2("Identification of markers"),
    shinyWidgets::dropdownButton(
      shiny::tagList(
        shiny::sliderInput(
          inputId = ns("logfc"),
          label = "logFC threshold",
          min = 0.00, max = 10.00, value = 0.50, step = 0.1
        ),
        shiny::sliderInput(
          inputId = ns("min_pct"),
          label = "Minimum gene frequency",
          min = 0.01, max = 1.00, value = 0.10, step = 0.01
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
    shiny::actionButton(ns("markers_button"), "Find markers!"),
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
      # width = "100%",
      # width = "30%",
      options = list(
        `actions-box` = TRUE,
        title = "Select/deselect subgroups",
        # actionsBox = TRUE,
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
    shiny::splitLayout(
      cellWidths = "50%",
      shiny::selectizeInput(
        inputId = ns("metadata"),
        label = "Metadata",
        choices = NULL
      ),
      shiny::verticalLayout(
        shiny::tags$b("Select the highlighted groups"),
        shinyWidgets::pickerInput(
          inputId = ns("select_groups"),
          # label = "Select the highlighted groups",
          choices = "",
          # width = "100%",
          # width = "30%",
          options = list(
            `actions-box` = TRUE,
            title = "Select/deselect groups",
            # actionsBox = TRUE,
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
    shiny::verticalLayout(
      shiny::selectizeInput(
        inputId = ns("gene_expr"),
        choices = NULL,
        label = "Gene name(s)",
        width = "95%",
        multiple = TRUE
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
      )
    ),
    shiny::splitLayout(
      cellWidths = c("40px", "40px"),
      gear_umaps(ns, "gene"),
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
        shiny::tagList("", a("https://github.com/Core-Bioinformatics/ClustAssess", href = "https://github.com/Core-Bioinformatics/ClustAssess", target = "_blank")),
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
      label = "Select the number of clusters (k) for the first comparison",
      choices = ""
    ),
    shiny::selectInput(
      inputId = ns("jsi_k_2"),
      label = "Select the number of clusters (k) for the second comparison",
      choices = ""
    ),
    shiny::plotOutput(ns("barcode_heatmap"), height = "auto", width = "98%")
  )
}

#' Writing objects
#'
#' @description to be completed
#'
#' @export
ui_comparisons <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Comparison",
    shiny::fluidRow(
      shiny::h1("Compare your current configuration", style = "margin-bottom:10px ")
    ),
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
    ui_comparison_markers(ns("markers")),
    style = "margin-left: 25px;margin-top:72px;"
  )
}
####### SERVER #######
server_comparison_markers <- function(id, k_choices) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
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

      server_comparison_markers_panels("group_left", k_choices)
      server_comparison_markers_panels("group_right", k_choices)

      shiny::observe({
        current_button_value <- as.integer(shiny::isolate(input$enable_markers))
        shiny::req(pkg_env$enable_markers_button() != current_button_value)
        pkg_env$enable_markers_button(current_button_value)
        # print(pkg_env$pressed_button)
        shinyjs::html("marker_text", "Preparing the objects for the analysis...")

        expr_matrix <- rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(pkg_env$genes_of_interest[pkg_env$used_genes], NULL))
        rownames(expr_matrix) <- pkg_env$used_genes
        # TODO check if you really need to read the rank matrix i.e. check if calculating the rank matrix on the fly would be faster
        add_env_variable("rank_matrix", rhdf5::h5read("expression.h5", "rank_of_interest", index = list(pkg_env$genes_of_interest[pkg_env$used_genes], NULL)))
        add_env_variable("expr_matrix", expr_matrix)
        # waiter::waiter_hide()
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

        markers_result <- calculate_markers(
          expression_matrix = pkg_env$expr_matrix, # expression matrix
          cells1 = cells_index_left,
          cells2 = cells_index_right,
          rank_matrix = pkg_env$rank_matrix, # rank matrix
          norm_method = ifelse(input$norm_type, "LogNormalize", ""),
          min_pct_threshold = input$min_pct,
          logfc_threshold = input$logfc
        )
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
          write.csv(markers_val(), file)
        }
      )


      shiny::observe({
        shiny::req(markers_val())
        shinyjs::show("markers_download_button")
      }) %>% shiny::bindEvent(markers_val())
    }
  )
}

server_comparison_markers_panels <- function(id, k_choices) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      # also add categorical metadata
      available_choices <- c(names(pkg_env$metadata_unique), k_choices)

      shiny::updateSelectInput(
        session = session,
        inputId = "select_k_markers",
        choices = available_choices,
        selected = available_choices[1]
      )

      shiny::observe({
        shiny::req(input$select_k_markers %in% available_choices)

        if (is.na(as.numeric(input$select_k_markers))) {
          available_subgroups <- pkg_env$metadata_unique[[input$select_k_markers]]
        } else {
          available_subgroups <- seq_len(as.numeric(input$select_k_markers))
        }

        shinyWidgets::updatePickerInput(
          session = session,
          inputId = "select_clusters_markers",
          choices = available_subgroups,
          selected = available_subgroups[1]
        )
      }) %>% shiny::bindEvent(input$select_k_markers)
    }
  )
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
        shiny::req(input$metadata)
        is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
        is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

        if (is_cluster || is_ecc) {
          if (is_cluster) {
            k_value <- as.numeric(strsplit(input$metadata, "_")[[1]][2])
            shinyjs::show(id = "select_groups")
            # shiny::freezeReactiveValue(input, "select_groups")
            shinyWidgets::updatePickerInput(
              session,
              inputId = "select_groups",
              choices = seq_len(k_value),
              selected = seq_len(k_value),
              clearOptions = TRUE
            )
          } else {
            shinyjs::hide(id = "select_groups")
          }

          changed_metadata(TRUE)

          return()
        }

        mtd_names <- pkg_env$metadata_unique[[input$metadata]]
        if (is.null(mtd_names)) {
          shinyjs::hide(id = "select_groups")
        } else {
          shinyjs::show(id = "select_groups")
          # shiny::freezeReactiveValue(input, "select_groups")
          shinyWidgets::updatePickerInput(
            session,
            inputId = "select_groups",
            choices = mtd_names,
            selected = mtd_names
          )
        }

        changed_metadata(TRUE)
      }) %>% shiny::bindEvent(input$metadata)

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
            color_info <- pkg_env$stab_obj$mbs[[k]]
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
          plot_data()
          input$select_groups
          input$metadata_pt_size
          input$metadata_axis_size
          input$metadata_text_size
          input$metadata_labels
          input$metadata_pt_type
          input$metadata_legend_size
          plt_height()

          shiny::isolate({
            if (changed_metadata() && !is.null(plot_data()$unique_values)) {
              matched_elems <- match(input$select_groups, plot_data()$unique_values)
              matched_elems <- matched_elems[!is.na(matched_elems)]
              shiny::req(
                length(matched_elems) == length(plot_data()$unique_values),
                length(matched_elems) == length(input$select_groups),
                cancelOutput = TRUE
              )
              changed_metadata(FALSE)
            }

            if (is.null(plot_data()$unique_values)) {
              old_par <- par(mai = c(0.1, 0, 0.1, 0))
              text_height <- strheight("TE\nXT\n", units = "inches", cex = input$metadata_legend_size)
            } else {
              old_par <- par(mar = c(0, 0, 0, 0))
              predicted_width <- strwidth(c(" ", plot_data()$unique_values), units = "inches", cex = input$metadata_legend_size) * ppi
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

              text_height <- strheight(
                paste(
                  rep("TEXT", number_rows + 1),
                  collapse = "\n"
                ),
                units = "inches",
                cex = input$metadata_legend_size
              )
            }
            par(old_par)
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
              text_size = input$metadata_text_size,
              axis_size = input$metadata_axis_size,
              labels = input$metadata_labels,
              groups_highlight = input$select_groups
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
            plot_data()
            plt_height()
            input$select_groups
            input$metadata_legend_size

            shiny::isolate({
              if (!is.null(plot_data()$unique_values)) {
                matched_elems <- match(input$select_groups, plot_data()$unique_values)
                matched_elems <- matched_elems[!is.na(matched_elems)]
                if (changed_metadata()) {
                  shiny::req(
                    length(matched_elems) == length(plot_data()$unique_values),
                    length(matched_elems) == length(input$select_groups),
                    cancelOutput = TRUE
                  )
                }

                color_values <- plot_data()$color_values[matched_elems]
                unique_values <- plot_data()$unique_values[matched_elems]
                # changed_metadata(FALSE)
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
          filetypes[[input$filetype_metadata]](file, width = input$width_metadata, height = input$height_metadata)
          color_plot2(
            embedding = pkg_env$stab_obj$umap,
            color_info = plot_data()$color_info,
            color_values = plot_data()$color_values,
            unique_values = plot_data()$unique_values,
            plt_height = input$height_metadata * ppi, # - metadata_legend_height(),
            plt_width = input$width_metadata * ppi,
            # predicted_height = (metadata_legend_height() - 1) / ppi,
            pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
            pt_size = input$metadata_pt_size,
            text_size = input$metadata_text_size,
            axis_size = input$metadata_axis_size,
            legend_text_size = input$metadata_legend_size,
            labels = input$metadata_labels,
            groups_highlight = input$select_groups,
            display_legend = TRUE
          )
          dev.off()
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
          shiny::req(input$gene_expr, cancelOutput = TRUE)
          plt_height()
          relaxation <- input$relaxation
          expr_threshold <- input$expr_threshold
          input$gene_pt_type
          input$gene_legend_size
          input$gene_axis_size
          input$gene_pt_size

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
              grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n)
            }
            if (length(input$gene_expr) > 1) {
              unique_values <- c("other", "cells above threshold")
              color_values <- c("lightgray", "red")
              used_matrix <- matrixStats::colSums2(used_matrix > expr_threshold) >= (length(input$gene_expr) - relaxation)
            } else if (expr_threshold > 0) {
              unique_values <- c("other", "cells above threshold")
              color_values <- c("lightgray", "red")
              used_matrix <- used_matrix > expr_threshold
            }

          old_par <- par(mai = c(0.1, 0, 0.1, 0))
          text_height <- strheight("TE\nXT\n", units = "inches", cex = input$gene_legend_size)
          par(old_par)
          gene_legend_height(text_height * ppi)

          color_plot2(
            embedding = pkg_env$stab_obj$umap, 
            color_info = used_matrix,
            plt_height = plt_height(),
            plt_width = plt_height(),
            display_legend = FALSE,
            unique_values = unique_values,
            color_values = color_values,
            pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
            pt_size = input$gene_pt_size,
            axis = input$gene_axis_size,
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
              color_values <- function(n) { grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n) }
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

      output$download_gene <- shiny::downloadHandler(
        filename = function() {
          paste0(input$filename_gene, ".", tolower(input$filetype_gene))
        },
        content = function(file) {
          shiny::req(input$expr_threshold, input$gene_expr, input$width_gene, input$height_gene, expr_matrix())
          filetypes[[input$filetype_gene]](file, width = input$width_gene, height = input$height_gene)
          unique_values <- NULL
          used_matrix <- expr_matrix()
          color_values <- function(n) {
            grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n)
          }
          if (length(input$gene_expr) > 1) {
            unique_values <- c("other", "cells above threshold")
            color_values <- c("lightgray", "red")
            used_matrix <- matrixStats::colSums2(used_matrix > input$expr_threshold) >= (length(input$gene_expr) - input$relaxation)
          } else if (input$expr_threshold > 0) {
            unique_values <- c("other", "cells above threshold")
            color_values <- c("lightgray", "red")
            used_matrix <- used_matrix > input$expr_threshold
          }

          color_plot2(
            embedding = pkg_env$stab_obj$umap,
            color_info = used_matrix,
            plt_height = input$height_gene * ppi, #- gene_legend_height(),
            plt_width = input$width_gene * ppi,
            predicted_height = (gene_legend_height() - 1) / ppi,
            # color_values = function(n) { paletteer::paletteer_c("grDevices::OrRd", n)},
            unique_values = unique_values,
            color_values = color_values,
            pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
            pt_size = input$gene_pt_size,
            text_size = input$gene_text_size,
            legend_text_size = input$gene_legend_size,
            axis_size = input$gene_axis_size,
            display_legend = TRUE
          )
          dev.off()
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
        clustering_1 <- as.matrix(pkg_env$stab_obj$mbs[[as.character(input$jsi_k_1)]])
        df_1 <- data.frame(clustering_1)
        df_1$cell <- rownames(df_1)
        clustering_2 <- as.matrix(pkg_env$stab_obj$mbs[[as.character(input$jsi_k_2)]])
        df_2 <- data.frame(clustering_2)
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
              mat[as.character(n), as.character(m)] <- bulkAnalyseR::jaccard_index(cluster_1, cluster_2)
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

        ggplot2::ggplot(df_mat, ggplot2::aes(Var1, Var2)) +
          ggplot2::geom_tile(ggplot2::aes(fill = value)) +
          ggplot2::geom_text(ggplot2::aes(label = round(value, 2))) +
          ggplot2::scale_fill_gradient2(
            low = scales::muted("darkred"),
            mid = "white",
            high = scales::muted("midnightblue"),
            midpoint = 0
          ) +
          ggplot2::scale_x_continuous(breaks = pretty(df_mat$Var1, n = length(all_clusters_2))) +
          ggplot2::scale_y_continuous(breaks = pretty(df_mat$Var2, n = length(all_clusters_1))) +
          ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white"),
            axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, size = 10, face = "bold"),
            axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
            axis.title = ggplot2::element_text(size = 14, face = "bold"),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))
          ) +
          ggplot2::xlab("Clusters in Configuration 2") +
          ggplot2::ylab("Clusters in Configuration 1") +
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

#' Writing objects
#'
#' @description to be completed
#'
#' @export
server_comparisons <- function(id, chosen_config, chosen_method) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      isolated_chosen_config <- shiny::isolate(chosen_config())
      ftype <- isolated_chosen_config$chosen_feature_type
      fsize <- isolated_chosen_config$chosen_set_size

      isolated_chosen_method <- shiny::isolate(chosen_method())
      cl_method <- isolated_chosen_method$method_name
      k_values <- isolated_chosen_method$n_clusters

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
          inputId = glue::glue("gene_panel_{panels}-gene_expr"),
          choices = c(names(pkg_env$genes_of_interest), names(pkg_env$genes_others)),
          # selected = NULL,
          selected = names(pkg_env$genes_of_interest[1]),
          server = TRUE,
          options = list(
            maxOptions = 7,
            create = TRUE,
            persist = TRUE
          )
        )
      }
      server_comparison_metadata_panel("metadata_panel_left")
      server_comparison_metadata_panel("metadata_panel_right")
      server_comparison_gene_panel("gene_panel_left")
      server_comparison_gene_panel("gene_panel_right")
      server_comparison_jsi("jsi_plot", k_values)
      server_comparison_markers("markers", k_values)

      gc()
    }
  )
}



#       temp_list <- rhdf5::h5read("stability.h5",'/')
#       temp_list$feature_stability <- NULL
#       obj_fsets <- names(temp_list$feature_ordering$stable)

#       slim_obj <- function(x){
#         temp_list[[fset]][[x]]$nn_con_comps <- NULL
#         temp_list[[fset]][[x]]$pca <- NULL
#         temp_list[[fset]][[x]]$clustering_stability$split_by_resolution <- NULL
#         temp_list[[fset]][[x]]$feature_list <- NULL
#       }

#       for (fset in obj_fsets){
#         lapply(temp_list$feature_ordering$stable[[1]], slim_obj)
#       }

#       add_env_variable("stab_obj_1", temp_list)
#       rm(temp_list)
#       gc()
#       # Render UI
#       ns <- shiny::NS(id)
#       options <- names(pkg_env$stab_obj_1$feature_ordering$stable)
#       output$compare_sel_fset_render <- shiny::renderUI({
#         ns <- session$ns
#         shiny::selectInput(
#           inputId = ns("compare_sel_fset"),
#           label = "Select feature - set:",
#           choices = options,
#           selected = options[1],
#           multiple = FALSE
#         )
#       })
#       outputOptions(output, "compare_sel_fset_render", suspendWhenHidden=FALSE)

#       options_2 <- pkg_env$stab_obj_1$feature_ordering$stable[[1]]
#       output$compare_sel_steps_render <- shiny::renderUI({
#         ns <- session$ns
#         shiny::selectInput(
#           inputId = ns("compare_sel_steps"),
#           label = "Select feature - size:",
#           choices = options_2,
#           selected = options_2[1],
#           multiple = FALSE
#         )
#       })
#       outputOptions(output, "compare_sel_steps_render", suspendWhenHidden=FALSE)

#       options_3 <- names(pkg_env$stab_obj_1[[1]][[1]]$clustering_stability$split_by_k$mbs)
#       output$clustering_method_choice_render <- shiny::renderUI({
#         ns <- session$ns
#         shiny::radioButtons(
#           inputId = ns("clustering_method_choice"),
#           label = "Select clustering Method:",
#           choices = options_3,
#           selected = options_3[1],
#           width='50%'
#         )
#       })
#       outputOptions(output, "clustering_method_choice_render", suspendWhenHidden=FALSE)


#       #Set the dropdown menu depending on the available clusters

#       chosen_config <- pkg_env$stab_obj_1[[unlist(chosen_config[1])]][[unlist(chosen_config[2])]]
#       embedding_fixed <- chosen_config$umap
#       chosen_method <- chosen_config$clustering_stability$split_by_k$mbs[[chosen_method]]
#       output$k_selection_fixed <- shiny::renderUI({
#         shiny::selectInput(ns("k_fixed"), "Select a number of clusters:", choices = names(chosen_method), selected = names(chosen_method)[1], multiple = FALSE)
#       })
#       output$render_meta_1 <- shiny::renderUI({
#         shiny::selectInput(ns("col_by_fixed"), "Colour by:", choices = c('ECC','Clusters',colnames(metadata$metadata)), selected='ECC', multiple = FALSE)
#       })
#       output$render_meta_2 <- shiny::renderUI({
#         shiny::selectInput(ns("col_by_fixed_2"), "Colour by:", choices = c('ECC','Clusters',colnames(metadata$metadata)), selected='Clusters', multiple = FALSE)
#       })
#       output$render_meta_3 <- shiny::renderUI({
#         shiny::selectInput(ns("col_by_choice"), "Colour by:", choices = c('ECC','Clusters',colnames(metadata$metadata)), selected='ECC', multiple = FALSE)
#       })
#       output$render_meta_4 <- shiny::renderUI({
#         shiny::selectInput(ns("col_by_choice_2"), "Colour by:", choices = c('ECC','Clusters',colnames(metadata$metadata)), selected='Clusters', multiple = FALSE)
#       })
#       output$k_selection <- shiny::renderUI({
#         shiny::selectInput(ns("k"), "Select a number of clusters:", choices = names(pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k$mbs[[input$clustering_method_choice]]), multiple = FALSE)
#       })

#       #Set all UMAPs: 1,3 belong to the fixed selection, 2,4 are the UMAPs that change depending on the users choice
#       stable_graph_fixed <- paste0(toupper(chosen_config$stable_config$base_embedding),'_',chosen_config$stable_config$graph_type)
#       stable_nn_fixed <- as.character(chosen_config$stable_config$n_neighbours)
#       stable_ecc_fixed <- chosen_config$nn_stability$n_neigh_ec_consistency[[stable_graph_fixed]][[stable_nn_fixed]]

#       umap_1 <- shiny::reactive({
#         if (input$col_by_fixed=='ECC'){
#           ECC <- stable_ecc_fixed
#           ggplot2::ggplot(data.frame(
#             chosen_config$umap),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = ECC)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::scale_color_viridis_c() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }else if(input$col_by_fixed=='Clusters'){
#           Clusters <- as.factor(as.matrix(chosen_method[[input$k_fixed]]))
#           ggplot2::ggplot(data.frame(
#             embedding_fixed),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = Clusters)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }else{
#           Feat <- metadata$metadata[[input$col_by_fixed]]
#           ggplot2::ggplot(data.frame(
#             embedding_fixed),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = Feat)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }
#       })
#       umap_3 <- shiny::reactive({
#         if (input$col_by_fixed_2=='ECC'){
#           ECC <- stable_ecc_fixed
#           ggplot2::ggplot(data.frame(
#             embedding_fixed),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = ECC)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::scale_color_viridis_c() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }else if(input$col_by_fixed_2=='Clusters'){
#           if(is.null(input$k_fixed)) {
#             return(ggplot2::ggplot() + ggplot2::theme_void())
#           }else{
#             Clusters <- as.factor(as.matrix(chosen_method[[input$k_fixed]]))
#             ggplot2::ggplot(data.frame(
#               embedding_fixed),
#               ggplot2::aes(x = .data$X1,
#                            y = .data$X2,
#                            color = Clusters)) +
#               ggplot2::xlab('UMAP 1') +
#               ggplot2::ylab('UMAP 2') +
#               ggplot2::geom_point() +
#               ggplot2::theme_bw() +
#               ggplot2::coord_fixed()
#           }
#         }else{
#           Feat <- metadata$metadata[[input$col_by_fixed_2]]
#           ggplot2::ggplot(data.frame(
#             embedding_fixed),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = Feat)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }
#       })
#       umap_2 <- shiny::reactive({
#         if(is.null(input$k)) {
#           return(ggplot2::ggplot() + ggplot2::theme_void())
#         }
#         if (input$col_by_choice=='ECC'){
#           ECC <- pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$nn_stability$n_neigh_ec_consistency[[paste0(toupper(pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$stable_config$base_embedding),'_',pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$stable_config$graph_type)]][[as.character(pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$stable_config$n_neighbours)]]
#           ggplot2::ggplot(data.frame(
#             pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = ECC)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::scale_color_viridis_c() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }else if(input$col_by_choice=='Clusters'){
#           Clusters <- as.factor(as.matrix(pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k$mbs[[input$clustering_method_choice]][[input$k]]))
#           ggplot2::ggplot(data.frame(
#             pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = Clusters)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }else{
#           Feat <- metadata$metadata[[input$col_by_choice]]
#           ggplot2::ggplot(data.frame(
#             pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = Feat)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }
#       })
#       umap_4 <- shiny::reactive({
#         if(is.null(input$k)) {
#           return(ggplot2::ggplot() + ggplot2::theme_void())
#         }
#         if (input$col_by_choice_2=='ECC'){
#           ECC <- stable_ecc_choice()
#           ggplot2::ggplot(data.frame(
#             pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = ECC)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::scale_color_viridis_c() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }else if(input$col_by_choice_2=='Clusters'){
#           if(is.null(input$k)) {
#             return(ggplot2::ggplot() + ggplot2::theme_void())
#           }else{
#             Clusters <- as.factor(as.matrix(pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k$mbs[[input$clustering_method_choice]][[input$k]]))
#             ggplot2::ggplot(data.frame(
#               pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
#               ggplot2::aes(x = .data$X1,
#                            y = .data$X2,
#                            color = Clusters)) +
#               ggplot2::xlab('UMAP 1') +
#               ggplot2::ylab('UMAP 2') +
#               ggplot2::geom_point() +
#               ggplot2::theme_bw() +
#               ggplot2::coord_fixed()
#           }
#         }else{
#           Feat <- metadata$metadata[[input$col_by_choice_2]]
#           ggplot2::ggplot(data.frame(
#             pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
#             ggplot2::aes(x = .data$X1,
#                          y = .data$X2,
#                          color = Feat)) +
#             ggplot2::xlab('UMAP 1') +
#             ggplot2::ylab('UMAP 2') +
#             ggplot2::geom_point() +
#             ggplot2::theme_bw() +
#             ggplot2::coord_fixed()
#         }
#       })
#       #Output all of the UMAPs
#       output$umap_fixed <- shiny::renderPlot({
#         if(is.null(umap_1())) {
#           return(ggplot2::ggplot() + ggplot2::theme_void())
#         }
#         umap_1()
#       })
#       output$umap_choice <- shiny::renderPlot({
#         if(is.null(umap_2())) {
#           return(ggplot2::ggplot() + ggplot2::theme_void())
#         }
#         umap_2()
#       })
#       output$umap_fixed_2 <- shiny::renderPlot({
#         if(is.null(umap_3())) {
#           return(ggplot2::ggplot() + ggplot2::theme_void())
#         }
#         umap_3()
#       })
#       output$umap_choice_2 <- shiny::renderPlot({
#         if(is.null(umap_4())) {
#           return(ggplot2::ggplot() + ggplot2::theme_void())
#         }
#         umap_4()
#       })

#       umap_filetype_1 <- shiny::reactive({
#         if (input$umap_filetype_1=='PDF'){
#           filename <- paste0(input$filename_umap_1,'.pdf')
#           return(filename)
#         }else if (input$umap_filetype_1=='PNG'){
#           filename <- paste0(input$filename_umap_1,'.png')
#           return(filename)
#         }else{
#           filename <- paste0(input$filename_umap_1,'.svg')
#           return(filename)
#         }
#       })
#       output$download_umap_1 <- shiny::downloadHandler(
#         filename = function() {umap_filetype_1()},
#         content = function(file) {
#           ggplot2::ggsave(file, umap_1(),width = input$width_umap_1,
#                           height = input$height_umap_1,
#                           units = "in",
#                           limitsize = FALSE)
#         }
#       )
#       umap_filetype_2 <- shiny::reactive({
#         if (input$umap_filetype_2=='PDF'){
#           filename <- paste0(input$filename_umap_2,'.pdf')
#           return(filename)
#         }else if (input$umap_filetype_2=='PNG'){
#           filename <- paste0(input$filename_umap_2,'.png')
#           return(filename)
#         }else{
#           filename <- paste0(input$filename_umap_2,'.svg')
#           return(filename)
#         }
#       })
#       output$download_umap_2 <- shiny::downloadHandler(
#         filename = function() {umap_filetype_2()},
#         content = function(file) {
#           ggplot2::ggsave(file, umap_2(),width = input$width_umap_2,
#                           height = input$height_umap_2,
#                           units = "in",
#                           limitsize = FALSE)
#         }
#       )
#       umap_filetype_3 <- shiny::reactive({
#         if (input$umap_filetype_3=='PDF'){
#           filename <- paste0(input$filename_umap_3,'.pdf')
#           return(filename)
#         }else if (input$umap_filetype_3=='PNG'){
#           filename <- paste0(input$filename_umap_3,'.png')
#           return(filename)
#         }else{
#           filename <- paste0(input$filename_umap_3,'.svg')
#           return(filename)
#         }
#       })
#       output$download_umap_3 <- shiny::downloadHandler(
#         filename = function() {umap_filetype_3()},
#         content = function(file) {
#           ggplot2::ggsave(file, umap_3(),width = input$width_umap_3,
#                           height = input$height_umap_3,
#                           units = "in",
#                           limitsize = FALSE)
#         }
#       )

#       umap_filetype_4 <- shiny::reactive({
#         if (input$umap_filetype_4=='PDF'){
#           filename <- paste0(input$filename_umap_4,'.pdf')
#           return(filename)
#         }else if (input$umap_filetype_4=='PNG'){
#           filename <- paste0(input$filename_umap_4,'.png')
#           return(filename)
#         }else{
#           filename <- paste0(input$filename_umap_4,'.svg')
#           return(filename)
#         }
#       })
#       output$download_umap_4 <- shiny::downloadHandler(
#         filename = function() {umap_filetype_4()},
#         content = function(file) {
#           ggplot2::ggsave(file, umap_4(),width = input$width_umap_4,
#                           height = input$height_umap_4,
#                           units = "in",
#                           limitsize = FALSE)
#         }
#       )

#       #JSI heatmap
#       barcode_heatmap <- shiny::reactive({
#         if(is.null(input$k)) {
#           return(ggplot2::ggplot() + ggplot2::theme_void())
#         }
#         df_fixed <- data.frame(as.matrix(chosen_method[[input$k_fixed]]))
#         df_fixed$cell <- rownames(df_fixed)
#         clustering_choice <- as.matrix(pkg_env$stab_obj_1[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k$mbs[[input$clustering_method_choice]][[input$k]])
#         df_choice <- data.frame(clustering_choice)
#         df_choice$cell <- rownames(df_choice)
#         all_clusters_1 <- unique(df_fixed[,1])
#         all_clusters_2 <- unique(df_choice[,1])

#         mat = matrix(, nrow = length(all_clusters_2),
#                      ncol = length(all_clusters_1))
#         colnames(mat) <- as.factor(sort(all_clusters_1))
#         rownames(mat) <- as.factor(sort(all_clusters_2))
#         if (input$heatmap_type=='JSI'){
#           for (m in all_clusters_1){
#             cluster_1 <- rownames(df_fixed[df_fixed[,1]==m,])
#             for (n in all_clusters_2){
#               cluster_2 <- rownames(df_choice[df_choice[,1]==n,])
#               mat[as.character(n),as.character(m)] <- bulkAnalyseR::jaccard_index(cluster_1, cluster_2)
#               label <- 'JSI'
#             }
#           }
#         }else{
#           for (m in all_clusters_1){
#             cluster_1 <- rownames(df_fixed[df_fixed[,1]==m,])
#             for (n in all_clusters_2){
#               cluster_2 <- rownames(df_choice[df_choice[,1]==n,])
#               mat[as.character(n),as.character(m)] <- length(intersect(cluster_1, cluster_2))
#               label <- 'Shared cells'
#             }
#           }
#         }

#         df_mat <- reshape::melt(mat)

#         ggplot2::ggplot(df_mat, ggplot2::aes(X1, X2)) +
#           ggplot2::geom_tile(ggplot2::aes(fill = value)) +
#           ggplot2::geom_text(ggplot2::aes(fill = value, label = round(value, 2))) +
#           ggplot2::scale_fill_gradient2(low = scales::muted("darkred"),
#                                         mid = "white",
#                                         high = scales::muted("midnightblue"),
#                                         midpoint = 0) +
#           ggplot2::scale_x_continuous(breaks = pretty(df_mat$X1, n = length(all_clusters_2))) +
#           ggplot2::scale_y_continuous(breaks = pretty(df_mat$X2, n = length(all_clusters_1))) +
#           ggplot2::theme(
#             panel.background=ggplot2::element_rect(fill="white"),
#             axis.text.x = ggplot2::element_text(hjust = 1,vjust=1,size = 10,face = "bold"),
#             axis.text.y = ggplot2::element_text(size = 10,face = "bold"),
#             axis.title = ggplot2::element_text(size=14,face="bold"),
#             axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
#             axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))) +
#           ggplot2::xlab("Clusters in Configuration 2") +
#           ggplot2::ylab("Clusters in Configuration 1") +
#           ggplot2::labs(fill=label)
#       })
#       output$barcode_heatmap <- shiny::renderPlot({
#         barcode_heatmap()
#       },height=500)
#       heatmap_filetype <- shiny::reactive({
#         if (input$heatmap_filetype=='PDF'){
#           filename <- paste0(input$filename_heatmap,'.pdf')
#           return(filename)
#         }else if (input$heatmap_filetype=='PNG'){
#           filename <- paste0(input$filename_heatmap,'.png')
#           return(filename)
#         }else{
#           filename <- paste0(input$filename_heatmap,'.svg')
#           return(filename)
#         }
#       })
#       output$download_heatmap <- shiny::downloadHandler(
#         filename = function() {heatmap_filetype()},
#         content = function(file) {
#           ggplot2::ggsave(file,barcode_heatmap(),width = input$width_heatmap,
#                           height = input$height_heatmap,
#                           units = "in",
#                           limitsize = FALSE)
#         }
#       )
#     }
#   )
# }
