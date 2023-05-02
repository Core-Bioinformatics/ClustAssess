proportion_widths <- 55
#### UI ####
gear_overall <- function(ns, id) {
  shiny::tagList(
    shiny::sliderInput(
      inputId = ns(paste0(id, "_boxplot_width")),
      label = "Boxplot width",
      min = 0.00, max = 1.00, value = 0.50
    ),
    shiny::sliderInput(
      inputId = ns(paste0(id, "_inter_distance")),
      label = "Space between groups",
      min = 1, max = 15, value = 1, step = 1
    ),
    shiny::sliderInput(
      inputId = ns(paste0(id, "_intra_distance")),
      label = "Space inside groups",
      min = 1, max = 15, value = 1, step = 1
    ),
    shiny::sliderInput(
      inputId = ns(paste0(id, "_text_size")),
      label = "Text size",
      min = 0.10, max = 10.00, value = 1.00, step = 0.1
    ),
    shiny::sliderInput(
      inputId = ns(paste0(id, "_text_offset")),
      label = "Boxplot text offset",
      min = 0.005, max = 0.15, value = 0.01, step = 0.005
    ),
  )
}

gear_umaps <- function(ns, id) {
  shinyWidgets::dropdownButton(
    shiny::tagList(
      shiny::sliderInput(
        inputId = ns(paste0(id, "_text_size")),
        label = "Text size",
        min = 0.10, max = 10.00, value = 1.00, step = 0.1
      ),
      shiny::sliderInput(
        inputId = ns(paste0(id, "_axis_size")),
        label = "Axis labels size",
        min = 0.10, max = 10.00, value = 1.00, step = 0.1
      ),
      shiny::sliderInput(
        inputId = ns(paste0(id, "_pt_size")),
        label = "Point size",
        min = 0.05, max = 5.00, value = 0.10, step = 0.05
      ),
      shinyWidgets::radioGroupButtons(
        inputId = ns(paste0(id, "_pt_type")),
        label = "Point type",
        choices = c("Pixel", "Circle")
      ),
      shinyWidgets::prettySwitch(
        inputId = ns(paste0(id, "_labels")),
        label = "Show labels",
        status = "success",
        fill = TRUE
      )
    ),
    circle = TRUE,
    status = "info",
    size = "sm",
    icon = shiny::icon("gear")
  )
}


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
      shiny::h2("ECC per individual resolution values")
    ),
    shinyWidgets::dropdownButton(
      shinyWidgets::sliderTextInput(
        inputId = ns("by_step_resolution"),
        label = "Resolution",
        choices = c("")
      ),
      gear_overall(ns, "by_step_res"),
      circle = TRUE,
      status = "info",
      size = "sm",
      icon = shiny::icon("gear")
    ),
    shiny::splitLayout(
      cellWidths = paste0(c(proportion_widths - 3, 100 - proportion_widths), "%"), #c("52%", "45%"),
      shiny::tagList(
        shiny::plotOutput(ns("boxplot_ecc"), height = "auto"),
        shiny::splitLayout(
          cellWidths = c("40px", "90%"),
          shinyWidgets::circleButton(ns("info_ecs_incr_res"),
            icon = shiny::icon("info"),
            size = "sm",
            status = "success"
          ),
          # shiny::h2("Incremental ECS per individual resolution values")
        ),
        shinyWidgets::dropdownButton(
          shinyWidgets::sliderTextInput(
            inputId = ns("incremental_resolution"),
            label = "Resolution",
            choices = c("")
          ),
          gear_overall(ns, "incremental_res"),
          circle = TRUE,
          status = "info",
          size = "sm",
          icon = shiny::icon("gear")
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
        shiny::plotOutput(ns("umap_ecc"), height = "auto", width = "100%"),
        shiny::plotOutput(ns("umap_ecc_legend"), height = "auto", width = "100%"),
        shiny::tableOutput(ns("table_ecc_info"))


        # shiny::textOutput(ns("boxplot_ecc_info")),
        # shiny::plotOutput(ns("umap_ecc"), height = "350px", width = "90%")
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
    shinyWidgets::dropdownButton(
      gear_overall(ns, "by_step"),
      circle = TRUE,
      status = "info",
      size = "sm",
      icon = shiny::icon("gear")
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
    shinyWidgets::dropdownButton(
      gear_overall(ns, "incremental"),
      circle = TRUE,
      status = "info",
      size = "sm",
      icon = shiny::icon("gear")
    ),
    shiny::plotOutput(ns("overall_boxplot_incremental"), height = "auto")
  )
}

ui_dimensionality_distribution_plots <- function(id, draw_line) {
  ns <- shiny::NS(id)
  style <- ifelse(draw_line, "border-right:5px solid", "")

  # shiny::wellPanel(
    shinyWidgets::panel(
    style = style,
    shiny::selectizeInput(ns("feature_type"), "Feature names", NULL),
    shiny::selectizeInput(ns("feature_steps"), "Feature set size", NULL),
    shiny::hr(),
    shiny::verticalLayout(
      shiny::selectizeInput(
        inputId = ns("gene_expr"),
        choices = NULL,
        label = "Gene name(s)",
        width = "95%",
        multiple = TRUE
      ),
      shiny::sliderInput(
        inputId = ns("expr_threshold"),
        label = "Gene expression threshold",
        min = 0, max = 10, value = 0,
        width = "95%"
      )
    ),
    gear_umaps(ns, "gene"),
    # shiny::actionButton(ns("download_gene_action"), "Download", icon = shiny::icon("download")),
    shiny::downloadButton(ns("download_gene"), "Download"),
    shiny::splitLayout(
      shiny::numericInput(ns("gene_height"), "Plot height (in)", value = 7),
      shiny::numericInput(ns("gene_width"), "Plot width (in)", value = 7)
    ),
    # shiny::uiOutput(ns("umap_gene_generator")),
    shiny::plotOutput(ns("umap_gene"), height = "auto"),
    shiny::plotOutput(ns("umap_gene_legend"), height = "auto"),
    shiny::selectizeInput(
      inputId = ns("metadata"),
      label = "Metadata",
      choices = NULL
    ),
    # shiny::splitLayout(
      # cellWidths = c("40px", "90%"),
    # shiny::splitLayout(
      # cellWidths = "40%",
      shiny::fluidRow(
        shiny::column(6,

        gear_umaps(ns, "metadata"),
          # shiny::actionButton(ns("download_metadata_action"), "Download", icon = shiny::icon("download")),
        shiny::downloadButton(ns("download_metadata"), "Download"),
        shiny::splitLayout(
          shiny::numericInput(ns("metadata_height"), "Plot height (in)", value = 7),
          shiny::numericInput(ns("metadata_width"), "Plot width (in)", value = 7)
        ),
        ),
        shiny::column(6,
        shinyWidgets::pickerInput(
                inputId = ns("select_groups"),
                choices = "",
                inline = FALSE,
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
        ),
      ),
    # ),
    # shiny::uiOutput(ns("umap_metadata_generator"))
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

#' Dimensionality reduction - ui side
#'
#' @description to be completed
#'
#' @export
ui_dimensionality_reduction <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Dimensionality Reduction",
    shiny::splitLayout(
      cellWidths = c("40px", "90%"),
      shinyWidgets::circleButton(ns("info_title"),
        icon = shiny::icon("info"),
        size = "sm",
        status = "success"
      ),
      shiny::h2("Assessing the stability of the dimensionality reduction"),
    ),
    shiny::uiOutput(ns("test")),
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
      # TODO add the ecc umap on the right with select
      # TODO less urgent add option to zoom in on the ecc umap
      # boxplot_ecc_df <- shiny::reactive(
      #   plot_feature_per_resolution_stability_boxplot(
      #     feature_object_list = pkg_env$stab_obj,
      #     resolution = input$resolution,
      #     return_df = TRUE
      #   )
      # )

      # # ecc_mouse_hover <- shiny::reactiveVal(c(NA, NA))
      # ecc_mouse_hover <- shiny::reactive({
      #   if (is.null(input$ecc_hover)) {
      #     return(c(NA, NA))
      #   }
      #   x <- input$ecc_hover$x
      #   # y <- input$ecc_hover$y
      #   min_step <- min(as.numeric(boxplot_ecc_df()$step_index))
      #   max_step <- max(as.numeric(boxplot_ecc_df()$step_index))
      #   nconfigs <- length(pkg_env$feature_ordering$original)
      #   chosen_step <- NULL
      #   chosen_config <- NULL

      #   for (steps_intervals in seq(from = min_step - 0.5, to = max_step - 0.5, by = 1)) {
      #     if (x < steps_intervals + 1) {
      #       chosen_step <- steps_intervals + 0.5
      #       break
      #     }
      #   }

      #   if (is.null(chosen_step)) {
      #     chosen_step <- max_step
      #   }

      #   for (i in seq_len(nconfigs)) {
      #     if (x < (steps_intervals + i / nconfigs)) {
      #       chosen_config <- i
      #       break
      #     }
      #   }

      #   if (is.null(chosen_config)) {
      #     chosen_config <- nconfigs
      #   }

      #   return(c(chosen_step, chosen_config))
      # }) %>% shiny::bindEvent(input$ecc_hover)


      # output$boxplot_ecc_info <- shiny::renderText({
      #   if (is.na(ecc_mouse_hover()[1])) {
      #     return("")
      #   }

      #   feature_type <- pkg_env$feature_types[ecc_mouse_hover()[2]]
      #   nfeatures <- pkg_env$feature_ordering$original[[feature_type]][ecc_mouse_hover()[1]]
      #   glue::glue("The ECC summary for {feature_type} - {nfeatures} and resolution {input$resolution}.")
      #   # Click if you want to visualise the UMAP distribution.")
      # })

      # output$table_ecc_info <- shiny::renderTable({
      #   if (is.na(ecc_mouse_hover()[1])) {
      #     return()
      #   }
      #   summary_ecc <- fivenum(pkg_env$stab_obj$by_steps[[ecc_mouse_hover()[2]]][[ecc_mouse_hover()[1]]][[input$resolution]]$ecc)
      #   data.frame(
      #     min = summary_ecc[1],
      #     Q1 = summary_ecc[2],
      #     median = summary_ecc[3],
      #     Q3 = summary_ecc[4],
      #     max = summary_ecc[5]
      #   )
      # })

      # shiny::observe({
      #   output$boxplot_ecc <- shiny::renderPlot({
      #     plot_feature_per_resolution_stability_boxplot(
      #       feature_object_list = pkg_env$stab_obj,
      #       resolution = input$resolution,
      #       text_size = input$text_size,
      #       boxplot_width = input$boxplot_width,
      #       dodge_width = input$dodge_width
      #     ) + ggplot2::theme(legend.position = "bottom")
      #   })

      #   output$boxplot_incr <- shiny::renderPlot({
      #     plot_feature_per_resolution_stability_incremental(
      #       feature_object_list = pkg_env$stab_obj,
      #       resolution = input$resolution,
      #       text_size = input$text_size,
      #       boxplot_width = input$boxplot_width,
      #       dodge_width = input$dodge_width
      #     )
      #   })
      # }) %>% shiny::bindEvent(input$resolution)


      # shiny::observe({
      #   hover_coords <- ecc_mouse_hover()

      #   output$umap_ecc <- shiny::renderPlot({
      #     # if (is.na(hover_coords[1])) {
      #     #   return(empty_ggplot())
      #     # }

      #     config_name <- pkg_env$feature_types[hover_coords[2]]
      #     step_index <- pkg_env$feature_ordering$original[[config_name]][hover_coords[1]]
      #     ecc <- pkg_env$stab_obj$by_steps[[config_name]][[step_index]][[input$resolution]]$ecc

      #     color_ggplot(pkg_env$stab_obj$embedding_list[[config_name]][[step_index]], ecc) +
      #       # ggplot2::guides(color = ggplot2::guide_legend(title = "ecc")) +
      #       ggplot2::scale_color_viridis_c(name = "ECC") +
      #       ggplot2::ggtitle(sprintf("ECC distribution\n%s - %s\nResolution %s", config_name, step_index, input$resolution))
      #   })
      # }) %>% shiny::bindEvent(input$ecc_click)

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

      plt_height <- shiny::reactive(
         floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * (1 - proportion_widths / 100)))
      )

      legend_height <- shiny::reactive({
        # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
        pdf(file = NULL, width = plt_height(), height = plt_height())
        par(mai = c(0.1, 0, 0.1, 0))
        text_height <- strheight("TE\nXT\n", units = "inches", cex = input$by_step_res_text_size)
        dev.off()
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
           plt_height()
        },
        width = function() {
           plt_height()
        },
        {

          shiny::req(
            ecc_value()
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
            plt_width  = plt_height(),
            text_size = input$by_step_res_text_size
          )
      })

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
            text_size = input$by_step_res_text_size,
            color_values = NULL,
            unique_values = NULL,
            color_info = ecc_value()
          )
        }
      )

      output$table_ecc_info <- shiny::renderTable({
        shiny::req(ecc_value())

        data.frame(
          "Stat" = c("Min", "Q1", "Median", "Q2", "Max"),
          "ECC value" = stats::fivenum(ecc_value())
        )
      })

      output$boxplot_ecc <- shiny::renderPlot(
        {
          shiny::req(input$by_step_resolution != "", input$by_step_resolution != "0")
          print(input$by_step_resolution)

          shiny_plot_feature_stability_boxplot(
            resval = input$by_step_resolution,
            feature_ordering = pkg_env$feature_ordering,
            # plt_height = floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2)),
            plt_height = floor(pkg_env$height_ratio * pkg_env$dimension()[2]),
            plt_width = pkg_env$dimension()[1] * (proportion_widths - 5) / 100,
            text_size = input$by_step_res_text_size,
            width = input$by_step_res_boxplot_width,
            space_intra_groups = input$by_step_res_intra_distance,
            plt_title = paste("res", input$by_step_resolution),
            space_inter_groups = input$by_step_res_inter_distance,
            text_offset = input$by_step_res_text_offset
          )
        },
        height = function() {
          floor(pkg_env$height_ratio * pkg_env$dimension()[2])
        }
      )

      output$boxplot_incr <- shiny::renderPlot({
        shiny::req(input$by_step_resolution != "", input$by_step_resolution != "0")

        shiny_plot_feature_stability_incremental(
          resval = input$incremental_resolution,
          feature_ordering = pkg_env$feature_ordering,
          # plt_height = floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2)),
          plt_height = floor(pkg_env$height_ratio * pkg_env$dimension()[2]),
          plt_width = pkg_env$dimension()[1] * (proportion_widths - 5) / 100,
          text_size = input$incremental_res_text_size,
          plt_title = paste("res", input$incremental_resolution),
          width = input$incremental_res_boxplot_width,
          space_intra_groups = input$incremental_res_intra_distance,
          space_inter_groups = input$incremental_res_inter_distance,
          text_offset = input$incremental_res_text_offset
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
            space_inter_groups = input$by_step_inter_distance,
            text_offset = input$by_step_text_offset
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
            space_inter_groups = input$incremental_inter_distance,
            text_offset = input$incremental_text_offset
          )
        },
        height = function() {
          floor(pkg_env$height_ratio * pkg_env$dimension()[2])
        }
      )

      # output$overall_boxplot_incremental <- shiny::renderPlot({
      #   plot_feature_overall_stability_incremental(
      #     feature_object_list = pkg_env$stab_obj,
      #     text_size = input$text_size,
      #     boxplot_width = input$boxplot_width,
      #     dodge_width = input$dodge_width
      #   )
      # })

      shiny::observe(dr_individual_ecc_info(session)) %>% shiny::bindEvent(input$info_ecc_res)
      shiny::observe(dr_individual_incremental_info(session)) %>% shiny::bindEvent(input$info_ecs_incr_res)
      shiny::observe(dr_overall_ecc_info(session)) %>% shiny::bindEvent(input$info_ecc_overall)
      shiny::observe(dr_overall_incremental_info(session)) %>% shiny::bindEvent(input$info_ecs_incremental_overall)
    }
  )
}

server_dimensionality_distribution <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      shiny::observeEvent(input$feature_type, {
        shiny::updateSelectizeInput(session,
          inputId = "feature_steps",
          choices = pkg_env$feature_ordering$stable[[input$feature_type]],
          server = FALSE
        )
      })

      # TODO check the performance
      # TODO metadata highlight only a group of cells
      # TODO adapt the download

      # shiny::observe({
      #   height <- floor(input$dimension * pkg_env$height_ratio)
      #   ns <- session$ns
      #   # output$umap_gene_generator <- shiny::renderUI(
      #   #   shiny::plotOutput(ns("umap_gene"), height = paste0(height, "px"))
      #   # )

      #   # output$umap_metadata_generator <- shiny::renderUI(
      #   #   shiny::plotOutput(ns("umap_metadata"), height = paste0(height, "px"))
      #   # )
      # }) %>% shiny::bindEvent(input$dimension)

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
        shiny::updateSliderInput(session,
          inputId = "expr_threshold",
          max = round(max_level_expr(), 3),
          step = round(max_level_expr() / 10, 3)
        )
      }) %>% shiny::bindEvent(input$gene_expr)

      shiny::observe({
        mtd_names <- pkg_env$metadata_unique[[input$metadata]]
        if (is.null(mtd_names)) {
          shinyjs::hide(id = "select_groups")

        } else {
          shinyjs::show(id = "select_groups")
          shinyWidgets::updatePickerInput(
            session,
            inputId = "select_groups",
            choices = mtd_names,
            selected = mtd_names
          )
        }

      }) %>% shiny::bindEvent(input$metadata)


     

      # umap_gene_plot <- shiny::reactive({
      #   if (is.null(input$feature_type) || is.null(input$feature_steps)) {
      #     Sys.sleep(0.5) # wait a little
      #   }

      #   embedding <- pkg_env$stab_obj$embedding_list[[input$feature_type]][[input$feature_steps]]

      #   if (length(input$gene_expr) == 0) {
      #     return(empty_ggplot())
      #   }

      #   # if (is.null(embedding)) {
      #   #   return(empty_ggplot())
      #   # }

      #   if (length(input$gene_expr) == 1) {
      #     return(expression_ggplot(embedding, expr_matrix(), input$expr_threshold) + ggplot2::theme(aspect.ratio = 1))
      #   }

      #   voting_scheme_ggplot(
      #     embedding,
      #     matrixStats::colSums2(expr_matrix() >= input$expr_threshold),
      #     nrow(expr_matrix())
      #   ) #+ ggplot2::theme(aspect.ratio = 1)
      # })

      # umap_metadata_plot <- shiny::reactive({
      #   if (is.null(input$feature_type) || is.null(input$feature_steps)) {
      #     Sys.sleep(0.5) # wait a little
      #   }

      #   embedding <- pkg_env$stab_obj$embedding_list[[input$feature_type]][[input$feature_steps]]

      #   # if (is.null(embedding)) {
      #   #   return(empty_ggplot())
      #   # }

      #   metadata_ggplot(embedding, input$metadata) #+ ggplot2::theme(aspect.ratio = 1)
      # })

      # output$umap_gene <- shiny::renderCachedPlot(
      #   {
      #     umap_gene_plot()
      #   },
      #   cacheKeyExpr = {
      #     list(input$feature_type, input$feature_steps, input$gene_expr, input$expr_threshold)
      #   },
      #   sizePolicy = shiny::sizeGrowthRatio(
      #     height = floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2)),
      #     width = floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2)),
      #     growthRate = 1.2
      #   )
      # )

      # output$umap_metadata <- shiny::renderCachedPlot(
      #   {
      #     umap_metadata_plot()
      #   },
      #   cacheKeyExpr = {
      #     list(input$feature_type, input$feature_steps, input$metadata)
      #   },
      #   sizePolicy = shiny::sizeGrowthRatio(
      #     height = 500, width = 500, growthRate = 1.2
      #   )
      # )

      plt_height <- shiny::reactive(
        floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.43))
      )

      gene_legend_height <- shiny::reactive({
        # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
        pdf(NULL, width =  plt_height(), height = plt_height())
        par(mai = c(0.1, 0, 0.1, 0))
        text_height <- strheight("TE\nXT\n", units = "inches", cex = input$gene_text_size)
        dev.off()
        return((0.2 + text_height) * ppi)
      })

      metadata_legend_height <- shiny::reactive({
        unique_values <- pkg_env$metadata_unique[[input$metadata]]
        # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
        pdf(file = NULL, width = plt_height(), height = plt_height())
        if (is.null(unique_values)) {
          par(mai = c(0.1, 0, 0.1, 0))
          text_height <- strheight("TE\nXT\n", units = "inches", cex = input$metadata_text_size)
          dev.off()
          return((0.2 + text_height) * ppi)
        }

        par(mar = c(0, 0, 0, 0))
        predicted_width <- strwidth(c(" ", unique_values), units = "inches", cex = input$metadata_text_size) * ppi
        space_width <- predicted_width[1]
        predicted_width <- predicted_width[2:length(predicted_width)]

        number_columns <- min(
          max(
            plt_height() %/% (5 * space_width + max(predicted_width)),
            1),
          length(unique_values)
        )
        number_rows <- ceiling(length(unique_values) / number_columns)

        text_height <- strheight(paste(
          rep("TEXT", number_rows + 1),
          collapse = "\n"
          ),
          units = "inches",
          cex = input$metadata_text_size) * ppi

        dev.off()

        return(text_height)
      })

      output$umap_gene <- shiny::renderPlot(
        height = function() {
          plt_height()
        },
        width = function() {
          plt_height()
        },
        {
          shiny::req(input$feature_steps, input$expr_threshold, input$gene_expr, expr_matrix())
          
          unique_values <- NULL
          used_matrix <- expr_matrix()
          color_values <- function(n) { grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n) }
          if (length(input$gene_expr) > 1) {
            unique_values <- c("FALSE", "TRUE")
            color_values <- c("lightgray", "red")
            used_matrix <- matrixStats::colSums2(used_matrix >= input$expr_threshold) == length(input$gene_expr)
          } else if (input$expr_threshold > 0) {
            unique_values <- c("FALSE", "TRUE")
            color_values <- c("lightgray", "red")
            used_matrix <- used_matrix >= input$expr_threshold
          }

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
            # color_values = function(n) { paletteer::paletteer_c("grDevices::OrRd", n)},
            unique_values = unique_values,
            color_values = color_values,
            pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
            pt_size = input$gene_pt_size,
            text_size = input$gene_text_size
          )
        }
      )

      output$umap_gene_legend <- shiny::renderPlot(
        height = function() {
          gene_legend_height()
        },
        width = function() {
          plt_height()
        },
        {
          shiny::req(gene_legend_height(), input$feature_steps, expr_matrix())

          unique_values <- NULL
          color_values <- function(n) {
            pdf(file = NULL)
            colour_init <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n)
            dev.off()

            return(colour_init)
          }

          if (length(input$gene_expr) > 1 || input$expr_threshold > 0) {
            unique_values <- c("other", "cells above threshold")
            color_values <- c("lightgray", "red")
          }

          # pdf(file = NULL)
          only_legend_plot(
            plt_width = plt_height(),
            text_size = input$gene_text_size,
            color_values = color_values,
            unique_values = unique_values,
            color_info = expr_matrix()
          )
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
          shiny::req(input$feature_steps, input$metadata)
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
            text_size = input$metadata_text_size,
            axis_size = input$metadata_axis_size,
            labels = input$metadata_labels,
            groups_highlight = input$select_groups
          )
        }
      )

        output$umap_metadata_legend <- shiny::renderPlot(
          height = function() {
            metadata_legend_height()
          },
          width = function() {
            plt_height()
          },
          {
            shiny::req(metadata_legend_height(), input$metadata)
            # pdf(file = NULL)
            # dev.off()
            only_legend_metadata_plot(
              metadata_name = input$metadata,
              plt_width = plt_height(),
              groups = input$select_groups,
              text_size = input$metadata_text_size 
            )
          }
        )

        output$download_gene <-  shiny::downloadHandler(
          filename = function() {
            "feature_genes.pdf"
          },
          content = function(file) {
            # ggplot2::ggsave(file, to_save_plot, width = width, height = height)
            shiny::req(input$feature_steps, input$expr_threshold, input$gene_expr, input$gene_width, input$gene_height, expr_matrix())
            pdf(file, width = input$metadata_width, height = input$metadata_height)
            unique_values <- NULL
            used_matrix <- expr_matrix()
            color_values <- function(n) { grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n) }
            if (length(input$gene_expr) > 1) {
              unique_values <- c("other", "cells above threshold")
              color_values <- c("lightgray", "red")
              used_matrix <- matrixStats::colSums2(used_matrix >= input$expr_threshold) == length(input$gene_expr)
            } else if (input$expr_threshold > 0) {
              unique_values <- c("other", "cells above threshold")
              color_values <- c("lightgray", "red")
              used_matrix <- used_matrix >= input$expr_threshold
            }

            color_plot2(
              embedding = rhdf5::h5read("stability.h5", paste(
                "feature_stability",
                "embedding_list",
                input$feature_type,
                input$feature_steps,
                sep = "/"
              )),
              color_info = used_matrix,
              plt_height = input$gene_height * ppi,
              plt_width = input$gene_width * ppi,
              # color_values = function(n) { paletteer::paletteer_c("grDevices::OrRd", n)},
              unique_values = unique_values,
              color_values = color_values,
              pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
              pt_size = input$gene_pt_size,
              text_size = input$gene_text_size,
              display_legend = TRUE
            )
            dev.off()
          }
        )

        output$download_metadata <-  shiny::downloadHandler(
          filename = function() {
            "feature_metadata.pdf"
          },
          content = function(file) {
            # ggplot2::ggsave(file, to_save_plot, width = width, height = height)
            shiny::req(input$feature_steps, input$metadata, input$metadata_width, input$metadata_height)
            pdf(file, width = input$metadata_width, height = input$metadata_height)
            metadata_plot(
              embedding = rhdf5::h5read("stability.h5", paste(
                "feature_stability",
                "embedding_list",
                input$feature_type,
                input$feature_steps,
                sep = "/"
              )),
              metadata_name = input$metadata,
              plt_height = input$metadata_height * ppi,
              plt_width = input$metadata_width * ppi,
              pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
              pt_size = input$metadata_pt_size,
              text_size = input$metadata_text_size,
              axis_size = input$metadata_axis_size,
              labels = input$metadata_labels,
              groups_highlight = input$select_groups,
              display_legend = TRUE
            )
            dev.off()
          }
        )


      # shiny::observe({
      #   download_plot_modal(
      #     session,
      #     output,
      #     paste0(paste("gene_plot", input$feature_type, input$feature_steps, paste(input$gene_expr, collapse = "_"), sep = "_"), ".pdf"),
      #     "download_gene"
      #   )
      # }) %>% shiny::bindEvent(input$download_gene_action)

      # shiny::observe({
      #   if (input$feature_type == "" || is.null(input$feature_steps)) {
      #     Sys.sleep(0.5) # wait a little
      #   }

      #   if (input$feature_steps == "" || is.null(input$feature_steps)) {
      #     Sys.sleep(0.5)
      #   }

      #   download_plot_modal(
      #     session,
      #     output,
      #     paste0(paste("metadata_plot", input$feature_type, input$feature_steps, input$metadata, sep = "_"), ".pdf"),
      #     "download_metadata"
      #   )



      #   print("output creat")
      # }) %>% shiny::bindEvent(input$download_metadata_action)

      # shiny::observe(
      #   {
      #     output$download_gene <- download_plot_handler(
      #       session,
      #       input$filename,
      #       umap_gene_plot(),
      #       input$width,
      #       input$height
      #     )

      #     output$download_metadata <- download_plot_handler(
      #       session,
      #       input$filename,
      #       umap_metadata_plot(),
      #       input$width,
      #       input$height
      #     )
      #   },
      #   priority = 10
      # ) %>% shiny::bindEvent(input$filename)
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

      shiny::observe(dr_choice_info(session)) %>% shiny::bindEvent(input$info_choice)

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

#' Dimensionality reduction - server side
#'
#' @description to be completed
#'
#' @export
server_dimensionality_reduction <- function(id, parent_session) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      update_sliders(session)

      print(paste(Sys.time(), "Loading the stability object"))
      # add_env_variable("stab_obj", rhdf5::h5read("stability.h5", "feature_stability"))
      print(paste(Sys.time(), "Finished loading"))

      shiny::observe(dr_title_info(session)) %>% shiny::bindEvent(input$info_title)
      shiny::observe(dr_comparison_info(session)) %>% shiny::bindEvent(input$info_comparison)
      # shiny::observe(
      server_dimensionality_stability("stability")
      # ) %>% shiny::bindEvent(input[["stability-resolution"]], ignoreInit = TRUE, once = TRUE)

      # server_dimensionality_distribution("distribution_left")
      # shiny::observe({
      #   print("Starts loading")
        server_dimensionality_distribution("distribution_left")
        server_dimensionality_distribution("distribution_right")
      # }) %>% shiny::bindEvent(input[["distribution_right-metadata"]], ignoreInit = TRUE, once = TRUE)
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

    shiny::updateSelectizeInput(
      session,
      inputId = glue::glue("distribution_{panels}-gene_expr"),
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

    shiny::updateSelectizeInput(
      session,
      inputId = glue::glue("distribution_{panels}-metadata"),
      server = FALSE,
      choices = colnames(pkg_env$metadata),
      selected = colnames(pkg_env$metadata)[1]
    )
  }
}



metadata_download_module <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

    }
  )
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
  predicted_width <- strwidth(fgroups, units = "inches", cex = text_size)
  predicted_height <- strheight(fgroups[1], units = "inches", cex = text_size)
  number_columns <- plt_width %/% (1.5 * max(predicted_width))
  number_rows <- ceiling(length(fgroups) / number_columns)
  current_margins <- par("mai")
  old_margin <- current_margins[1]
  current_margins[1] <- current_margins[1] + (number_rows + 2) * predicted_height[1] * 1.2 
  par(mai = current_margins)
  
  col <- rhdf5::h5read("stability.h5", "feature_stability/colours")
  n_groups <- length(fgroups)
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

  boundaries <- boxplot(
    formula = ecc ~ ftype + fsize,
    data = rhdf5::h5read("stability.h5", ifelse(resval == "overall",
      "feature_stability/overall/by_step",
      paste0("feature_stability/by_steps/", resval)
    )),
    col = col,
    at = at_values,
    xaxt = "n",
    xlab = "# features",
    boxwex = width * (n_groups + (space_intra_groups - 1) * (n_groups - 1)) / n_groups,
    outline = FALSE,
    cex.lab = text_size,
    cex.main = text_size * 1.2,
    cex.axis = text_size,
    main = plt_title
  )

  abline(v = abline_coords, lty = "dashed", col = "grey")

  y_text <- boundaries$stats[nrow(boundaries$stats), ]
  y_text[is.na(y_text)] <- 0

  text(
    x = at_values,
    y = y_text + text_offset,
    name_values,
    cex = text_size,
    xpd = NA
  )

  legend(
    "bottomleft",
    legend = fgroups,
    col = col,
    pch = 15,
    cex = text_size,
    pt.cex = text_size * 2,
    bty = 'n',
    ncol = number_columns,
    # text.width = NA,
    xpd = TRUE,
    inset = c(0, -((old_margin + predicted_height[1] * 2) / plt_height))
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
  predicted_width <- strwidth(fgroups, units = "inches", cex = text_size)
  predicted_height <- strheight(fgroups[1], units = "inches", cex = text_size)
  number_columns <- plt_width %/% (1.2 * max(predicted_width))
  number_rows <- ceiling(length(fgroups) / number_columns)
  current_margins <- par("mai")
  old_margin <- current_margins[1]
  current_margins[1] <- current_margins[1] + (number_rows + 2) * predicted_height[1] * 1.2
  par(mai = current_margins)

  col <- rhdf5::h5read("stability.h5", "feature_stability/colours")
  n_groups <- length(fgroups)
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

  # par(mai = c(0.1, 1, 0.5 + text_offset * plot_height / 96, 0.1))
  boundaries <- boxplot(
    formula = ecc ~ ftype + fsize,
    data = rhdf5::h5read("stability.h5", ifelse(resval == "overall",
      "feature_stability/overall/incremental",
      paste0("feature_stability/incremental/", resval)
    )),
    col = col,
    at = at_values,
    xaxt = "n",
    xlab = "# features",
    boxwex = width * (n_groups + (space_intra_groups - 1) * (n_groups - 1)) / n_groups,
    outline = FALSE,
    cex.lab = text_size,
    cex.main = text_size * 1.2,
    cex.axis = text_size,
    main = plt_title
  )

  abline(v = abline_coords, lty = "dashed", col = "grey")

  y_text <- boundaries$stats[nrow(boundaries$stats), ]
  y_text[is.na(y_text)] <- 0

  text(
    x = at_values,
    y = y_text + text_offset,
    name_values,
    cex = text_size,
    xpd = NA
  )

  legend(
    "bottomleft",
    legend = fgroups,
    col = col,
    pch = 15,
    cex = text_size,
    pt.cex = text_size * 2,
    bty = 'n',
    ncol = number_columns,
    # text.width = NA,
    xpd = TRUE,
    inset = c(0, -((old_margin + predicted_height[1] * 2) / plt_height))
  )
}
