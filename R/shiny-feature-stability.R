#### UI ####
ui_dimensionality_stability <- function(id) {
  ns <- shiny::NS(id)
  plot_settings <- shinyWidgets::dropdownButton(
    shiny::uiOutput(ns("resolution_slider")),
    shinyWidgets::sliderTextInput(
      inputId = ns("resolution"),
      label = "Resolution",
      choices = c("0.8")
    ),
    shiny::sliderInput(
      inputId = ns("boxplot_width"),
      label = "Boxplot width",
      min = 0.00, max = 1.00, value = 0.40
    ),
    shiny::sliderInput(
      inputId = ns("dodge_width"),
      label = "Dodge width",
      min = 0.00, max = 1.00, value = 0.70
    ),
    shiny::sliderInput(
      inputId = ns("text_size"),
      label = "Text size",
      min = 0.10, max = 25.00, value = 4.00
    ),
    circle = TRUE,
    status = "info",
    size = "sm",
    icon = shiny::icon("gear")
  )

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
    plot_settings,
    shiny::splitLayout(
      cellWidths = c("52%", "45%"),
      shiny::plotOutput(ns("boxplot_ecc"),
        hover = shiny::hoverOpts(id = ns("ecc_hover"), delayType = "throttle"),
        click = ns("ecc_click"),
        height = "500px"
      ),
      shiny::div(
        style = "width: 90%",
        shiny::textOutput(ns("boxplot_ecc_info")),
        shiny::tableOutput(ns("table_ecc_info")),
        shiny::plotOutput(ns("umap_ecc"), height = "350px", width = "90%")
      )
    ),
    shiny::splitLayout(
      cellWidths = c("40px", "90%"),
      shinyWidgets::circleButton(ns("info_ecs_incr_res"),
        icon = shiny::icon("info"),
        size = "sm",
        status = "success"
      ),
      shiny::h2("Incremental ECS per individual resolution values")
    ),
    shiny::fluidRow(
      shiny::column(
        8,
        shiny::plotOutput(ns("boxplot_incr"))
      ),
      shiny::column(
        4,
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
    shiny::plotOutput(ns("overall_boxplot_ecc")),
    shiny::splitLayout(
      cellWidths = c("40px", "90%"),
      shinyWidgets::circleButton(ns("info_ecs_incremental_overall"),
        icon = shiny::icon("info"),
        size = "sm",
        status = "success"
      ),
      shiny::h2("Overall incremental stability")
    ),
    shiny::plotOutput(ns("overall_boxplot_incremental"))
  )
}

ui_dimensionality_distribution_plots <- function(id, draw_line) {
  ns <- shiny::NS(id)
  style <- ifelse(draw_line, "border-right:5px solid", "")

  shiny::wellPanel(
    style = style,
    shiny::uiOutput(ns("test_output")),
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
        min = 0, max = 10, value = 0.1,
        width = "95%"
      )
    ),
    shiny::actionButton(ns("download_gene_action"), "Download", icon = shiny::icon("download")),
    shiny::uiOutput(ns("umap_gene_generator")),
    # shiny::plotOutput(ns("umap_gene"), height = "1000px"),
    shiny::selectizeInput(
      inputId = ns("metadata"),
      label = "Metadata",
      choices = NULL
    ),
    shiny::actionButton(ns("download_metadata_action"), "Download", icon = shiny::icon("download")),
    shiny::uiOutput(ns("umap_metadata_generator"))
    # shiny::plotOutput(ns("umap_metadata"), height = "auto")
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
      boxplot_ecc_df <- shiny::reactive(
        plot_feature_per_resolution_stability_boxplot(
          feature_object_list = pkg_env$stab_obj,
          resolution = input$resolution,
          return_df = TRUE
        )
      )

      # ecc_mouse_hover <- shiny::reactiveVal(c(NA, NA))
      ecc_mouse_hover <- shiny::reactive({
        if (is.null(input$ecc_hover)) {
          return(c(NA, NA))
        }
        x <- input$ecc_hover$x
        # y <- input$ecc_hover$y
        min_step <- min(as.numeric(boxplot_ecc_df()$step_index))
        max_step <- max(as.numeric(boxplot_ecc_df()$step_index))
        nconfigs <- length(pkg_env$feature_ordering$original)
        chosen_step <- NULL
        chosen_config <- NULL

        for (steps_intervals in seq(from = min_step - 0.5, to = max_step - 0.5, by = 1)) {
          if (x < steps_intervals + 1) {
            chosen_step <- steps_intervals + 0.5
            break
          }
        }

        if (is.null(chosen_step)) {
          chosen_step <- max_step
        }

        for (i in seq_len(nconfigs)) {
          if (x < (steps_intervals + i / nconfigs)) {
            chosen_config <- i
            break
          }
        }

        if (is.null(chosen_config)) {
          chosen_config <- nconfigs
        }

        return(c(chosen_step, chosen_config))
      }) %>% shiny::bindEvent(input$ecc_hover)


      output$boxplot_ecc_info <- shiny::renderText({
        if (is.na(ecc_mouse_hover()[1])) {
          return("")
        }

        feature_type <- pkg_env$feature_types[ecc_mouse_hover()[2]]
        nfeatures <- pkg_env$feature_ordering$original[[feature_type]][ecc_mouse_hover()[1]]
        glue::glue("The ECC summary for {feature_type} - {nfeatures} and resolution {input$resolution}.")
        # Click if you want to visualise the UMAP distribution.")
      })

      output$table_ecc_info <- shiny::renderTable({
        if (is.na(ecc_mouse_hover()[1])) {
          return()
        }
        summary_ecc <- fivenum(pkg_env$stab_obj$by_steps[[ecc_mouse_hover()[2]]][[ecc_mouse_hover()[1]]][[input$resolution]]$ecc)
        data.frame(
          min = summary_ecc[1],
          Q1 = summary_ecc[2],
          median = summary_ecc[3],
          Q3 = summary_ecc[4],
          max = summary_ecc[5]
        )
      })

      shiny::observe({
        output$boxplot_ecc <- shiny::renderPlot({
          plot_feature_per_resolution_stability_boxplot(
            feature_object_list = pkg_env$stab_obj,
            resolution = input$resolution,
            text_size = input$text_size,
            boxplot_width = input$boxplot_width,
            dodge_width = input$dodge_width
          ) + ggplot2::theme(legend.position = "bottom")
        })

        output$boxplot_incr <- shiny::renderPlot({
          plot_feature_per_resolution_stability_incremental(
            feature_object_list = pkg_env$stab_obj,
            resolution = input$resolution,
            text_size = input$text_size,
            boxplot_width = input$boxplot_width,
            dodge_width = input$dodge_width
          )
        })
      }) %>% shiny::bindEvent(input$resolution)


      shiny::observe({
        hover_coords <- ecc_mouse_hover()

        output$umap_ecc <- shiny::renderPlot({
          # if (is.na(hover_coords[1])) {
          #   return(empty_ggplot())
          # }

          config_name <- pkg_env$feature_types[hover_coords[2]]
          step_index <- pkg_env$feature_ordering$original[[config_name]][hover_coords[1]]
          ecc <- pkg_env$stab_obj$by_steps[[config_name]][[step_index]][[input$resolution]]$ecc

          color_ggplot(pkg_env$stab_obj$embedding_list[[config_name]][[step_index]], ecc) +
            # ggplot2::guides(color = ggplot2::guide_legend(title = "ecc")) +
            ggplot2::scale_color_viridis_c(name = "ECC") +
            ggplot2::ggtitle(sprintf("ECC distribution\n%s - %s\nResolution %s", config_name, step_index, input$resolution))
        })
      }) %>% shiny::bindEvent(input$ecc_click)

      output$overall_boxplot_ecc <- shiny::renderPlot({
        plot_feature_overall_stability_boxplot(
          feature_object_list = pkg_env$stab_obj,
          text_size = input$text_size,
          boxplot_width = input$boxplot_width,
          dodge_width = input$dodge_width
        )
      })

      output$overall_boxplot_incremental <- shiny::renderPlot({
        plot_feature_overall_stability_incremental(
          feature_object_list = pkg_env$stab_obj,
          text_size = input$text_size,
          boxplot_width = input$boxplot_width,
          dodge_width = input$dodge_width
        )
      })

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

      shiny::observe({
        height <- floor(input$dimension * pkg_env$height_ratio)
        ns <- session$ns
        output$umap_gene_generator <- shiny::renderUI(
          shiny::plotOutput(ns("umap_gene"), height = paste0(height, "px"))
        )

        output$umap_metadata_generator <- shiny::renderUI(
          shiny::plotOutput(ns("umap_metadata"), height = paste0(height, "px"))
        )
      }) %>% shiny::bindEvent(input$dimension)

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

      umap_gene_plot <- shiny::reactive({
        if (is.null(input$feature_type) || is.null(input$feature_steps)) {
          Sys.sleep(0.5) # wait a little
        }

        embedding <- pkg_env$stab_obj$embedding_list[[input$feature_type]][[input$feature_steps]]

        if (length(input$gene_expr) == 0) {
          return(empty_ggplot())
        }

        # if (is.null(embedding)) {
        #   return(empty_ggplot())
        # }

        if (length(input$gene_expr) == 1) {
          return(expression_ggplot(embedding, expr_matrix(), input$expr_threshold) + ggplot2::theme(aspect.ratio = 1))
        }

        voting_scheme_ggplot(
          embedding,
          matrixStats::colSums2(expr_matrix() >= input$expr_threshold),
          nrow(expr_matrix())
        ) + ggplot2::theme(aspect.ratio = 1)
      })

      umap_metadata_plot <- shiny::reactive({
        if (is.null(input$feature_type) || is.null(input$feature_steps)) {
          Sys.sleep(0.5) # wait a little
        }

        embedding <- pkg_env$stab_obj$embedding_list[[input$feature_type]][[input$feature_steps]]

        # if (is.null(embedding)) {
        #   return(empty_ggplot())
        # }

        metadata_ggplot(embedding, input$metadata) + ggplot2::theme(aspect.ratio = 1)
      })

      output$umap_gene <- shiny::renderCachedPlot(
        {
          umap_gene_plot()
        },
        cacheKeyExpr = {
          list(input$feature_type, input$feature_steps, input$gene_expr, input$expr_threshold)
        },
        sizePolicy = shiny::sizeGrowthRatio(
          height = 500, width = 500, growthRate = 1.2
        )
      )

      output$umap_metadata <- shiny::renderCachedPlot(
        {
          umap_metadata_plot()
        },
        cacheKeyExpr = {
          list(input$feature_type, input$feature_steps, input$metadata)
        },
        sizePolicy = shiny::sizeGrowthRatio(
          height = 500, width = 500, growthRate = 1.2
        )
      )


      shiny::observe({
        download_plot_modal(
          session,
          output,
          paste0(paste("gene_plot", input$feature_type, input$feature_steps, paste(input$gene_expr, collapse = "_"), sep = "_"), ".pdf"),
          "download_gene"
        )
      }) %>% shiny::bindEvent(input$download_gene_action)

      shiny::observe({
        if (input$feature_type == "" || is.null(input$feature_steps)) {
          Sys.sleep(0.5) # wait a little
        }

        if (input$feature_steps == "" || is.null(input$feature_steps)) {
          Sys.sleep(0.5)
        }

        download_plot_modal(
          session,
          output,
          paste0(paste("metadata_plot", input$feature_type, input$feature_steps, input$metadata, sep = "_"), ".pdf"),
          "download_metadata"
        )



        print("output creat")
      }) %>% shiny::bindEvent(input$download_metadata_action)

      shiny::observe(
        {
          output$download_gene <- download_plot_handler(
            session,
            input$filename,
            umap_gene_plot(),
            input$width,
            input$height
          )

          output$download_metadata <- download_plot_handler(
            session,
            input$filename,
            umap_metadata_plot(),
            input$width,
            input$height
          )
        },
        priority = 10
      ) %>% shiny::bindEvent(input$filename)
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
      add_env_variable("stab_obj", rhdf5::h5read("stability.h5", "feature_stability"))

      shiny::observe(dr_title_info(session)) %>% shiny::bindEvent(input$info_title)
      shiny::observe(dr_comparison_info(session)) %>% shiny::bindEvent(input$info_comparison)
      shiny::observe(
        server_dimensionality_stability("stability")
      ) %>% shiny::bindEvent(input[["stability-resolution"]], ignoreInit = TRUE, once = TRUE)

      shiny::observe({
        server_dimensionality_distribution("distribution_left")
        server_dimensionality_distribution("distribution_right")
      }) %>% shiny::bindEvent(input[["distribution_right-metadata"]], ignoreInit = TRUE, once = TRUE)
      feature_choice <- server_dimensionality_choice("feature_choice", parent_session)

      feature_choice
    }
  )
}

update_sliders <- function(session) {
  shinyWidgets::updateSliderTextInput(
    session,
    "stability-resolution",
    choices = pkg_env$feature_ordering$resolution[[1]][[1]], # NOTE might need an updated if there are different resolution values between feature types
    selected = pkg_env$feature_ordering$resolution[[1]][[1]][1]
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
