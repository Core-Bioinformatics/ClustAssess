#### UI ####
ui_dimensionality_stability <- function(id) {
  ns <- shiny::NS(id)
  plot_settings <- shinyWidgets::dropdownButton(
    shiny::selectInput(
      inputId = ns("feature_types"),
      label = "Feature names",
      choices = feature_types,
      selected = feature_types,
      multiple = TRUE
    ),
    shinyWidgets::sliderTextInput(
      inputId = ns("resolution"),
      label = "Resolution",
      choices = feature_used_resolution_vals,
      selected = feature_used_resolution_vals[1]
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
      cellWidths = c("55%", "40%"),
      shiny::plotOutput(ns("boxplot_ecc"),
        hover = shiny::hoverOpts(id = ns("ecc_hover"), delayType = "throttle"),
        click = ns("ecc_click"),
        height = "500px",
      ),
      shiny::fluidRow(
        shiny::textOutput(ns("boxplot_ecc_info")),
        shiny::tableOutput(ns("table_ecc_info")),
        shiny::plotOutput(ns("umap_ecc"), height = "350px"),
        width = "40%"
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

ui_dimensionality_distribution_plots <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::selectInput(ns("feature_type"), "Feature names", feature_types),
    shiny::selectInput(ns("feature_steps"), "Feature set size", NULL),
    shiny::hr(),
    shiny::fluidRow(
      shiny::column(
        6,
        shiny::selectInput(
          inputId = ns("gene_expr"),
          label = "Gene name(s)",
          choices = rownames(expression_matrix)[1:3000],
          selected = rownames(expression_matrix)[1],
          multiple = TRUE
        )
      ),
      shiny::column(
        6,
        shiny::conditionalPanel(
          condition = "input.gene_expr.length > 1",
          shiny::sliderInput(
            inputId = ns("expr_threshold"),
            label = "Gene expression threshold",
            min = 0, max = 10, value = 0.1
          ),
          ns = ns
        )
      )
    ),
    shiny::actionButton(ns("download_gene_action"), "Download", icon = shiny::icon("download")),
    shiny::plotOutput(ns("umap_gene")),
    shiny::selectInput(
      inputId = ns("metadata"),
      label = "Metadata",
      choices = colnames(metadata),
      selected = colnames(metadata)[1]
    ),
    shiny::actionButton(ns("download_metadata_action"), "Download", icon = shiny::icon("download")),
    shiny::plotOutput(ns("umap_metadata"))
  )
}

ui_dimensionality_recommendation <- function(id) {
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
    shiny::radioButtons(
      ns("radio_feature_type"),
      label = "Choose the feature type for the downstream analysis:",
      choices = feature_types,
      selected = feature_types[1],
      width = "80%",
    ),
    shiny::radioButtons(
      inputId = ns("radio_feature_size"),
      label = "Choose the size of the feature set for the downstream analysis:",
      choices = names(stab_obj[[feature_types[1]]]),
      selected = NULL,
      width = "80%"
    ),
    shiny::uiOutput(ns("recommendation")),
    shiny::actionButton(ns("fix_feature_button"),
      "Fix the configuration!",
      style = "font-size:20px;",
      class = "btn-danger"
    ),
    style = "padding:50px; font-size:20px;"
  )
}


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
    shiny::fluidRow(
      shiny::column(
        6,
        ui_dimensionality_distribution_plots(ns("distribution_left")),
        style = "margin-bottom:30px;border-right:5px solid;padding:10px;"
      ),
      shiny::column(
        6,
        ui_dimensionality_distribution_plots(ns("distribution_right")),
        style = "margin-bottom:30px;padding:10px;"
      )
    ),
    ui_dimensionality_recommendation(ns("recommend")),
    style = "margin-bottom:30px;"
  )
}


#### SERVER ####

server_dimensionality_stability <- function(id, panel_id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      print(shiny::outputOptions(output))
      selected_configs <- shiny::reactive(input$feature_types)

      ecc_mouse_hover <- shiny::reactiveVal(c(NA, NA))
      shiny::observeEvent(input$ecc_hover, {
        if (is.null(input$ecc_hover)) {
          return()
        }
        x <- input$ecc_hover$x
        # y <- input$ecc_hover$y
        min_step <- min(as.numeric(boxplot_ecc_df()$step_index))
        max_step <- max(as.numeric(boxplot_ecc_df()$step_index))
        nconfigs <- length(selected_configs())
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

        ecc_mouse_hover(c(chosen_step, chosen_config))
      })

      # ecc_mouse_click <- shiny::reactiveVal(c(NA, NA))

      ecc_mouse_click <- shiny::reactive({
        if (is.null(input$ecc_click)) {
          return(c(NA, NA))
        }
        ecc_mouse_hover()
      })

      output$boxplot_ecc_info <- shiny::renderText({
        if (is.na(ecc_mouse_hover()[1])) {
          return("")
        }

        feature_type <- names(current_feature_obj()$steps_stability)[ecc_mouse_hover()[2]]
        nfeatures <- names(current_feature_obj()$steps_stability[[ecc_mouse_hover()[2]]])[ecc_mouse_hover()[1]]
        glue::glue("The ECC summary for {feature_type} - {nfeatures} and resolution {input$resolution}.")
        # Click if you want to visualise the UMAP distribution.")
      })

      output$table_ecc_info <- shiny::renderTable({
        if (is.na(ecc_mouse_hover()[1])) {
          return()
        }
        summary_ecc <- fivenum(current_feature_obj()$steps_stability[[ecc_mouse_hover()[2]]][[ecc_mouse_hover()[1]]][[input$resolution]]$ecc)
        data.frame(
          min = summary_ecc[1],
          Q1 = summary_ecc[2],
          median = summary_ecc[3],
          Q3 = summary_ecc[4],
          max = summary_ecc[5]
        )
      })

      # shiny::observe({
        output$stepchoosing <- shiny::renderUI({
          print("apelat pentru")
          print(input$feature_types)
          ns <- session$ns
          i <- 1
          ret_list <- list()
          for (conf in input$feature_types) {
            ret_list[[i]] <- shiny::selectInput(ns(paste0("nsteps_", i)),
              paste("Nsteps for", conf),
              choices = names(stab_obj$feature_importance$steps_stability[[conf]]),
              selected = names(stab_obj$feature_importance$steps_stability[[conf]]),
              multiple = TRUE
            )
            i <- i + 1
          }
          ret_list
        })
      # })# %>% shiny::bindEvent(input$feature_types)

      current_feature_obj <- shiny::reactive({
        if (is.null(selected_configs())[1]) {
          return(NULL)
        }

        for (i in seq_along(selected_configs())) {
          if (is.null(input[[glue::glue("nsteps_{i}")]])[1]) {
            return(NULL)
          }
        }

        steps_stability <- lapply(seq_along(selected_configs()), function(i) {
          steps_names <- names(stab_obj$feature_importance$steps_stability[[i]])
          chosen_steps <- steps_names %in% input[[paste0("nsteps_", i)]]
          chosen_steps <- steps_names[chosen_steps]
          stab_obj$feature_importance$steps_stability[[i]][chosen_steps]
        })
        names(steps_stability) <- selected_configs()

        incremental_stability <- lapply(seq_along(selected_configs()), function(i) {
          kept_steps <- c()
          chosen_steps <- input[[paste0("nsteps_", i)]]

          if (length(chosen_steps) == 0) {
            return(list())
          }

          for (j in seq_len(length(chosen_steps) - 1)) {
            checked_steps <- paste(chosen_steps[j], chosen_steps[j + 1], sep = "-")
            if (checked_steps %in% names(stab_obj$feature_importance$incremental_stability[[i]])) {
              kept_steps <- c(kept_steps, checked_steps)
            }
          }

          stab_obj$feature_importance$incremental_stability[[i]][kept_steps]
        })
        names(incremental_stability) <- selected_configs()

        list(
          steps_stability = steps_stability,
          incremental_stability = incremental_stability
        )
      })

      boxplot_ecc_df <- shiny::reactive({
        plot_feature_per_resolution_stability_boxplot(current_feature_obj(),
          resolution = input$resolution,
          return_df = TRUE
        )
      })

      output$boxplot_ecc <- shiny::renderPlot({
        if (length(selected_configs()) == 0) {
          return(empty_ggplot())
        }

        if (is.null(current_feature_obj())) {
          return(empty_ggplot())
        }

        plot_feature_per_resolution_stability_boxplot(current_feature_obj(),
          resolution = input$resolution,
          text_size = input$text_size,
          boxplot_width = input$boxplot_width,
          dodge_width = input$dodge_width
        ) + ggplot2::theme(legend.position = "bottom")
      })

      output$boxplot_incr <- shiny::renderPlot({
        if (length(selected_configs()) == 0) {
          return(empty_ggplot())
        }

        if (is.null(current_feature_obj())) {
          return(empty_ggplot())
        }

        for (conf in selected_configs()) {
          if (length(names(current_feature_obj()$incremental_stability[[conf]])) == 0) {
            return(empty_ggplot())
          }
        }
        plot_feature_per_resolution_stability_incremental(current_feature_obj(),
          resolution = input$resolution,
          text_size = input$text_size,
          boxplot_width = input$boxplot_width,
          dodge_width = input$dodge_width
        )
      })

      output$umap_ecc <- shiny::renderPlot({
        if (is.na(ecc_mouse_click()[1])) {
          return(empty_ggplot())
        }

        if (is.null(current_feature_obj())) {
          return(empty_ggplot())
        }

        config_name <- selected_configs()[ecc_mouse_click()[2]]
        step_index <- names(current_feature_obj()$steps_stability[[config_name]])[ecc_mouse_click()[1]]
        ecc <- current_feature_obj()$steps_stability[[config_name]][[step_index]][[input$resolution]]$ecc

        color_ggplot(stab_obj[[1]]$embedding_list[[config_name]][[step_index]], ecc) +
          # ggplot2::guides(color = ggplot2::guide_legend(title = "ecc")) +
          ggplot2::scale_color_viridis_c(name = "ECC") +
          ggplot2::ggtitle(glue::glue("ECC distribution {config_name} - {step_index} for resolution {input$resolution}"))
      })

      output$overall_boxplot_ecc <- shiny::renderPlot({
        if (length(selected_configs()) == 0) {
          return(empty_ggplot())
        }

        if (is.null(current_feature_obj())) {
          return(empty_ggplot())
        }

        plot_feature_overall_stability_boxplot(current_feature_obj(),
          text_size = input$text_size,
          boxplot_width = input$boxplot_width,
          dodge_width = input$dodge_width
        )
      })

      output$overall_boxplot_incremental <- shiny::renderPlot({
        if (length(selected_configs()) == 0) {
          return(empty_ggplot())
        }

        if (is.null(current_feature_obj())) {
          return(empty_ggplot())
        }

        plot_feature_overall_stability_incremental(current_feature_obj(),
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
        shiny::updateSelectInput(session,
          inputId = "feature_steps",
          choices = names(stab_obj[[1]]$steps_stability[[input$feature_type]])
        )
      })

      shiny::observeEvent(input$gene_expr, {
        shiny::updateSliderInput(session,
          inputId = "expr_threshold",
          max = round(max(expression_matrix[input$gene_expr, ]), 3),
          step = round(max(expression_matrix[input$gene_expr, ]) / 10, 3)
        )
      })

      umap_gene_plot <- shiny::reactive({
        if (is.null(input$feature_type) || is.null(input$feature_steps)) {
          return(empty_ggplot())
        }

        embedding <- stab_obj$feature_importance$embedding_list[[input$feature_type]][[input$feature_steps]]

        if (length(input$gene_expr) == 0) {
          return(empty_ggplot())
        }

        if (is.null(embedding)) {
          return(empty_ggplot())
        }

        if (length(input$gene_expr) == 1) {
          expr_val <- expression_matrix[input$gene_expr, ]

          return(expression_ggplot(embedding, expr_val))
        }

        voting_scheme_ggplot(embedding, input$gene_expr, input$expr_threshold)
      })

      umap_metadata_plot <- shiny::reactive({
        if (is.null(input$feature_type) || is.null(input$feature_steps)) {
          return(empty_ggplot())
        }

        embedding <- stab_obj$feature_importance$embedding_list[[input$feature_type]][[input$feature_steps]]

        if (is.null(embedding)) {
          return(empty_ggplot())
        }

        metadata_ggplot(embedding, input$metadata)
      })

      output$umap_gene <- shiny::renderPlot({
        umap_gene_plot()
      })
      output$umap_metadata <- shiny::renderPlot(umap_metadata_plot())

      shiny::observe({
        download_modal(
          session,
          paste0(paste("gene_plot", input$feature_type, input$feature_steps, paste(input$gene_expr, collapse = "_"), sep = "_"), ".pdf"),
          "download_gene"
        )
      }) %>% shiny::bindEvent(input$download_gene_action)

      shiny::observe({
        download_modal(
          session,
          paste0(paste("metadata_plot", input$feature_type, input$feature_steps, input$metadata, sep = "_"), ".pdf"),
          "download_metadata"
        )
      }) %>% shiny::bindEvent(input$download_metadata_action)

      shiny::observe(
        output$download_gene <- download_handler(
          session,
          input$filename,
          umap_gene_plot(),
          input$width,
          input$height
        )
      ) %>% shiny::bindEvent(input$filename)

      shiny::observe(
        output$download_metadata <- download_handler(
          session,
          input$filename,
          umap_metadata_plot(),
          input$width,
          input$height
        )
      ) %>% shiny::bindEvent(input$filename)
    }
  )
}

server_dimensionality_recommendation <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      shiny::observeEvent(input$radio_feature_type, {
        shiny::updateRadioButtons(
          session,
          label = glue::glue("Choose the size of the feature set {input$radio_feature_type} for the downstream analysis: (We recommend {names(stab_obj[[input$radio_feature_type]])[1]})"),
          inputId = "radio_feature_size",
          choices = names(stab_obj[[input$radio_feature_type]]),
          selected = names(stab_obj[[input$radio_feature_type]])[1]
        )
      })

      shiny::observe(dr_choice_info(session)) %>% shiny::bindEvent(input$info_choice)

      v <- shiny::reactive(list(
        chosen_feature_type = input$radio_feature_type,
        chosen_set_size = input$radio_feature_size
      )) %>% shiny::bindEvent(input$fix_feature_button)

      return(v)
    }
  )
}

server_dimensionality_reduction <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      shiny::observe(dr_title_info(session)) %>% shiny::bindEvent(input$info_title)
      shiny::observe(dr_comparison_info(session)) %>% shiny::bindEvent(input$info_comparison)
      server_dimensionality_stability("stability", "dim_reduc")
      server_dimensionality_distribution("distribution_left")
      server_dimensionality_distribution("distribution_right")
      recommendation <- server_dimensionality_recommendation("recommend")

      recommendation
    }
  )
}
