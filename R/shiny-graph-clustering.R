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
    # shiny::uiOutput(ns("clustering_options")),
    shiny::radioButtons(
        inputId = ns("radio_cluster_method"),
        label = "Choose the clustering method for the downstream analysis:",
        choices = "",
        # selected = clustering_options[1],
        width = "100%"
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
    shiny::h3("Boxplot distribution of the resolution and ncluster-wise stability of the clustering methods"),
    shiny::tabsetPanel(
      id = ns("boxplot_tabset"),
      boxplot_settings(ns),
      shiny::tabPanel(
        "Vary by resolution",
        shiny::plotOutput(ns("boxplot_per_res"), height = "auto")
      ),
      shiny::tabPanel(
        "Vary by k",
        shiny::plotOutput(ns("boxplot_per_k"), height = "auto")
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
    shiny::fluidRow(
      shiny::column(1, gear_umaps(ns, "clustering_umap")),
      shiny::column(11,
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
    shiny::h3("Boxplot distribution of the overall stability"),
    shiny::splitLayout(
      shiny::plotOutput(ns("boxplot_overall_resolution"), height = "auto"),
      shiny::plotOutput(ns("boxplot_overall_k"), height = "auto")
    )
  )
}

ui_graph_clustering_k_stab <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Correspondence between the resolution value and the number of clusters"),
    shiny::uiOutput(ns("k_resolution_generator")),
    shiny::h3("The stability of the number of clusters"),
    shiny::uiOutput(ns("k_stability_generator")),
  )
}


ui_comparison_markers <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::selectInput(
      inputId = ns("select_method_markers"),
      label = "Select the clustering method",
      choices = ""
    ),
    # shiny::uiOutput(ns("select_n_clusters_generator")),
    shiny::selectInput(
      inputId = ns("select_k_markers"),
      label = "Select the number of clusters (k)",
      choices = ""
    ),
    shinyWidgets::pickerInput(
          inputId = ns("select_clusters_markers"),
          label = "Select the group of k",
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
    shiny::actionButton(ns("show_config"), "Show config", type = "info", class = "btn-info", disabled = TRUE),
    ui_graph_clustering_per_value_boxplot(ns("per_value_boxplot")),
    ui_graph_clustering_per_value_umap(ns("per_value_umap")),
    ui_graph_clustering_overall_boxplot(ns("overall_boxplot")),
    ui_graph_clustering_k_stab(ns("k_stab")),
    ui_graph_clustering_choice(ns("cluster_method_choice")),

    shiny::h2("Identification of markers"),
    shiny::fluidRow(
      shiny::column(
          6,
          ui_comparison_markers(ns("group_left"))
      ),
      shiny::column(
          6,
          ui_comparison_markers(ns("group_right"))
      )
  ),
  shiny::actionButton(ns("markers_button"),
      "Find markers!"
  ),
  DT::dataTableOutput(ns("markers")),
  shiny::downloadButton(ns("markers_download_button"),
      "Download markers!"
  ),
  )
}

###### SERVER ######
server_graph_clustering_choice <- function(id, parent_session) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      # clustering_options <- names(pkg_env$stab_obj$clustering_stability[[1]])
      clustering_options <- names(pkg_env$stab_obj$structure_list)
      # print(names(pkg_env$stab_ob()))
      shiny::updateRadioButtons(
        session = session,
        inputId = "radio_cluster_method",
        choices = clustering_options,
        selected = clustering_options[1],
      )

      user_choice <- shiny::reactive(input$radio_cluster_method) %>% shiny::bindEvent(input$fix_cluster_button)

      shiny::observe({
        shiny::showTab("tabset_id", "Comparison", select = TRUE, session = parent_session)
      }) %>% shiny::bindEvent(input$fix_cluster_button, ignoreInit = TRUE)

      return(user_choice)
    }
  )
}

server_graph_clustering_per_value_boxplot <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      output$boxplot_per_res <- shiny::renderPlot(
        height = function() {
          floor(pkg_env$height_ratio * pkg_env$dimension()[2])
        },
        {
          # shiny::req(pkg_env$lock_k())
          shiny_plot_clustering_per_value_stability(pkg_env$stab_obj$ecc_by_res,
            value_type = "resolution",
            width = input$boxplot_width,
            space_inter_groups = input$space_inter_groups,
            space_intra_groups = input$space_intra_groups,
            text_size = input$text_size
          )
        }
      )

      output$boxplot_per_k <- shiny::renderPlot(
        height = function() {
          floor(pkg_env$height_ratio * pkg_env$dimension()[2])
        },
        {
          shiny_plot_clustering_per_value_stability(pkg_env$stab_obj$ecc_by_k,
            value_type = "k",
            width = input$boxplot_width,
            space_inter_groups = input$space_inter_groups,
            space_intra_groups = input$space_intra_groups,
            text_size = input$text_size
          )
        }
      )

      shiny::outputOptions(output, "boxplot_per_k", suspendWhenHidden = FALSE)
      shiny::outputOptions(output, "boxplot_per_res", suspendWhenHidden = FALSE)
    }
  )
}

server_graph_clustering_overall_boxplot <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {

      output$boxplot_overall_resolution <- shiny::renderPlot(
        height = function() {
          floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2))
        },
        {
        # shiny::req(pkg_env$lock_k())
        shiny_plot_clustering_overall_stability(pkg_env$stab_obj$clustering_stability,
          value_type = "resolution"
        )
      })

      output$boxplot_overall_k <- shiny::renderPlot(
        height = function() {
          floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] / 2))
        },
        {
        # shiny::req(pkg_env$lock_k())
        
        shiny_plot_clustering_overall_stability(pkg_env$stab_obj$clustering_stability,
          value_type = "k"
        )
      })
    }
  )
}

server_graph_clustering_per_value_umap <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      # output$select_cluster_method_generator <- shiny::renderUI(
      shiny::updateSelectInput(
        session = session,
        inputId = "select_method",
        choices = names(pkg_env$stab_obj$structure_list),
        selected = names(pkg_env$stab_obj$structure_list)[1],
      )
      # )

      shiny::observe({
        shiny::req(input$select_method != "")
        # print(input$select_method)
        shiny::updateSelectInput(
          session = session,
          inputId = "select_k",
          choices = pkg_env$stab_obj$structure_list[[input$select_method]],
          selected = pkg_env$stab_obj$structure_list[[input$select_method]][1]
        )
      })

      shiny::observe({
        shiny::req(input$select_k != "")
        shinyWidgets::updatePickerInput(
          session = session,
          inputId = "select_clusters",
          choices = seq_len(as.integer(input$select_k))
        )
      }) %>% shiny::bindEvent(input$select_k)

      # output$select_n_clusters_generator <- shiny::renderUI({
      #   shiny::req(input$select_method)
      #   sorted_k <- stringr::str_sort(names(pkg_env$stab_obj$clustering_stability$split_by_k[[input$select_method]]), numeric = TRUE)
      #   shiny::selectInput(
      #     inputId = session$ns("select_k"),
      #     label = "Select the number of clusters (k)",
      #     choices = sorted_k,
      #     selected = sorted_k[1]
      #   )
      # })

      # render_plot_by_height("umap_k", session)
      # render_plot_by_height("umap_ecc", session)
      plt_height <- shiny::reactive(
        floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.5))
        
      )

      k_legend_height <- shiny::reactive({
        # shiny::req(input$select_k != "")
        shiny::req(input$select_k %in% pkg_env$stab_obj$structure_list[[input$select_method]])

        # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
        pdf(file = NULL, width = plt_height(), height = plt_height())
        unique_values <- seq_len(as.integer(input$select_k))

        predicted_width <- strwidth(c(" ", unique_values), units = "inches", cex = input$clustering_umap_text_size) * 96
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
          cex = input$clustering_umap_text_size) * ppi

        dev.off()

        return(text_height)
      })

      ecc_legend_height <- shiny::reactive({
        # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
        pdf(file = NULL, width = plt_height(), height = plt_height())
        par(mai = c(0.1, 0, 0.1, 0))
        text_height <- strheight("TE\nXT\n", units = "inches", cex = input$clustering_umap_text_size)
        dev.off()
        return((0.2 + text_height) * ppi)
      })

      output$umap_k <- shiny::renderPlot(
        height = function() {
          plt_height()
        },
        width = function() {
          plt_height()
        },
        {
        # shiny::req(pkg_env$lock_choice(), pkg_env$lock_k(), pkg_env$lock_stable)
        shiny::req(input$select_k %in% pkg_env$stab_obj$structure_list[[input$select_method]])

        color_plot2(
          embedding = pkg_env$stab_obj$umap,
          color_info = factor(rhdf5::h5read("stability.h5", paste(
            # input$
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
          labels = input$clustering_umap_labels,
          groups_highlight = input$select_clusters
        )
      })

         output$umap_k_legend <- shiny::renderPlot(
          height = function() {
            k_legend_height()
          },
          width = function() {
            plt_height()
          },
          {
            shiny::req(k_legend_height())
            # shiny::req(input$select_k %in% pkg_env$stab_obj$structure_list[[input$select_method]])
            unique_values <- seq_len(as.integer(input$select_k))
            if (!is.null(input$select_clusters)) {
              unique_values <- as.integer(input$select_clusters)
            }

            only_legend_plot(
              unique_values = unique_values,
              color_values = rhdf5::h5read("stability.h5", paste0("colors/", input$select_k))[unique_values],
              color_info = NULL,
              plt_width = plt_height(),
              text_size = input$clustering_umap_text_size 
            )

          }
        )

        ecc_value <- shiny::reactive({
          shiny::req(input$select_k != "")
          
          formatted_k <- sprintf("%06d", as.integer(input$select_k))
           rhdf5::h5read("stability.h5", paste(
              # input$
              pkg_env$lock_stable$feature_set,
              pkg_env$lock_stable$n_features,
              "clustering_stability",
              "split_by_k",
              "ecc",
              paste(formatted_k, input$select_method, sep = ";"),
              sep = "/"
            ))
        })

        output$umap_ecc <- shiny::renderPlot(
          height = function() {
            plt_height()
          },
          width = function() {
            plt_height()
          },
          {
          # shiny::req(pkg_env$lock_choice(), pkg_env$lock_k(), pkg_env$lock_stable)
          shiny::req(input$select_k != "")
          formatted_k <- sprintf("%06d", as.integer(input$select_k))

          ecc_order <- rhdf5::h5read("stability.h5", paste(
            # input$
            pkg_env$lock_stable$feature_set,
            pkg_env$lock_stable$n_features,
            "clustering_stability",
            "split_by_k",
            "ecc_order",
            paste(formatted_k, input$select_method, sep = ";"),
            sep = "/"
          ))

          color_plot2(
            embedding = pkg_env$stab_obj$umap[ecc_order, ],
            color_info = ecc_value(),
            color_values = NULL,
            unique_values = NULL,
            plt_height = plt_height(),
            plt_width = plt_height(),
            pch = ifelse(input$clustering_umap_pt_type == "Pixel", ".", 19),
            pt_size = input$clustering_umap_pt_size,
            text_size = input$clustering_umap_text_size
          )
        })

        output$umap_ecc_legend <- shiny::renderPlot(
          height = function() {
            ecc_legend_height()
          },
          width = function() {
            plt_height()
          },
          {
            # shiny::req(k_legend_height())
            shiny::req(input$select_k %in% pkg_env$stab_obj$structure_list[[input$select_method]])
            only_legend_plot(
              unique_values = NULL,
              color_values = NULL,
              color_info = ecc_value(),
              plt_width = plt_height(),
              text_size = input$clustering_umap_text_size 
            )

          }
        )

      # output$umap_ecc <- shiny::renderPlot({
      #   shiny::req(pkg_env$lock_choice(), pkg_env$lock_k())
      #   color_plot(
      #     pkg_env$stab_obj$umap,
      #     as.numeric(pkg_env$stab_obj$clustering_stability$split_by_k[[input$select_method]][[input$select_k]]$ecc),
      #     # color_palette = "viridis::viridis",
      #     ncolors = 50
      #   ) # + ggplot2::theme(
      #   #   legend.position = "bottom",
      #   #   # aspect.ratio = 1,
      #   #   legend.key.width = ggplot2::unit(floor(pkg_env$dimension()[1] / 16), "points")
      #   # ) +
      #   #   ggplot2::scale_color_viridis_c("ECC")
      # })
    }
  )
}

server_graph_clustering_k_stab <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      render_plot_by_height("k_resolution", session)
      render_plot_by_height("k_stability", session)

      output$k_resolution <- shiny::renderPlot({
        # shiny::req(pkg_env$lock_k())
        shiny_plot_k_resolution_corresp(pkg_env$stab_obj$clustering_stability,
          colour_information = "ecc"
        )
      })

      output$k_stability <- shiny::renderPlot({
        # shiny::req(pkg_env$lock_k())
        shiny_plot_k_n_partitions(pkg_env$stab_obj$clustering_stability,
          colour_information = "ecc"
        )
      })
    }
  )
}


server_comparison_markers <- function(id) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      shiny::updateSelectInput(
        session = session,
        inputId = "select_method_markers",
        choices = names(pkg_env$stab_obj$structure_list),
        selected = names(pkg_env$stab_obj$structure_list)[1],
      )
      # )

      shiny::observe({
        shiny::req(input$select_method_markers %in% names(pkg_env$stab_obj$structure_list))
        # print(input$select_method)
        shiny::updateSelectInput(
          session = session,
          inputId = "select_k_markers",
          choices = pkg_env$stab_obj$structure_list[[input$select_method_markers]],
          selected = pkg_env$stab_obj$structure_list[[input$select_method_markers]][1]
        )
      })

      shiny::observe({
        shiny::req(input$select_k_markers %in% pkg_env$stab_obj$structure_list[[input$select_method_markers]])
        shinyWidgets::updatePickerInput(
          session = session,
          inputId = "select_clusters_markers",
          choices = seq_len(as.integer(input$select_k_markers))
        )
      }) %>% shiny::bindEvent(input$select_k_markers)

    }
  )
}

#' Graph clustering - server side
#'
#' @description to be completed
#'
#' @export
server_graph_clustering <- function(id, feature_choice, parent_session) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      isolated_fchoice <- shiny::isolate(feature_choice())
      ftype <- isolated_fchoice$chosen_feature_type
      fsize <- isolated_fchoice$chosen_set_size

      # add_env_variable("dimension", window_height)

      print(paste(Sys.time(), "Loading the stability object"))
      if (exists("stable_config")) {
        stable_config <- NULL
      }
      add_env_variable("stab_obj", list(
        ecc_by_k = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "ecc", sep = "/")),
        ecc_by_res = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_resolution", "ecc", sep = "/")),
        structure_list = rhdf5::h5read("stability.h5", paste(ftype, fsize, "clustering_stability", "split_by_k", "structure_list", sep = "/")),
        umap = rhdf5::h5read("stability.h5", paste(ftype, fsize, "umap", sep = "/"))
      ))
      stable_config <- rhdf5::h5read("stability.h5", paste(ftype, fsize, "stable_config", sep = "/"))
      add_env_variable("lock_choice", shiny::reactive(input$"cluster_method_choice-radio_cluster_method"))
      add_env_variable("lock_k", shiny::reactive(input$"per_value_umap-select_k"))
      add_env_variable("lock_stable", stable_config)
      print(paste(Sys.time(), "Finished loading"))

      shiny::observe(
        shinyjs::enable("show_config")
      ) %>% shiny::bindEvent(input$"cluster_method_choice-radio_cluster_method", ignoreInit = TRUE, once = TRUE)

      clustering_choice <- server_graph_clustering_choice("cluster_method_choice", parent_session)


      


      # shiny::observe({
      server_graph_clustering_per_value_umap("per_value_umap")
      # }) %>% shiny::bindEvent(input$"cluster_method_choice-radio_cluster_method", ignoreInit = TRUE, once = TRUE)

      # shiny::observe({
      # shiny::req(input$"per_value_umap-select_k")
      server_graph_clustering_per_value_boxplot("per_value_boxplot")
      # server_graph_clustering_overall_boxplot("overall_boxplot")
      # server_graph_clustering_k_stab("k_stab")
      # }) %>% shiny::bindEvent(input$"per_value_umap-select_k", ignoreInit = TRUE, once = TRUE)

      shiny::observe({
        # shiny::req(input$"cluster_method_choice-radio_cluster_method")
        shiny::showModal(
          stable_config_info(stable_config),
          session = session
        )
      }) %>% shiny::bindEvent(input$show_config)

  
  
      server_comparison_markers("group_left")
      server_comparison_markers("group_right")
      used_genes <- rhdf5::h5read("stability.h5", paste(ftype, "feature_list", sep = "/"))[seq_len(as.numeric(fsize))]
      expr_matrix <- rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(pkg_env$genes_of_interest[used_genes], NULL))
      rank_matrix <- rhdf5::h5read("expression.h5", "rank_of_interest", index = list(pkg_env$genes_of_interest[used_genes], NULL))
      rownames(expr_matrix) <- used_genes

      markers_val <- shiny::reactive({
        print("start")
        print(system.time({
          mb1 <-  factor(rhdf5::h5read("stability.h5", paste(
            # input$
            ftype,
            fsize,
            "clustering_stability",
            "split_by_k",
            "mbs",
            input$"group_left-select_method_markers",
            input$"group_left-select_k_markers",
            sep = "/"
          )))
        
      mb2 <-  factor(rhdf5::h5read("stability.h5", paste(
            # input$
            ftype,
            fsize,
            "clustering_stability",
            "split_by_k",
            "mbs",
            input$"group_right-select_method_markers",
            input$"group_right-select_k_markers",
            sep = "/"
          )))

      cells_index_left <- which(mb1 %in% input$"group_left-select_clusters_markers" )
      cells_index_right <- which(mb2 %in% input$"group_right-select_clusters_markers")

      }))
      # cells_left <- colnames(expr_matrix())[cells_index_left]
      # cells_right <- colnames(expr_matrix())[cells_index_right]
      

      calculate_markers(
        expression_matrix = expr_matrix,
        cells1 = cells_index_left,
        cells2 = cells_index_right,
        rank_matrix = rank_matrix
      )

      # print(mks)
      # write.csv(mks, "test.csv")

      # return(mks)

      }) %>% shiny::bindEvent(input$markers_button)


    output$markers <- DT::renderDataTable({
      shiny::req(markers_val())
      markers_val()
    }, rownames = FALSE)

    output$markers_download_button <- shiny::downloadHandler(
      filename = function() { "markers.csv" },
      content = function(file) {
         write.csv(markers_val(), file)
      }
    )

    shinyjs::hide("markers_download_button")

    shiny::observe({
      shiny::req(markers_val())
      shinyjs::show("markers_download_button")
    }) %>% shiny::bindEvent(markers_val())


      # gc()
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
    status = "info",
    size = "sm",
    icon = shiny::icon("gear")
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
  split_names <- lapply(names(ecc_list), function(x) { strsplit(x, ";")[[1]]})
  k_or_res_values <- sapply(split_names, function(x) { x[1] })
  unique_values <- as.numeric(unique(k_or_res_values))
  cl_method <- sapply(split_names, function(x) { x[2] })
  unique_cl_methods <- unique(cl_method)
  n_methods <- length(unique_cl_methods)
  cl_method <- match(cl_method, unique_cl_methods)
  
  at_values <- rep(0, length(cl_method))
  text_coords <- rep(0, length(unique_values))
  abline_coords <- rep(0, length(unique_values)-1)
  
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
      abline_coords[count_diff] <- (at_values[i] + at_values[i-1]) / 2
      updated_count <- FALSE
      text_coords[count_diff] <- mean(at_values[start_index:(i-1)])
      start_index <- i
    }
  }
  
  text_coords[count_diff+1] <- mean(at_values[start_index:length(cl_method)])
  
  
  boxplot(
    ecc_list,
    outline = FALSE,
    at = at_values,
    col = rhdf5::h5read("stability.h5", paste0("colors/", n_methods))[cl_method],
    boxwex = width * (n_methods + (space_intra_groups - 1) * (n_methods -1)) / n_methods,
    xaxt = "n",
    xlab = NA,
    cex.axis = text_size,
    cex.lab = text_size
  )
  abline(v = abline_coords, lty = "dashed", col = "grey")
  title(xlab = "k", ylab = "ecc", cex.lab = text_size)
  axis(side = 1, at = text_coords, labels = unique_values, las = 2, cex.axis = text_size)

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

  colorsi <- rainbow(max_levels, s = 0.5)


  return({
    boxplot(
      ecc ~ method + x_value,
      data = melted_df,
      col = colorsi,
      pch = ".",
      outline = FALSE,
      at = at_values,
      xaxt = "n",
      boxwex = width * (max_levels + (space_intra_groups - 1) * (max_levels-1)) / max_levels,
      xlab = NA,
      cex.axis = text_size,
      cex.lab = text_size
    )
    axis(side = 1, at = middle_values, labels = cluster_values, las = 2, cex.axis = text_size)
    title(xlab = value_type, cex.lab = text_size)
    abline(v = abline_coords, lty = "dashed", col = "grey")
  })
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

  colorsi <- rainbow(nlevels(melted_df$method), s = 0.5)
  return({
    boxplot(
      ecc ~ method,
      data = melted_df,
      col = colorsi,
      pch = "."
    )
  })

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
    ) +
    ggplot2::theme(legend.position = "bottom")

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
    ) +
    ggplot2::theme(legend.position = "bottom")



  return(main_plot)
}
