### UI ###
ui_graph_clustering <- function(id) {
  ns <- shiny::NS(id)

  shiny::tabPanel(
    "Graph Clustering",
    shiny::uiOutput(ns("selected_conf")),
    shiny::uiOutput(ns("header1")),
    shiny::fluidRow(
      shiny::column(
        8,
        shiny::plotOutput(ns("boxplot_per_value")),
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
    shiny::plotOutput(ns("k_stability"))
  )
}

### SERVER ###

server_graph_clustering <- function(id, recommendation) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      chosen_object <- shiny::reactive(
        stab_obj[[recommendation()[[1]]]][[recommendation()[[2]]]]
      )

      stable_config <- shiny::reactive(
        chosen_object()$stable_config
      )

      clustering_importance <- shiny::reactive(
        chosen_object()$clustering_importance
      )

      output$selected_conf <- shiny::renderUI({
        shiny::h2("Please fix your choice of feature type and the size of the set from the 'Dimensionality reduction' tab!")
      })

      shiny::observe({
        output$selected_conf <- shiny::renderUI({
          shiny::fluidRow(
            shiny::div(glue::glue("Feature type: {stable_config()[[1]]}")),
            div(glue::glue("Feature set size: {stable_config()[[2]]}")),
            div(glue::glue("Base embedding: {stable_config()[[3]]}")),
            div(glue::glue("Graph type: {stable_config()[[4]]}")),
            div(glue::glue("Number of neighbours: {stable_config()[[5]]}"))
          )
        }) %>% shiny::bindEvent(stable_config())

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
      }) %>% shiny::bindEvent(recommendation())

      output$boxplot_per_value <- shiny::renderPlot({
        plot_clustering_per_value_stability(clustering_importance(),
          value_type = "resolution"
        )
      })

      output$boxplot_overall_resolution <- shiny::renderPlot({
        plot_clustering_overall_stability(clustering_importance(),
          value_type = "resolution"
        )
      })

      output$boxplot_overall_k <- shiny::renderPlot({
        plot_clustering_overall_stability(clustering_importance(),
          value_type = "k"
        )
      })

      output$k_resolution <- shiny::renderPlot({
        plot_k_resolution_corresp(clustering_importance(),
          colour_information = "ecc"
        )
      })

      output$k_stability <- shiny::renderPlot({
        plot_k_n_partitions(clustering_importance(),
          colour_information = "ecc"
        )
      })
    }
  )
}
