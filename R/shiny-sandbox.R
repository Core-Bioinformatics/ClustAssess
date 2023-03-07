####### UI #######
ui_sandbox <- function(id){
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Sandbox",
    shiny::fluidRow(
      shiny::h1('Compare your current configuration'),
      shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info_Sandoboxs"), status = 'success', icon = shiny::icon('info'),size='sm'),
                             shiny::h3(shiny::strong('Compare any two configurations')),
                             shiny::div(style="white-space: pre-wrap; /* css-3 */
                                                                      white-space: -moz-pre-wrap; /* Mozilla, since 1999 */
                                                                      white-space: -pre-wrap; /* Opera 4-6 */
                                                                      white-space: -o-pre-wrap; /* Opera 7 */
                                                                      word-wrap: break-word; /* Internet Explorer 5.5+ */",
                                        shiny::h5('In this plot you can compare any configuration to another. Here, you can also colour each configuration by the ECC, clusters, as well as any other metadata features that you have previously specified. You should choose the number of clusters that makes the most sense to you. Please note that this tab should not be used to select a final configuration. We encourage you to go throught the inital tabs to find the most suitable configuration for your study.')),
                             shiny::h5('For more information please go to:'),
                             shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href="https://github.com/Core-Bioinformatics/ClustAssess",target="_blank")),
                             placement = "right",
                             arrow = F,
                             maxWidth = '700px'),
      shiny::splitLayout(cellWidths = c('50%','50%'),
                         shiny::div(style="width:90%;",shiny::verticalLayout(shiny::h1('Configuration 1'),
                                                                             shiny::splitLayout(cellWidths = c('50%','50%'),
                                                                                                shiny::verticalLayout(shiny::selectInput(ns("sandbox_sel_fset_1"), "Select feature - set", choices =names(stab_obj$feature_stability$by_steps), selected =names(stab_obj$feature_stability$by_steps)[1], multiple = FALSE),
                                                                                                                      shiny::selectInput(ns("sandbox_sel_steps_1"), "Select feature - size:", choices = stab_obj$feature_ordering$stable[[1]], selected = stab_obj$feature_ordering$stable[[1]][1], multiple = FALSE)),
                                                                                                shiny::verticalLayout(shiny::radioButtons(ns('sandbox_clustering_method_1'),'Select clustering Method',choices=names(stab_obj[[1]][[1]]$clustering_stability$split_by_k),selected=names(stab_obj[[1]][[1]]$clustering_stability$split_by_k)[1], width='50%'),
                                                                                                                      shiny::uiOutput(ns('sandbox_k_selection_1')))),
                                                                             shiny::plotOutput(ns("sandbox_umap_1")),
                                                                             shiny::splitLayout(cellWidths = c('50%','50%'),
                                                                                                shiny::div(style="width:90%;",shiny::verticalLayout(shiny::selectInput(ns("sandbox_col_1"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), multiple = FALSE),
                                                                                                                                                    shiny::selectInput(ns("sandbox_col_3"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), selected='Clusters', multiple = FALSE))),
                                                                                                shiny::verticalLayout(shiny::h1(),
                                                                                                                      shinyWidgets::dropdownButton(
                                                                                                                        label = "",
                                                                                                                        icon = shiny::icon("download"),
                                                                                                                        status = "success",
                                                                                                                        size='sm',
                                                                                                                        shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                        shiny::textInput(ns("sandbox_filename_umap_1"), "File name:", width = "80%"),
                                                                                                                        shiny::numericInput(ns("sandbox_width_umap_1"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                        shiny::numericInput(ns("sandbox_height_umap_1"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                        shiny::selectInput(ns('sandbox_umap_filetype_1'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                        shiny::downloadButton(ns('sandbox_download_umap_1'), label="Download Plot")
                                                                                                                        
                                                                                                                      ),
                                                                                                                      shiny::h1(),
                                                                                                                      shiny::h1(),
                                                                                                                      shinyWidgets::dropdownButton(
                                                                                                                        label = "",
                                                                                                                        icon = shiny::icon("download"),
                                                                                                                        status = "success",
                                                                                                                        size='sm',
                                                                                                                        shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                        shiny::textInput(ns("sandbox_filename_umap_3"), "File name:", width = "80%"),
                                                                                                                        shiny::numericInput(ns("sandbox_width_umap_3"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                        shiny::numericInput(ns("sandbox_height_umap_3"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                        shiny::selectInput(ns('sandbox_umap_filetype_3'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                        shiny::downloadButton(ns('sandbox_download_umap_3'), label="Download Plot")
                                                                                                                        
                                                                                                                      ))),
                                                                             shiny::plotOutput(ns("sandbox_umap_3"))),style="border-right:5px solid;"),
                         shiny::div(style="width:90%;",shiny::verticalLayout(shiny::uiOutput(ns('k_selection')),
                                                                             shiny::h1('Configuration 2'),
                                                                             shiny::splitLayout(cellWidths = c("50%","50%"),
                                                                                                shiny::verticalLayout(shiny::selectInput(ns("sandbox_sel_fset_2"), "Select feature - set", choices =names(stab_obj$feature_stability$by_steps), selected =names(stab_obj$feature_stability$by_steps)[1], multiple = FALSE),
                                                                                                                      shiny::selectInput(ns("sandbox_sel_steps_2"), "Select feature - size:", choices = stab_obj$feature_ordering$stable[[1]], selected = stab_obj$feature_ordering$stable[[1]][1], multiple = FALSE)),
                                                                                                shiny::verticalLayout(shiny::radioButtons(ns('sandbox_clustering_method_2'),'Select clustering Method',choices=names(stab_obj[[1]][[1]]$clustering_stability$split_by_k),selected=names(stab_obj[[1]][[1]]$clustering_stability$split_by_k)[1], width='50%'),
                                                                                                                      shiny::uiOutput(ns('sandbox_k_selection_2')))),
                                                                             shiny::plotOutput(ns("sandbox_umap_2")),
                                                                             shiny::splitLayout(cellWidths = c('33%','33%','33%'),
                                                                                                shiny::div(style="width:90%;",shiny::verticalLayout(shiny::selectInput(ns("sandbox_col_2"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), multiple = FALSE),
                                                                                                                                                    shiny::selectInput(ns("sandbox_col_4"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), selected='Clusters', multiple = FALSE))),
                                                                                                shiny::div(style="width:90%;",shiny::verticalLayout(shiny::h1(),
                                                                                                                                                    shinyWidgets::dropdownButton(
                                                                                                                                                      label = "",
                                                                                                                                                      icon = shiny::icon("download"),
                                                                                                                                                      status = "success",
                                                                                                                                                      size='sm',
                                                                                                                                                      shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                                                      shiny::textInput(ns("sandbox_filename_umap_2"), "File name:", width = "80%"),
                                                                                                                                                      shiny::numericInput(ns("sandbox_width_umap_2"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::numericInput(ns("sandbox_height_umap_2"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::selectInput(ns('sandbox_umap_filetype_2'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                                                      shiny::downloadButton(ns('sandbox_download_umap_2'), label="Download Plot")
                                                                                                                                                      
                                                                                                                                                    ),
                                                                                                                                                    shiny::h1(),
                                                                                                                                                    shiny::h1(),
                                                                                                                                                    shinyWidgets::dropdownButton(
                                                                                                                                                      label = "",
                                                                                                                                                      icon = shiny::icon("download"),
                                                                                                                                                      status = "success",
                                                                                                                                                      size='sm',
                                                                                                                                                      shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                                                      shiny::textInput(ns("sandbox_filename_umap_4"), "File name:", width = "80%"),
                                                                                                                                                      shiny::numericInput(ns("sandbox_width_umap_4"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::numericInput(ns("sandbox_height_umap_4"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::selectInput(ns('sandbox_umap_filetype_4'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                                                      shiny::downloadButton(ns('sandbox_download_umap_4'), label="Download Plot")
                                                                                                                                                      
                                                                                                                                                    )))),
                                                                             shiny::plotOutput(ns("sandbox_umap_4")))
                         ),
      ),
      shiny::h1('Jaccard Simmilarity Index (JSI)/Cells per cluster'),
      shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info"), status = 'success', icon = shiny::icon('info'),size='sm'),
                             shiny::h3(shiny::strong('Jaccard Simmilarity Index (JSI) between clusters')),
                             shiny::br(),
                             shiny::h5('This plot aims to showcase the behaviour of the individual clusters on the different partitions. JSI is calculated for the cell barcodes for every cluster, in both configurations, in a pair-wise manner.'),
                             shiny::h1('\n'),
                             shiny::h5('For more information please go to:'),
                             shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href="https://github.com/Core-Bioinformatics/ClustAssess",target="_blank")),
                             placement = "right",
                             arrow = F,
                             maxWidth = '700px'),
      shiny::verticalLayout(shinyWidgets::dropdownButton(
        label = "",
        icon = shiny::icon("download"),
        status = "success",
        size='sm',
        shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
        shiny::textInput(ns("sandbox_filename_heatmap"), "File name:", width = "80%"),
        shiny::numericInput(ns("sandbox_width_heatmap"), "Width (in):", 7, 3, 100, 0.1),
        shiny::numericInput(ns("sandbox_height_heatmap"), "Height (in):", 7, 3, 100, 0.1),
        shiny::selectInput(ns('sandbox_heatmap_filetype'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
        shiny::downloadButton(ns('sandbox_download_heatmap'), label="Download Plot")
        
      ),
      shinyWidgets::dropdownButton(
        label = "",
        icon = shiny::icon("cog"),
        status = "success",
        size='sm',
        shiny::radioButtons(ns('sandbox_heatmap_type'),'Calculate similarity',choices=c('JSI','Cells per cluster'))
      ),
      
      shiny::plotOutput(ns('sandbox_barcode_heatmap')))
    ),style = "margin-left: 25px;",
    shiny::tags$head(shiny::tags$style(shiny::HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              ")
    )))
}


####### SERVER #######

server_sandbox <- function(id){
  shiny::moduleServer(
    id,
    function(input, output, session) {
      #Set the dropdown menu depending on the available clusters
      ns <- shiny::NS(id)
      #Set the dropdown menu depending on the available clusters
      output$sandbox_k_selection_1 <- shiny::renderUI({
        shiny::selectInput(ns("sandbox_k_1"), "Select a number of clusters:", choices = names(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_1]]), multiple = FALSE)
      })
      output$sandbox_k_selection_2 <- shiny::renderUI({
        shiny::selectInput(ns("sandbox_k_2"), "Select a number of clusters:", choices = names(stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_2]]), multiple = FALSE)
      })
      #Set the chosen configurations
      sandbox_1 <- shiny::reactive({
        embedding <- stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$umap
      })
      sandbox_2 <- shiny::reactive({
        embedding <- stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$umap
      })
      s_umap_1 <- shiny::reactive({
        if(is.null(input$sandbox_k_1)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        if (input$sandbox_col_1=='ECC'){
          ECC <- stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$nn_stability$n_neigh_ec_consistency[[paste0(toupper(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$stable_config$base_embedding),'_',stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$stable_config$graph_type)]][[as.character(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$stable_config$n_neighbours)]]
          ggplot2::ggplot(data.frame(
            sandbox_1()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$sandbox_col_1=='Clusters'){
          Clusters <-as.factor(as.matrix(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_1]][[input$sandbox_k_1]]$mb))
          ggplot2::ggplot(data.frame(
            sandbox_1()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata[[input$sandbox_col_1]]
          ggplot2::ggplot(data.frame(
            sandbox_1()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      s_umap_2 <- shiny::reactive({
        if(is.null(input$sandbox_k_2)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        if (input$sandbox_col_2=='ECC'){
          ECC <- stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$nn_stability$n_neigh_ec_consistency[[paste0(toupper(stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$stable_config$base_embedding),'_',stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$stable_config$graph_type)]][[as.character(stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$stable_config$n_neighbours)]]
          ggplot2::ggplot(data.frame(
            sandbox_2()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$sandbox_col_2=='Clusters'){
          Clusters <-as.factor(as.matrix(stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_2]][[input$sandbox_k_2]]$mb))
          ggplot2::ggplot(data.frame(
            sandbox_2()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata[[input$sandbox_col_2]]
          ggplot2::ggplot(data.frame(
            sandbox_2()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      s_umap_3 <- shiny::reactive({
        if(is.null(input$sandbox_k_1)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        if (input$sandbox_col_3=='ECC'){
          ECC <- stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$nn_stability$n_neigh_ec_consistency[[paste0(toupper(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$stable_config$base_embedding),'_',stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$stable_config$graph_type)]][[as.character(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$stable_config$n_neighbours)]]
          ggplot2::ggplot(data.frame(
            sandbox_1()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$sandbox_col_3=='Clusters'){
          if(is.null(input$sandbox_k_1)) {
            return(ggplot2::ggplot() + ggplot2::theme_void())
          }else{
            Clusters <-as.factor(as.matrix(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_1]][[input$sandbox_k_1]]$mb))
            ggplot2::ggplot(data.frame(
              sandbox_1()),
              ggplot2::aes(x = .data$X1,
                           y = .data$X2,
                           color = Clusters)) +
              ggplot2::geom_point() +
              ggplot2::theme_bw() 
          }
        }else{
          Feat <- metadata[[input$sandbox_col_3]]
          ggplot2::ggplot(data.frame(
            sandbox_1()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      s_umap_4 <- shiny::reactive({
        if(is.null(input$sandbox_k_2)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        if (input$sandbox_col_4=='ECC'){
          ECC <- stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$nn_importance$n_neigh_ec_consistency[[2]][[1]]
          ggplot2::ggplot(data.frame(
            sandbox_2()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$sandbox_col_4=='Clusters'){
          if(is.null(input$sandbox_k_2)) {
            return(ggplot2::ggplot() + ggplot2::theme_void())
          }else{
            Clusters <- as.factor(as.matrix(stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_2]][[input$sandbox_k_2]]$mb))
            ggplot2::ggplot(data.frame(
              sandbox_2()),
              ggplot2::aes(x = .data$X1,
                           y = .data$X2,
                           color = Clusters)) +
              ggplot2::geom_point() +
              ggplot2::theme_bw() 
          }
        }else{
          Feat <- metadata[[input$sandbox_col_4]]
          ggplot2::ggplot(data.frame(
            sandbox_2()),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      
      output$sandbox_umap_1 <-shiny::renderPlot({
        s_umap_1()
      })
      output$sandbox_umap_2 <-shiny::renderPlot({
        s_umap_2()
      })
      output$sandbox_umap_3 <-shiny::renderPlot({
        if(is.null(s_umap_3())) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        s_umap_3()
      })
      output$sandbox_umap_4 <-renderPlot({
        if(is.null(s_umap_4())) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        s_umap_4()
      })
      
      
      sandbox_umap_filetype_1 <- shiny::reactive({
        if (input$sandbox_umap_filetype_1=='PDF'){
          filename <- paste0(input$sandbox_filename_umap_1,'.pdf')
          return(filename)
        }else if (input$sandbox_umap_filetype_1=='PNG'){
          filename <- paste0(input$sandbox_filename_umap_1,'.png')
          return(filename)
        }else{
          filename <- paste0(input$sandbox_filename_umap_1,'.svg')
          return(filename)
        }
      })
      output$sandbox_download_umap_1 <- shiny::downloadHandler(
        filename = function() {sandbox_umap_filetype_1()},
        content = function(file) {
          ggplot2::ggsave(file,s_umap_1(),width = input$sandbox_width_umap_1,
                          height = input$sandbox_height_umap_1,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      
      
      sandbox_umap_filetype_2 <- shiny::reactive({
        if (input$sandbox_umap_filetype_2=='PDF'){
          filename <- paste0(input$sandbox_filename_umap_2,'.pdf')
          return(filename)
        }else if (input$sandbox_umap_filetype_2=='PNG'){
          filename <- paste0(input$sandbox_filename_umap_2,'.png')
          return(filename)
        }else{
          filename <- paste0(input$sandbox_filename_umap_2,'.svg')
          return(filename)
        }
      })
      output$sandbox_download_umap_2 <- shiny::downloadHandler(
        filename = function() {sandbox_umap_filetype_2()},
        content = function(file) {
          ggplot2::ggsave(file,s_umap_2(),width = input$sandbox_width_umap_2,
                          height = input$sandbox_height_umap_2,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      
      
      sandbox_umap_filetype_3 <- shiny::reactive({
        if (input$sandbox_umap_filetype_3=='PDF'){
          filename <- paste0(input$sandbox_filename_umap_3,'.pdf')
          return(filename)
        }else if (input$sandbox_umap_filetype_3=='PNG'){
          filename <- paste0(input$sandbox_filename_umap_3,'.png')
          return(filename)
        }else{
          filename <- paste0(input$sandbox_filename_umap_3,'.svg')
          return(filename)
        }
      })
      output$sandbox_download_umap_3 <- shiny::downloadHandler(
        filename = function() {sandbox_umap_filetype_3()},
        content = function(file) {
          ggplot2::ggsave(file,s_umap_3(),width = input$sandbox_width_umap_3,
                          height = input$sandbox_height_umap_3,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      
      sandbox_umap_filetype_4 <- shiny::reactive({
        if (input$sandbox_umap_filetype_4=='PDF'){
          filename <- paste0(input$sandbox_filename_umap_4,'.pdf')
          return(filename)
        }else if (input$sandbox_umap_filetype_4=='PNG'){
          filename <- paste0(input$sandbox_filename_umap_4,'.png')
          return(filename)
        }else{
          filename <- paste0(input$sandbox_filename_umap_4,'.svg')
          return(filename)
        }
      })
      output$sandbox_download_umap_4 <- shiny::downloadHandler(
        filename = function() {sandbox_umap_filetype_4()},
        content = function(file) {
          ggplot2::ggsave(file,s_umap_4(),width = input$sandbox_width_umap_4,
                          height = input$sandbox_height_umap_4,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      
      
      #JSI heatmap
      sandbox_barcode_heatmap <- shiny::reactive({
        if(is.null(input$sandbox_k_2)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        clustering_1 <- as.matrix(stab_obj[[input$sandbox_sel_fset_1]][[input$sandbox_sel_steps_1]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_1]][[input$sandbox_k_1]]$mb)
        df_1 <- data.frame(clustering_1) 
        df_1$cell <- rownames(df_1)
        clustering_2 <- clustering <- as.matrix(stab_obj[[input$sandbox_sel_fset_2]][[input$sandbox_sel_steps_2]]$clustering_stability$split_by_k[[input$sandbox_clustering_method_2]][[input$sandbox_k_2]]$mb)
        df_2 <- data.frame(clustering_2) 
        df_2$cell <- rownames(df_2)
        all_clusters_1 <- unique(df_1[,1])
        all_clusters_2 <- unique(df_2[,1])
        
        mat = matrix(, nrow = length(all_clusters_2), 
                     ncol = length(all_clusters_1))
        colnames(mat) <- sort(all_clusters_1)
        rownames(mat) <- sort(all_clusters_2)
        if (input$sandbox_heatmap_type=='JSI'){
          for (m in all_clusters_1){
            cluster_1 <- rownames(df_1[df_1[,1]==m,])
            for (n in all_clusters_2){
              cluster_2 <- rownames(df_2[df_2[,1]==n,])
              mat[as.character(n),as.character(m)] <- bulkAnalyseR::jaccard_index(cluster_1, cluster_2)
            }
          }
        }else{
          for (m in all_clusters_1){
            cluster_1 <- rownames(df_1[df_1[,1]==m,])
            for (n in all_clusters_2){
              cluster_2 <- rownames(df_2[df_2[,1]==n,])
              mat[as.character(n),as.character(m)] <- length(intersect(cluster_1, cluster_2))
            }
          }
        }
        df_mat <- reshape::melt(mat)
        
        ggplot2::ggplot(df_mat, ggplot2::aes(X1, X2)) + 
          ggplot2::geom_tile(ggplot2::aes(fill = value)) + 
          ggplot2::geom_text(ggplot2::aes(fill = value, label = round(value, 2))) + 
          ggplot2::scale_fill_gradient2(low = scales::muted("darkred"), 
                                        mid = "white", 
                                        high = scales::muted("midnightblue"), 
                                        midpoint = 0) + 
          ggplot2::scale_x_continuous(breaks = pretty(df_mat$X1, n = length(all_clusters_2))) +
          ggplot2::scale_y_continuous(breaks = pretty(df_mat$X2, n = length(all_clusters_1))) +
          ggplot2::theme(
            panel.background=ggplot2::element_rect(fill="white"), 
            axis.text.x = ggplot2::element_text(hjust = 1,vjust=1,size = 10,face = "bold"),
            axis.text.y = ggplot2::element_text(size = 10,face = "bold"),
            axis.title=ggplot2::element_text(size=14,face="bold"),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))) + 
          xlab("Clusters in Configuration 2") +
          ylab("Clusters in Configuration 1") +
          labs(fill="JSI")
      })
      output$sandbox_barcode_heatmap <- shiny::renderPlot({
        sandbox_barcode_heatmap()
      })
      
      sandbox_heatmap_filetype <- shiny::reactive({
        if (input$sandbox_heatmap_filetype=='PDF'){
          filename <- 'heatmap.pdf'
          return(filename)
        }else if (input$sandbox_heatmap_filetype=='PNG'){
          filename <- 'heatmap.png'
          return(filename)
        }else{
          filename <- 'heatmap.svg'
          return(filename)
        }
      })
      output$sandbox_download_heatmap <- shiny::downloadHandler(
        filename = function() {sandbox_heatmap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,sandbox_barcode_heatmap(),width = 30,
                          height = 20,
                          units = "cm")
        }
      )
      sandbox_heatmap_filetype <- shiny::reactive({
        if (input$sandbox_heatmap_filetype=='PDF'){
          filename <- paste0(input$sandbox_filename_heatmap,'.pdf')
          return(filename)
        }else if (input$sandbox_heatmap_filetype=='PNG'){
          filename <- paste0(input$sandbox_filename_heatmap,'.png')
          return(filename)
        }else{
          filename <- paste0(input$sandbox_filename_heatmap,'.svg')
          return(filename)
        }
      })
      output$sandbox_download_heatmap <- shiny::downloadHandler(
        filename = function() {sandbox_heatmap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,sandbox_barcode_heatmap(),width = input$sandbox_width_heatmap,
                          height = input$sandbox_height_heatmap,
                          units = "in",
                          limitsize = FALSE)
        }
      )
    }
  )
}