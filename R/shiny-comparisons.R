####### UI #######

ui_comparisons <- function(id){
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Config Comparison",
    shiny::fluidRow(
      shiny::splitLayout(cellWidths = c('50%','50%'),
                         shiny::div(style="width:90%; overflow-wrap: break-word;",shiny::h1('Compare your current configuration',style="margin-bottom:10px "),
                                    shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info_UMAPs"), status = 'success', icon = shiny::icon('info'),size='sm'),
                                                           shiny::h3(shiny::strong('Compare your current configuration')),
                                                           shiny::div(style="white-space: pre-wrap; /* css-3 */
                                                                      white-space: -moz-pre-wrap; /* Mozilla, since 1999 */
                                                                      white-space: -pre-wrap; /* Opera 4-6 */
                                                                      white-space: -o-pre-wrap; /* Opera 7 */
                                                                      word-wrap: break-word; /* Internet Explorer 5.5+ */",
                                                                      shiny::h5('In this plot you can compare the configuration you have selected in the previous tabs to any other configuration. Here, you can also colour each configuration by the ECC, clusters, as well as any other metadata features that you have previously specified. You should choose the number of clusters that makes the most sense to you.')),
                                                           shiny::h5('For more information please go to:'),
                                                           shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href="https://github.com/Core-Bioinformatics/ClustAssess",target="_blank")),
                                                           placement = "right",
                                                           arrow = F,
                                                           maxWidth = '700px'),
                                    shiny::uiOutput(ns('k_selection_fixed')),
                                    shiny::verticalLayout(shiny::h1('Configuration 1'),
                                                          shiny::plotOutput(ns("umap_fixed")),
                                                          shiny::splitLayout(cellWidths = c('50%','50%'),
                                                                             shiny::div(style="width:90%;",shiny::verticalLayout(shiny::selectInput(ns("col_by_fixed"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), multiple = FALSE),
                                                                                                                                 shiny::selectInput(ns("col_by_fixed_2"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), selected='Clusters', multiple = FALSE))),
                                                                             shiny::div(style="width:90%;",shiny::verticalLayout(shiny::h1(),
                                                                                                                                 shinyWidgets::dropdownButton(
                                                                                                                                   label = "",
                                                                                                                                   icon = shiny::icon("download"),
                                                                                                                                   status = "success",
                                                                                                                                   size='sm',
                                                                                                                                   shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                                   shiny::textInput(ns("filename_umap_1"), "File name:", width = "80%"),
                                                                                                                                   shiny::numericInput(ns("width_umap_1"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                                   shiny::numericInput(ns("height_umap_1"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                                   shiny::selectInput(ns('umap_filetype_1'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                                   shiny::downloadButton(ns('download_umap_1'), label="Download Plot")
                                                                                                                                   
                                                                                                                                 ),
                                                                                                                                 shiny::h1(),
                                                                                                                                 shiny::h1(),
                                                                                                                                 shinyWidgets::dropdownButton(
                                                                                                                                   label = "",
                                                                                                                                   icon = shiny::icon("download"),
                                                                                                                                   status = "success",
                                                                                                                                   size='sm',
                                                                                                                                   shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                                   shiny::textInput(ns("filename_umap_3"), "File name:", width = "80%"),
                                                                                                                                   shiny::numericInput(ns("width_umap_3"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                                   shiny::numericInput(ns("height_umap_3"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                                   shiny::selectInput(ns('umap_filetype_3'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                                   shiny::downloadButton(ns('download_umap_3'), label="Download Plot")
                                                                                                                                   
                                                                                                                                 )))),
                                                          shiny::plotOutput(ns("umap_fixed_2"))),style="border-right:5px solid;"),
                         shiny::div(style="width:90%;",shiny::verticalLayout(shiny::splitLayout(cellWidths = c("50%","50%"),
                                                                                                shiny::verticalLayout(shiny::selectInput(ns("compare_sel_fset"), "Select feature - set", choices = names(stab_obj$feature_stability$by_steps), selected = names(stab_obj$feature_stability$by_steps)[1], multiple = FALSE),
                                                                                                                      shiny::selectInput(ns("compare_sel_steps"), "Select feature - size:", choices = stab_obj$feature_ordering$stable[[1]], selected = stab_obj$feature_ordering$stable[[1]][1], multiple = FALSE)),
                                                                                                shiny::verticalLayout(shiny::radioButtons(ns('clustering_method_choice'),'Select clustering Method',choices=names(stab_obj[[1]][[1]]$clustering_stability$split_by_k),selected=names(stab_obj[[1]][[1]]$clustering_stability$split_by_k)[1], width='50%'),
                                                                                                                      shiny::uiOutput(ns('k_selection')))),
                                                                             shiny::h1('Configuration 2'),
                                                                             shiny::plotOutput(ns("umap_choice")),
                                                                             shiny::splitLayout(cellWidths = c('33%','33%','33%'),
                                                                                                shiny::div(style="width:90%;",shiny::verticalLayout(shiny::selectInput(ns("col_by_choice"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), multiple = FALSE),
                                                                                                                                                    shiny::selectInput(ns("col_by_choice_2"), "Colour by:", choices = c('ECC',colnames(metadata), 'Clusters'), selected='Clusters', multiple = FALSE))),
                                                                                                shiny::div(style="width:90%;",shiny::verticalLayout(shiny::h1(),
                                                                                                                                                    shinyWidgets::dropdownButton(
                                                                                                                                                      label = "",
                                                                                                                                                      icon = shiny::icon("download"),
                                                                                                                                                      status = "success",
                                                                                                                                                      size='sm',
                                                                                                                                                      shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                                                      shiny::textInput(ns("filename_umap_2"), "File name:", width = "80%"),
                                                                                                                                                      shiny::numericInput(ns("width_umap_2"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::numericInput(ns("height_umap_2"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::selectInput(ns('umap_filetype_2'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                                                      shiny::downloadButton(ns('download_umap_2'), label="Download Plot")
                                                                                                                                                      
                                                                                                                                                    ),
                                                                                                                                                    shiny::h1(),
                                                                                                                                                    shiny::h1(),
                                                                                                                                                    shinyWidgets::dropdownButton(
                                                                                                                                                      label = "",
                                                                                                                                                      icon = shiny::icon("download"),
                                                                                                                                                      status = "success",
                                                                                                                                                      size='sm',
                                                                                                                                                      shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
                                                                                                                                                      shiny::textInput(ns("filename_umap_4"), "File name:", width = "80%"),
                                                                                                                                                      shiny::numericInput(ns("width_umap_4"), "Width (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::numericInput(ns("height_umap_4"), "Height (in):", 7, 3, 100, 0.1),
                                                                                                                                                      shiny::selectInput(ns('umap_filetype_4'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                                                                      shiny::downloadButton(ns('download_umap_4'), label="Download Plot")
                                                                                                                                                      
                                                                                                                                                    )))),
                                                                             shiny::plotOutput(ns("umap_choice_2")))
                         ),
      ),
      shiny::h1('Jaccard Simmilarity Index (JSI)/Cells per cluster'),
      shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info"), status = 'success', icon = shiny::icon('info'),size='sm'),
                             shiny::h3(shiny::strong('Jaccard Simmilarity Index (JSI) between clusters')),
                             shiny::br(),
                             shiny::h5('This plot aims to showcase the behaviour of the individual clusters on the different partitions. JSI is calculated for the cell barcodes for every cluster, in both configurations, in a pair-wise manner.'),
                             shiny::h1('\n'),
                             shiny::h5('For more information please go to:'),
                             shiny::tagList("", a("https://github.com/Core-Bioinformatics/ClustAssess", href="https://github.com/Core-Bioinformatics/ClustAssess",target="_blank")),
                             placement = "right",
                             arrow = F,
                             maxWidth = '700px'),
      shinyWidgets::dropdownButton(
        label = "",
        icon = shiny::icon("download"),
        status = "success",
        size='sm',
        shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
        shiny::textInput(ns("filename_heatmap"), "File name:", width = "80%"),
        shiny::numericInput(ns("width_heatmap"), "Width (in):", 7, 3, 100, 0.1),
        shiny::numericInput(ns("height_heatmap"), "Height (in):", 7, 3, 100, 0.1),
        shiny::selectInput(ns('heatmap_filetype'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
        shiny::downloadButton(ns('download_heatmap'), label="Download Plot")
        
      ),
      shinyWidgets::dropdownButton(
        label = "",
        icon = shiny::icon("cog"),
        status = "success",
        size='sm',
        shiny::radioButtons(ns('heatmap_type'),'Calculate similarity',choices=c('JSI','Cells per cluster'), width='100%')
      ),
      
      shiny::plotOutput(ns('barcode_heatmap'))
      
    ),style = "margin-left: 25px;")
}
####### SERVER #######

server_comparisons <- function(id,chosen_config,chosen_method){
  shiny::moduleServer(
    id,
    function(input, output, session) {
      #Set the dropdown menu depending on the available clusters
      ns <- shiny::NS(id)
      chosen_config <- stab_obj[[unlist(chosen_config[1])]][[unlist(chosen_config[2])]]
      embedding_fixed <- chosen_config$umap
      chosen_method <- chosen_config$clustering_stability$split_by_k[[chosen_method]]
      output$k_selection_fixed <- shiny::renderUI({
        shiny::selectInput(ns("k_fixed"), "Select a number of clusters:", choices = names(chosen_method), selected = names(chosen_method)[1], multiple = FALSE)
      })
      
      output$k_selection <- shiny::renderUI({
        shiny::selectInput(ns("k"), "Select a number of clusters:", choices = names(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k[[input$clustering_method_choice]]), multiple = FALSE)
      })
      
      #Set all UMAPs: 1,3 belong to the fixed selection, 2,4 are the UMAPs that change gepending on the users choice
      stable_graph_fixed <- paste0(toupper(chosen_config$stable_config$base_embedding),'_',chosen_config$stable_config$graph_type)
      stable_nn_fixed <- as.character(chosen_config$stable_config$n_neighbours)
      stable_ecc_fixed <- chosen_config$nn_stability$n_neigh_ec_consistency[[stable_graph_fixed]][[stable_nn_fixed]]
      
      umap_1 <- shiny::reactive({
        if (input$col_by_fixed=='ECC'){
          ECC <- stable_ecc_fixed
          ggplot2::ggplot(data.frame(
            chosen_config$umap),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_fixed=='Clusters'){
          Clusters <- as.factor(as.matrix(chosen_method[[input$k_fixed]]$mb))
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata[[input$col_by_fixed]]
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      umap_3 <- shiny::reactive({
        if (input$col_by_fixed_2=='ECC'){
          ECC <- stable_ecc_fixed
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_fixed_2=='Clusters'){
          if(is.null(input$k_fixed)) {
            return(ggplot2::ggplot() + ggplot2::theme_void())
          }else{
            Clusters <- as.factor(as.matrix(chosen_method[[input$k_fixed]]$mb))
            ggplot2::ggplot(data.frame(
              embedding_fixed),
              ggplot2::aes(x = .data$X1,
                           y = .data$X2,
                           color = Clusters)) +
              ggplot2::geom_point() +
              ggplot2::theme_bw()
          }
        }else{
          Feat <- metadata[[input$col_by_fixed_2]]
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      umap_2 <- shiny::reactive({
        if(is.null(input$k)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        if (input$col_by_choice=='ECC'){
          ECC <- stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$nn_stability$n_neigh_ec_consistency[[paste0(toupper(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$stable_config$base_embedding),'_',stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$stable_config$graph_type)]][[as.character(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$stable_config$n_neighbours)]]
          ggplot2::ggplot(data.frame(
            stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_choice=='Clusters'){
          Clusters <- as.factor(as.matrix(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k[[input$clustering_method_choice]][[input$k]]$mb))
          ggplot2::ggplot(data.frame(
            stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata[[input$col_by_choice]]
          ggplot2::ggplot(data.frame(
            stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      umap_4 <- shiny::reactive({
        if(is.null(input$k)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        if (input$col_by_choice_2=='ECC'){
          ECC <- stable_ecc_choice()
          ggplot2::ggplot(data.frame(
            stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_choice_2=='Clusters'){
          if(is.null(input$k)) {
            return(ggplot2::ggplot() + ggplot2::theme_void())
          }else{
            Clusters <- as.factor(as.matrix(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k[[input$clustering_method_choice]][[input$k]]$mb))
            ggplot2::ggplot(data.frame(
              stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
              ggplot2::aes(x = .data$X1,
                           y = .data$X2,
                           color = Clusters)) +
              ggplot2::geom_point() +
              ggplot2::theme_bw() 
          }
        }else{
          Feat <- metadata[[input$col_by_choice_2]]
          ggplot2::ggplot(data.frame(
            stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap),
            ggplot2::aes(x = .data$X1,
                         y = .data$X2,
                         color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      #Output all of the UMAPs
      output$umap_fixed <- shiny::renderPlot({
        umap_1()
      })
      output$umap_choice <- shiny::renderPlot({
        umap_2()
      })
      output$umap_fixed_2 <- shiny::renderPlot({
        umap_3()
      })
      output$umap_choice_2 <- shiny::renderPlot({
        umap_4()
      })
      
      umap_filetype_1 <- shiny::reactive({
        if (input$umap_filetype_1=='PDF'){
          filename <- paste0(input$filename_umap_1,'.pdf')
          return(filename)
        }else if (input$umap_filetype_1=='PNG'){
          filename <- paste0(input$filename_umap_1,'.png')
          return(filename)
        }else{
          filename <- paste0(input$filename_umap_1,'.svg')
          return(filename)
        }
      })
      output$download_umap_1 <- shiny::downloadHandler(
        filename = function() {umap_filetype_1()},
        content = function(file) {
          ggplot2::ggsave(file, umap_1(),width = input$width_umap_1,
                          height = input$height_umap_1,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      umap_filetype_2 <- shiny::reactive({
        if (input$umap_filetype_2=='PDF'){
          filename <- paste0(input$filename_umap_2,'.pdf')
          return(filename)
        }else if (input$umap_filetype_2=='PNG'){
          filename <- paste0(input$filename_umap_2,'.png')
          return(filename)
        }else{
          filename <- paste0(input$filename_umap_2,'.svg')
          return(filename)
        }
      })
      output$download_umap_2 <- shiny::downloadHandler(
        filename = function() {umap_filetype_2()},
        content = function(file) {
          ggplot2::ggsave(file, umap_2(),width = input$width_umap_2,
                          height = input$height_umap_2,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      umap_filetype_3 <- shiny::reactive({
        if (input$umap_filetype_3=='PDF'){
          filename <- paste0(input$filename_umap_3,'.pdf')
          return(filename)
        }else if (input$umap_filetype_3=='PNG'){
          filename <- paste0(input$filename_umap_3,'.png')
          return(filename)
        }else{
          filename <- paste0(input$filename_umap_3,'.svg')
          return(filename)
        }
      })
      output$download_umap_3 <- shiny::downloadHandler(
        filename = function() {umap_filetype_3()},
        content = function(file) {
          ggplot2::ggsave(file, umap_3(),width = input$width_umap_3,
                          height = input$height_umap_3,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      
      umap_filetype_4 <- shiny::reactive({
        if (input$umap_filetype_4=='PDF'){
          filename <- paste0(input$filename_umap_4,'.pdf')
          return(filename)
        }else if (input$umap_filetype_4=='PNG'){
          filename <- paste0(input$filename_umap_4,'.png')
          return(filename)
        }else{
          filename <- paste0(input$filename_umap_4,'.svg')
          return(filename)
        }
      })
      output$download_umap_4 <- shiny::downloadHandler(
        filename = function() {umap_filetype_4()},
        content = function(file) {
          ggplot2::ggsave(file, umap_4(),width = input$width_umap_4,
                          height = input$height_umap_4,
                          units = "in",
                          limitsize = FALSE)
        }
      )
      
      #JSI heatmap
      barcode_heatmap <- shiny::reactive({
        if(is.null(input$k)) {
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        df_fixed <- data.frame(as.matrix(chosen_method[[input$k_fixed]]$mb)) 
        df_fixed$cell <- rownames(df_fixed)
        clustering_choice <- as.matrix(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_stability$split_by_k[[input$clustering_method_choice]][[input$k]]$mb)
        df_choice <- data.frame(clustering_choice) 
        df_choice$cell <- rownames(df_choice)
        all_clusters_1 <- unique(df_fixed[,1])
        all_clusters_2 <- unique(df_choice[,1])
        
        mat = matrix(, nrow = length(all_clusters_2), 
                     ncol = length(all_clusters_1))
        colnames(mat) <- as.factor(sort(all_clusters_1))
        rownames(mat) <- as.factor(sort(all_clusters_2))
        if (input$heatmap_type=='JSI'){
          for (m in all_clusters_1){
            cluster_1 <- rownames(df_fixed[df_fixed[,1]==m,])
            for (n in all_clusters_2){
              cluster_2 <- rownames(df_choice[df_choice[,1]==n,])
              mat[as.character(n),as.character(m)] <- bulkAnalyseR::jaccard_index(cluster_1, cluster_2)
            }
          }
        }else{
          for (m in all_clusters_1){
            cluster_1 <- rownames(df_fixed[df_fixed[,1]==m,])
            for (n in all_clusters_2){
              cluster_2 <- rownames(df_choice[df_choice[,1]==n,])
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
            axis.title = ggplot2::element_text(size=14,face="bold"),
            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
            axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))) + 
          ggplot2::xlab("Clusters in Configuration 2") +
          ggplot2::ylab("Clusters in Configuration 1") +
          ggplot2::labs(fill="JSI")
      })
      output$barcode_heatmap <- shiny::renderPlot({
        barcode_heatmap()
      })
      heatmap_filetype <- shiny::reactive({
        if (input$heatmap_filetype=='PDF'){
          filename <- paste0(input$filename_heatmap,'.pdf')
          return(filename)
        }else if (input$heatmap_filetype=='PNG'){
          filename <- paste0(input$filename_heatmap,'.png')
          return(filename)
        }else{
          filename <- paste0(input$filename_heatmap,'.svg')
          return(filename)
        }
      })
      output$download_heatmap <- shiny::downloadHandler(
        filename = function() {heatmap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,barcode_heatmap(),width = input$width_heatmap,
                          height = input$height_heatmap,
                          units = "in",
                          limitsize = FALSE)
        }
      )
    }
  )
}
