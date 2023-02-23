ui_comparisons <- function(id){
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Config Comparison",
    shiny::fluidRow(
      shiny::splitLayout(cellWidths = c('50%','50%'),
                         shiny::div(style="width:90%;",shiny::h1('Compare your current configuration',style="margin-bottom:50px "),
                                    shiny::uiOutput(ns('k_selection_fixed')),
                                    shiny::verticalLayout(shiny::h1('Configuration 1'),
                                                          shiny::plotOutput(ns("umap_fixed")),
                                                          shiny::splitLayout(cellWidths = c('50%','50%'),
                                                                             shiny::div(style="width:90%;",shiny::verticalLayout(shiny::selectInput(ns("col_by_fixed"), "Colour by:", choices = c('ECC',colnames(metadata_feats), 'Clusters'), multiple = FALSE),
                                                                                                                                 shiny::selectInput(ns("col_by_fixed_2"), "Colour by:", choices = c('ECC',colnames(metadata_feats), 'Clusters'), selected='Clusters', multiple = FALSE))),
                                                                             shiny::div(style="width:90%;",shiny::verticalLayout(shiny::downloadButton(ns('download_umap_1'), label="",
                                                                                                                                                       style = "margin: 25px;"),
                                                                                                      shiny::downloadButton(ns('download_umap_3'), label="",
                                                                                                                            style = "margin: 25px;")))),
                                                          shiny::plotOutput(ns("umap_fixed_2"))),style="border-right:5px solid;"),
                         shiny::div(style="width:90%;",shiny::verticalLayout(shiny::splitLayout(cellWidths = c("50%","50%"),
                                                                                                shiny::verticalLayout(shiny::selectInput(ns("compare_sel_fset"), "Select feature - set", choices = stab_obj_feature_config, selected = stab_obj_feature_config[1], multiple = FALSE),
                                                                                                                      shiny::selectInput(ns("compare_sel_steps"), "Select feature - size:", choices = names(stab_obj[[stab_obj_feature_config[1]]]), selected = names(stab_obj[[stab_obj_feature_config[1]]])[1], multiple = FALSE)),
                                                                                                shiny::verticalLayout(shiny::radioButtons(ns('clustering_method_choice'),'Select clustering Method',choices=stab_obj_clustering_methods,selected=stab_obj_clustering_methods[1], width='50%'),
                                                                                                                      shiny::uiOutput(ns('k_selection')))),
                                                                             shiny::h1('Configuration 2'),
                                                                             shiny::plotOutput(ns("umap_choice")),
                                                                             shiny::splitLayout(cellWidths = c('33%','33%','33%'),
                                                                                                shiny::div(style="width:90%;",shiny::verticalLayout(shiny::selectInput(ns("col_by_choice"), "Colour by:", choices = c('ECC',colnames(metadata_feats), 'Clusters'), multiple = FALSE),
                                                                                                                                                    shiny::selectInput(ns("col_by_choice_2"), "Colour by:", choices = c('ECC',colnames(metadata_feats), 'Clusters'), selected='Clusters', multiple = FALSE))),
                                                                                                shiny::div(style="width:90%;",shiny::verticalLayout(shiny::downloadButton(ns('download_umap_2'), label="",
                                                                                                                         style = "margin: 25px;"),
                                                                                                                         shiny::downloadButton(ns('download_umap_4'), label="",
                                                                                                                         style = "margin: 25px;"))),
                                                                                                shiny::selectInput(ns('filetype'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='40%')),
                                                                             shiny::plotOutput(ns("umap_choice_2")))
                  ),
      ),
      shiny::h1('Jaccard Simmilarity Index (JSI) between clusters'),
      shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info"), status = 'success', icon = icon('info'),size='sm'),
                             shiny::h3(shiny::strong('Jaccard Simmilarity Index (JSI) between clusters')),
                             shiny::br(),
                             shiny::h5('This plot aims to showcase the behaviour of the individual clusters on the different partitions. JSI is calculated for the cell barcodes for every cluster, in both configurations, in a pair-wise manner.'),
                             shiny::h1('\n'),
                             shiny::h5('For more information please go to:'),
                             shiny::tagList("", github),
               placement = "right",
               arrow = F,
               maxWidth = '700px'),
      shiny::verticalLayout(shiny::splitLayout(cellWidths = c('15%','85%'),
                                               shiny::radioButtons(ns('heatmap_type'),'Calculate similarity',choices=c('JSI','Cells per cluster'), width='50%'),
                                               shiny::div(style="width:90%;",shiny::splitLayout(cellWidths = c('15%','85%'),
                                                                                                shiny::selectInput(ns('heatmap_filetype'),'Filetype',choices = c('PDF','PNG','SVG'),selected='PDF',width='100%'),
                                                                                                shiny::downloadButton(ns('download_heatmap'), label="",
                                                                                   style = "margin: 25px; ")))),
                            shiny::plotOutput(ns('barcode_heatmap')))
    ),style = "margin-left: 25px;",
    tags$head(tags$style(shiny::HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              "))))
}




####### SERVER #######

server_comparisons <- function(id,chosen_config,chosen_method){
  shiny::moduleServer(
    id,
    function(input, output, session) {
      #Set the dropdown menu depending on the available clusters
      ns <- shiny::NS(id)
      embedding_fixed <- chosen_config$umap
      chosen_method <- chosen_config$clustering_importance$split_by_k[[chosen_method]]
      output$k_selection_fixed <- shiny::renderUI({
        shiny::selectInput(ns("k_fixed"), "Select a number of clusters:", choices = names(chosen_method), multiple = FALSE)
      })
      output$k_selection <- shiny::renderUI({
        shiny::selectInput(ns("k"), "Select a number of clusters:", choices = names(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_importance$split_by_k[[input$clustering_method_choice]]), multiple = FALSE)
      })
      #Set the user's choice for comaprison
      choice <- shiny::reactive({
        embedding_choice <- stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$umap
        clustering_choice <- as.matrix(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_importance$split_by_k[[input$clustering_method_choice]][[input$k]]$partitions[[1]]$mb)
        combo <- list(embedding_choice,clustering_choice)
      })
      #Set all UMAPs: 1,3 belong to the fixed selection, 2,4 are the UMAPs that change gepending on the users choice
      umap_1 <- shiny::reactive({
        if (input$col_by_fixed=='ECC'){
          ECC <- chosen_config$nn_importance$n_neigh_ec_consistency[[2]][[1]]
          ggplot2::ggplot(data.frame(
            chosen_config$umap),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_fixed=='Clusters'){
          Clusters <- as.matrix(as.matrix(chosen_method[[input$k_fixed]]$partitions[[1]]$mb))
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata_feats[[input$col_by_fixed]]
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      umap_2 <- shiny::reactive({
        if(is.null(input$k)) {
          return(ggplot2::ggplot() + theme_void())
        }
        if (input$col_by_choice=='ECC'){
          ECC <- stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$nn_importance$n_neigh_ec_consistency[[2]][[1]]
          ggplot2::ggplot(data.frame(
            choice()[1]),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_choice=='Clusters'){
          Clusters <- as.matrix(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_importance$split_by_k[[input$clustering_method_choice]][[input$k]]$partitions[[1]]$mb)
          ggplot2::ggplot(data.frame(
            choice()[1]),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata_feats[[input$col_by_choice]]
          ggplot2::ggplot(data.frame(
            choice()[1]),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      umap_3 <- shiny::reactive({
        if (input$col_by_fixed_2=='ECC'){
          ECC <- chosen_config$nn_importance$n_neigh_ec_consistency[[2]][[1]]
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_fixed_2=='Clusters'){
          Clusters <- as.matrix(as.matrix(chosen_method[[input$k_fixed]]$partitions[[1]]$mb))
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata_feats[[input$col_by_fixed_2]]
          ggplot2::ggplot(data.frame(
            embedding_fixed),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = Feat)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }
      })
      umap_4 <- shiny::reactive({
        if(is.null(input$k)) {
          return(ggplot() + theme_void())
        }
        if (input$col_by_choice_2=='ECC'){
          ECC <- stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$nn_importance$n_neigh_ec_consistency[[2]][[1]]
          ggplot2::ggplot(data.frame(
            choice()[1]),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = ECC)) +
            ggplot2::geom_point() +
            ggplot2::scale_color_viridis_c() +
            ggplot2::theme_bw()
        }else if(input$col_by_choice_2=='Clusters'){
          Clusters <- as.matrix(stab_obj[[input$compare_sel_fset]][[input$compare_sel_steps]]$clustering_importance$split_by_k[[input$clustering_method_choice]][[input$k]]$partitions[[1]]$mb)
          ggplot2::ggplot(data.frame(
            choice()[1]),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
                color = Clusters)) +
            ggplot2::geom_point() +
            ggplot2::theme_bw()
        }else{
          Feat <- metadata_feats[[input$col_by_choice_2]]
          ggplot2::ggplot(data.frame(
            choice()[1]),
            ggplot2::aes(x = .data$UMAP_1,
                y = .data$UMAP_2,
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
      
      #Options to download the UMAPs in different formats
      umap_filetype <- shiny::reactive({
        if (input$filetype=='PDF'){
          filename <- 'umap.pdf'
          return(filename)
        }else if (input$filetype=='PNG'){
          filename <- 'umap.png'
          return(filename)
        }else{
          filename <- 'umap.svg'
          return(filename)
        }
      })
      output$download_umap_1 <- shiny::downloadHandler(
        filename = function() {umap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,umap_1(),width = 20,
                 height = 20,
                 units = "cm")
        }
      )
      output$download_umap_2 <- shiny::downloadHandler(
        filename = function() {umap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,umap_2(),width = 20,
                 height = 20,
                 units = "cm")
        }
      )
      output$download_umap_3 <- shiny::downloadHandler(
        filename = function() {umap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,umap_3(),width = 20,
                 height = 20,
                 units = "cm")
        }
      )
      output$download_umap_4 <- shiny::downloadHandler(
        filename = function() {umap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,umap_4(),width = 20,
                 height = 20,
                 units = "cm")
        }
      )
      
      #JSI heatmap
      barcode_heatmap <- shiny::reactive({
        if(is.null(input$k)) {
          return(ggplot() + theme_void())
        }
        df_fixed <- data.frame(as.matrix(chosen_method[[input$k_fixed]]$partitions[[1]]$mb)) 
        df_fixed$cell <- rownames(df_fixed)
        clustering_choice <- choice()[2]
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
        
        ggplot2::ggplot(df_mat, ggplot2::aes(df_mat$X1, df_mat$X2)) + 
          ggplot2::geom_tile(ggplot2::aes(fill = value)) + 
          ggplot2::geom_text(ggplot2::aes(fill = df_mat$value, label = round(df_mat$value, 2))) +
          ggplot2::scale_fill_gradient2(low = muted("darkred"), 
                               mid = "white", 
                               high = muted("midnightblue"), 
                               midpoint = 0) + 
          ggplot2::scale_x_continuous(breaks = pretty(df_mat$X1, n = length(all_clusters_2))) +
          ggplot2::scale_y_continuous(breaks = pretty(df_mat$X2, n = length(all_clusters_1))) +
          ggplot2::theme(
            panel.background=ggplot2::element_rect(fill="white"), 
            axis.text.x = ggplot2::element_text(hjust = 1,vjust=1,size = 10,face = "bold"),
            axis.text.y = ggplot2::element_text(size = 10,face = "bold"),
            axis.title = ggplot2::element_text(size=14,face="bold"),
            axis.title.y = ggplot2::element_text(margin = margin(r = 20, l = 30)),
            axis.title.x = ggplot2::element_text(margin = margin(t = 20, b = 30))) + 
          ggplot2::xlab("Clusters in Configuration 2") +
          ggplot2::ylab("Clusters in Configuration 1") +
          ggplot2::labs(fill="JSI")
      })
      output$barcode_heatmap <- shiny::renderPlot({
        barcode_heatmap()
      })
      heatmap_filetype <- shiny::reactive({
        if (input$heatmap_filetype=='PDF'){
          filename <- 'heatmap.pdf'
          return(filename)
        }else if (input$heatmap_filetype=='PNG'){
          filename <- 'heatmap.png'
          return(filename)
        }else{
          filename <- 'heatmap.svg'
          return(filename)
        }
      })
      output$download_heatmap <- shiny::downloadHandler(
        filename = function() {heatmap_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,barcode_heatmap(),width = 30,
                 height = 20,
                 units = "cm")
        }
      )
    }
  )
}


