
ui_graph_construction <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    "Graph Construction",
    h1("Relationship between the number of neighbours and the number of connected components"),
    shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info"), status = 'success', icon = shiny::icon('info'),size='sm'),
                           shiny::h3(shiny::strong('Relationship between the number of neighbours and the number of connected components')),
             shiny::br(),
             shiny::h5('This plot describes the covariation between the number neighbours and the number of connected components obtained using both PCA and UMAP reductions as base for graph building. As the number of neighbours increases, the number of connected components decreases (this is an expected result, as increasing the number of neighbours result in a better connected graph). Please note that the number of connected components provides a lower bound on the number of clusters we can obtain by downstream community detection algorithms such as Louvain and Leiden.'),
             shiny::h1('\n'),
             shiny::h5('For more information please go to:'),
             shiny::tagList("", github),
             placement = "right",
             arrow = F,
             maxWidth = '700px'),
    shinyWidgets::dropdownButton(
      label = "",
      icon = shiny::icon("download"),
      status = "success",
      size='sm',
      shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
      shiny::textInput(ns("filename_con_comps"), "File name:", width = "80%"),
      shiny::numericInput(ns("width_con_comps"), "Width (in):", 7, 3, 100, 0.1),
      shiny::numericInput(ns("height_con_comps"), "Height (in):", 7, 3, 100, 0.1),
      shiny::selectInput(ns('filetype_con_comps'),'Filetype',choices = c('PDF','PNG','SVG'),width='50%',selected='PDF'),
      shiny::downloadButton(ns('download_con_comps'), label="Download Plot")
      
    ),
    shiny::splitLayout(cellWidths = c("15%", "85%"),shiny::div(style="width:90%;",shiny::verticalLayout(shiny::selectInput(ns("palette_plot_conn_comps"), "Colour scheme:",
                                                                                               c("Colour scheme 1" = "col_1",
                                                                                                 "Colour scheme 2" = "col_2",
                                                                                                 "Colour scheme 3" = "col_3"),size=3,selectize=F))),
                       shiny::plotOutput(ns("plot_conn_comps"))),
    
    h1("Relationship between the number of neighbours and then number of clusters"),
    shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info_2"), status = 'success', icon = shiny::icon('info'),size='sm'),
                           shiny::h3(shiny::strong('Relationship between the number of neighbours and the number of connected components')),
                           shiny::br(),
                           shiny::h5('This plot summarizes the relationship between the number of nearest neighbours and k, the number of clusters.he plot also illustrates the difference between the two graph types: SNN has a tighter distribution of k (over multiple iterations) compared to NN. If the initial object contains graphs based on both UMAP and PCA embedding, this plot will showcase the impact of this choice, as well.'),
                           shiny::h5('This plot is interactive:'),
                           shiny::h5('- You can select a desired set of configurations, and chose between three different colour schemes.'),
                           shiny::h5('- You can hover over the plot to show detailed information on the different boxes.'),
                           shiny::h5('- You can click anywhere in the plot to obtain a UMAP representation coloured with the stability (ECC) across different seeds for a specific number of nearest neighbours.'),
                           shiny::h1('\n'),
                           shiny::h5('For more information please go to:'),
             shiny::tagList("", github),
             placement = "right",
             arrow = F,
             maxWidth = '700px'),
    shinyWidgets::dropdownButton(
      label = "",
      icon = shiny::icon("download"),
      status = "success",
      size='sm',
      shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
      shiny::textInput(ns("filename_neigh_k"), "File name:", width = "80%"),
      shiny::numericInput(ns("width_neigh_k"), "Width (in):", 7, 3, 100, 0.1),
      shiny::numericInput(ns("height_neigh_k"), "Height (in):", 7, 3, 100, 0.1),
      shiny::selectInput(ns('filetype_neigh_k'),'Filetype',choices = c('PDF','PNG','SVG'),width='50%',selected='PDF'),
      shiny::downloadButton(ns('download_neigh_k'), label="Download Plot")
      
    ),
    shiny::fluidRow(shiny::tags$head(shiny::tags$style(shiny::HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
                    shiny::splitLayout(cellWidths = c("15%", "55%","30%"),
                                       shiny::div(style="width:90%;",shiny::verticalLayout(shiny::checkboxGroupInput(ns('sel_conn_comps'),'Select configurations',choices=config_choices,selected=config_choices, width='50%'),
                                                                                           shiny::selectInput(ns("palette_plot_neigh_k"), "Colour scheme:",
                                                                           c("Colour scheme 1" = "col_1",
                                                                             "Colour scheme 2" = "col_2",
                                                                             "Colour scheme 3" = "col_3"),size=3,selectize=F),
                                                                           shiny::verbatimTextOutput(ns("hover_info")),
                                                                           shiny::h1('\n'))),
                                       shiny::plotOutput(ns("plot_neigh_k"),hover = shiny::hoverOpts(id = ns("nn_hover"), delayType = "throttle"),click = shiny::clickOpts(id=ns("ecc_click")), width="100%", height="700px"),
                                       shiny::div(style="width:90%;",shiny::verticalLayout(shiny::verbatimTextOutput(ns("click_info")),shiny::plotOutput(ns('umap_ecc'))))
             )
    ),
    shiny::h1("Stability"),
    shinyWidgets::dropMenu(shinyWidgets::dropdownButton(ns("Info_3"), status = 'success', icon = shiny::icon('info'),size='sm'),
                    shiny::h3(shiny::strong('Stability of different configurations across a varying number of nearest neighbours')),
                    shiny::br(),
                    shiny::h5('The stability of these parameters can be also evaluated using the Element-Centric Consistency applied on the partition list obtained over the multiple runs. The following summary plot underlines the evolution of the consistency as the number of neighbours increases.'),
                    shiny::h5('This plot is interactive:'),
                    shiny::h5('- You can select a desired set of configurations, and chose between three different colour schemes.'),
                    shiny::h5('- You can hover over the plot to show detailed information on the different boxes.'),
                    shiny::h5('- You can click anywhere in the plot to obtain a UMAP representation coloured with the stability (ECC) across different seeds for a specific number of nearest neighbours.'),
                    shiny::h1('\n'),
                    shiny::h5('For more information please go to:'),
                    shiny::tagList("", github),
             placement = "right",
             arrow = F,
             maxWidth = '700px'),
    shinyWidgets::dropdownButton(
      label = "",
      icon = shiny::icon("download"),
      status = "success",
      size='sm',
      shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
      shiny::textInput(ns("filename_neigh_stability"), "File name:", width = "80%"),
      shiny::numericInput(ns("width_neigh_stability"), "Width (in):", 7, 3, 100, 0.1),
      shiny::numericInput(ns("height_neigh_stability"), "Height (in):", 7, 3, 100, 0.1),
      shiny::selectInput(ns('filetype_neigh_stability'),'Filetype',choices = c('PDF','PNG','SVG'),width='50%',selected='PDF'),
      shiny::downloadButton(ns('download_neigh_stability'), label="Download Plot")
      
    ),
    shiny::fluidRow(
      shiny::splitLayout(cellWidths = c("15%", "55%","30%"),
                         shiny::verticalLayout(shiny::div(style="width:90%;",shiny::checkboxGroupInput(ns('sel_stab'),'Select configurations',choices=config_choices,selected=config_choices)),
                                               shiny::selectInput(ns("palette_neigh_stability"), "Colour scheme:",
                                             c("Colour scheme 1" = "col_1",
                                               "Colour scheme 2" = "col_2",
                                               "Colour scheme 3" = "col_3"),size=3,selectize=F),
                                             shiny::verbatimTextOutput(ns("hover_info_nnstab")),
                                             shiny::h1('\n')),
                         shiny::plotOutput(ns("plot_neigh_stability"),hover = shiny::hoverOpts(id = ns("nn_stability_hover"), delayType = "throttle"),click = shiny::clickOpts(id=ns("neigh_stability_click")), width="100%", height="700px"),
                         shiny::div(style="width:90%;",shiny::verticalLayout(shiny::verbatimTextOutput(ns("click_info_2")),shiny::plotOutput(ns('neigh_umap_ecc')))))
    ),style="margin-left:10px "
  )
}


##### SERVER #####

server_graph_construction <- function(id,chosen_config){
  shiny::moduleServer(
    id,
    function(input, output, session) {
      #Reactive for first plot
      plot_conn_comps_Input <- shiny::reactive({
        if (input$palette_plot_conn_comps=='col_1'){
          option <- 'RdBu'
        }else if (input$palette_plot_conn_comps=='col_2'){
          option <- 'RdYIGn'
        }else{
          option <-'Greys'
        }
        ClustAssess::plot_connected_comps_evolution(chosen_config$nn_conn_comps) +
          ggplot2::scale_fill_brewer(palette = option)
      })
      #1 Plot Relationship between the number of neighbours and then number of connected components
      output$plot_conn_comps <- shiny::renderPlot({
        plot_conn_comps_Input()
      }, height = 400, width = 1000)
      #generate filenames
      con_comps_filetype <- shiny::reactive({
        if (input$filetype_con_comps=='PDF'){
          filename <- paste0(input$filename_con_comps,'.pdf')
          return(filename)
        }else if (input$filetype_con_comps=='PNG'){
          filename <- paste0(input$filename_con_comps,'.png')
          return(filename)
        }else{
          filename <- paste0(input$filename_con_comps,'.svg')
          return(filename)
        }
      })
      #Download Plots
      output$download_con_comps <- shiny::downloadHandler(
        filename = function() {con_comps_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,plot_conn_comps_Input(),width = input$width_con_comps,
                 height = input$height_con_comps,
                 units = "in",
                 limitsize = FALSE)
        }
      )
      #Add hovering output for second plot, and add the text box
      ecc_mouse_hover <- shiny::reactiveVal(c(NA, NA,NA))
      output$hover_info <- shiny::renderPrint({
        if(is.null(input$nn_hover)) {
          return(cat('Hover over the plot to see more'))
        }
        x <- input$nn_hover$x
        col <- round(x)
        obj <- chosen_config$nn_importance$n_neigh_k_corresp[input$sel_conn_comps]
        obj <- reshape2::melt(obj)
        val <- sort(unique(obj$L2))[col]
        configs <- unique(obj$L1)
        cat(paste0('Selected NN:','\t',val,'\n'))
        for (config in configs){
          filt_obj <- obj[obj$L1==config,]
          
          maxi <- max(filt_obj$value[filt_obj$L2==val])
          mini <- min(filt_obj$value[filt_obj$L2==val])
          iqr <- stats::IQR(filt_obj$value[filt_obj$L2==val])
          
          cat(paste0(config,'\n'))
          cat(paste0('Max: ',maxi,'\n'))
          cat(paste0('Min: ',mini,'\n'))
          cat(paste0('IQR: ',iqr,'\n'))
        }
      })
      #Reactive for second plot
      plot_n_neigh_k_correspondence_Input <- shiny::reactive({
        if (input$palette_plot_neigh_k=='col_1'){
          option <- 'RdBu'
        }else if (input$palette_plot_neigh_k=='col_2'){
          option <- 'RdYIGn'
        }else{
          option <-'Greys'
        }
        if (length(input$sel_conn_comps)==0){
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }else{
          ClustAssess::plot_n_neigh_k_correspondence(chosen_config$nn_importance) +
            ggplot2::scale_fill_brewer(palette = option)
        }
      })
      #2 Plot Relationship between the number of neighbours and then number of clusters
      output$plot_neigh_k <- shiny::renderPlot({
        plot_n_neigh_k_correspondence_Input()
      })
      neigh_k_filetype <- shiny::reactive({
        if (input$filetype_neigh_k=='PDF'){
          filename <- paste0(input$filename_neigh_k,'.pdf')
          return(filename)
        }else if (input$filetype_neigh_k=='PNG'){
          filename <- paste0(input$filename_neigh_k,'.png')
          return(filename)
        }else{
          filename <- paste0(input$filename_neigh_k,'.svg')
          return(filename)
        }
      })
      #Download Plots
      output$download_neigh_k <- shiny::downloadHandler(
        filename = function() {neigh_k_filetype()},
        content = function(file) {
          ggplot2::ggsave(file,plot_n_neigh_k_correspondence_Input(),width = input$width_neigh_k,
                 height = input$height_neigh_k,
                 units = "in")
        }
      )
      #Modify the function for the third plot
      plot_n_neigh_ecs_mod <- function(nn_ecs_object,
                                       boxplot_width = 0.5,palette_option) {
        obj <- nn_ecs_object$n_neigh_ec_consistency[input$sel_stab]
        melted_obj <- reshape2::melt(obj)
        colnames(melted_obj) <- c("ECC", "n_neigh", "config_name")
        
        melted_obj$n_neigh <- factor(melted_obj$n_neigh)
        melted_obj$n_neigh <- factor(melted_obj$n_neigh, levels(melted_obj$n_neigh)[stringr::str_order(levels(melted_obj$n_neigh), numeric = TRUE)])
        
        ggplot2::ggplot(
          melted_obj,
          ggplot2::aes(
            x = .data$n_neigh,
            y = .data$ECC,
            fill = .data$config_name
          )
        ) +
          #ggplot2::geom_boxplot(width = boxplot_width) +
          ggplot2::geom_boxplot() +
          ggplot2::theme_classic() +
          ggplot2::scale_fill_brewer(palette = palette_option) +
          ggplot2::labs(
            x = "# of nearest neighbors",
            y = "EC consistency",
            fill = "configuration"
          ) +
          ggplot2::ggtitle("Distribution of ECC across different seeds for different # neighbors")
      }
      #Plot the UMAP + add some text on the chosen selection
      output$click_info<- shiny::renderPrint({
        if(is.null(input$ecc_click)) {
          return(cat('Click the plot to see more'))
        }
        col <- round(input$ecc_click$x)
        obj <- chosen_config$nn_importance$n_neigh_k_corresp[input$sel_conn_comps]
        obj <- reshape2::melt(obj)
        val <- sort(unique(obj$L2))[col]
        configs <- unique(obj$L1)
        if (input$ecc_click$x<0.5){
          selection <- sort(config_choices)[1]
        }else{
          steps <- length(configs)
          input_steps <- (input$ecc_click$x - (0.5+(col-1)))*steps
          selected_config <- ceiling(input_steps)
          selection <- sort(configs)[selected_config]
        }
        cat(paste0('UMAP:','\n',selection,'\n',val,' NNs'))
      })
      output$umap_ecc <- shiny::renderPlot({
        if(is.null(input$ecc_click)) {
          return(ggplot() + theme_void())
        }
        #get click info for ecc
        col <- round(input$ecc_click$x)
        obj <- chosen_config$nn_importance$n_neigh_k_corresp[input$sel_conn_comps]
        obj <- reshape2::melt(obj)
        val <- sort(unique(obj$L2))[col]
        configs <- unique(obj$L1)
        if (input$ecc_click$x<0.5){
          selection <- sort(config_choices)[1]
        }else{
          steps <- length(configs)
          input_steps <- (input$ecc_click$x - (0.5+(col-1)))*steps
          selected_config <- ceiling(input_steps)
          selection <- sort(config_choices)[selected_config]
        }
        ecc <- chosen_config$nn_importance$n_neigh_ec_consistency[input$sel_stab]
        ecc <- ecc[[selection]][[col]]
        
        ggplot2::ggplot(data.frame(
          chosen_config$umap),
          aes(x = .data$UMAP_1,
              y = .data$UMAP_2,
              color = ecc)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_viridis_c() +
          ggplot2::theme_bw()
      })
      #Reactive for third plot
      plot_neigh_stability_Input <- shiny::reactive({
        if (input$palette_neigh_stability=='col_1'){
          option <- 'RdBu'
        }else if (input$palette_neigh_stability=='col_2'){
          option <- 'RdYIGn'
        }else{
          option <-'Greys'
        }
        if (length(input$sel_conn_comps)==0){
          return(ggplot2:ggplot() + ggplot2::theme_void())
        }else{
          plot_n_neigh_ecs(chosen_config$nn_importance) +
            ggplot2::scale_fill_brewer(palette = option)
        }
      })
      #3. Plot Stability
      output$plot_neigh_stability <- shiny::renderPlot({
        plot_neigh_stability_Input()
      })
      #generate filenames
      neigh_stability_filetype <- shiny::reactive({
        if (input$filetype_neigh_stability=='PDF'){
          filename <- paste0(input$filename_neigh_stability,'.pdf')
          return(filename)
        }else if (input$filetype_neigh_stability=='PNG'){
          filename <- filename <- paste0(input$filename_neigh_stability,'.png')
          return(filename)
        }else{
          filename <- filename <- paste0(input$filename_neigh_stability,'.svg')
          return(filename)
        }
      })
      #Download Plots
      output$download_neigh_stability <- shiny::downloadHandler(
        filename = function() {neigh_stability_filetype()},
        content = function(file) {
          ggsave(file,plot_neigh_stability_Input(),width = input$width_neigh_stability,
                 height = input$height_neigh_stability,
                 units = "in")
        }
      )
      #Hover info
      output$hover_info_nnstab <- shiny::renderPrint({
        if(is.null(input$nn_stability_hover)) {
          return(cat('Hover over the plot to see more'))
        }
        x <- input$nn_stability_hover$x
        col <- round(x)
        obj <- chosen_config$nn_importance$n_neigh_ec_consistency[input$sel_stab]
        obj <- reshape2::melt(obj)
        val <- sort(unique(obj$L2))[col]
        configs <- unique(obj$L1)
        cat(paste0('Selected NN:','\t',val,'\n'))
        for (config in configs){
          filt_obj <- obj[obj$L1==config,]
          
          maxi <- round(max(filt_obj$value[filt_obj$L2==val]),digits=2)
          mini <- round(min(filt_obj$value[filt_obj$L2==val]),digits=2)
          iqr <- round(IQR(filt_obj$value[filt_obj$L2==val]),digits=2)
          
          cat(paste0(config,'\n'))
          cat(paste0('Max: ',maxi,'\n'))
          cat(paste0('Min: ',mini,'\n'))
          cat(paste0('IQR: ',iqr,'\n'))
        }
      })
      #Plot UMAP + add some text on the chosen selection
      output$click_info_2<- shiny::renderPrint({
        if(is.null(input$neigh_stability_click)) {
          return(cat('Click the plot to see more'))
        }
        col <- round(input$neigh_stability_click$x)
        obj <- chosen_config$nn_importance$n_neigh_ec_consistency[input$sel_stab]
        obj <- reshape2::melt(obj)
        val <- sort(unique(obj$L2))[col]
        configs <- unique(obj$L1)
        if (input$neigh_stability_click$x<0.5){
          selection <- sort(config_choices)[1]
        }else{
          steps <- length(configs)
          input_steps <- (input$neigh_stability_click$x - (0.5+(col-1)))*steps
          selected_config <- ceiling(input_steps)
          selection <- sort(configs)[selected_config]
        }
        cat(paste0('UMAP:','\n',selection,'\n',val,' NNs'))
      })
      output$neigh_umap_ecc <- shiny::renderPlot({
        if(is.null(input$neigh_stability_click)) {
          return(ggplot() + theme_void())
        }
        #get click info for ecc
        col <- round(input$neigh_stability_click$x)
        obj <- chosen_config$nn_importance$n_neigh_ec_consistency[input$sel_stab]
        obj <- reshape2::melt(obj)
        val <- sort(unique(obj$L2))[col]
        configs <- unique(obj$L1)
        if (input$neigh_stability_click$x<0.5){
          selection <- sort(config_choices)[1]
        }else{
          steps <- length(configs)
          input_steps <- (input$neigh_stability_click$x - (0.5+(col-1)))*steps
          selected_config <- ceiling(input_steps)
          selection <- sort(configs)[selected_config]
        }
        ecc <- chosen_config$nn_importance$n_neigh_ec_consistency[input$sel_stab]
        ecc <- ecc[[selection]][[col]]
        
        ggplot2::ggplot(data.frame(
          chosen_config$umap),
          aes(x = .data$UMAP_1,
              y = .data$UMAP_2,
              color = ecc)) +
          ggplot2::geom_point() +
          ggplot2::scale_color_viridis_c() +
          ggplot2::theme_bw()
      })
    }
  )
}
