####### UI #######

#' Writing objects
#'
#' @description to be completed
#'
#' @export
ui_landing_page <- function(id){
  ns <- shiny::NS(id)
  shiny::tabPanel("Home",
                  icon = shiny::icon("home"),
                  titlePanel("ClustAssess"),
                  shiny::headerPanel(
                    list(shiny::tags$head(shiny::tags$style(".container {
                                        position: relative;
                                        text-align: center;
                                        color: white;
                                      }
                                      
                                      .bottom-left {
                                        position: absolute;
                                        bottom: 8px;
                                        left: 126px;
                                      }")),
                         shiny::HTML('<div class="container" style="width:100%">
                          <img src="https://raw.githubusercontent.com/Core-Bioinformatics/ClustAssess/release-0.4.0/docs/articles/ClustAssess_files/figure-html/ClustAssess_starry_night.png" alt="starry" style="width:100%"/>
                          <!--<div class="bottom-left">Automated pipeline for assessing the robustness of single-cell clustering</div>-->
                       </div>')
                    )),
                  shiny::fluidRow(
                    shiny::column(6, shinyLP::panel_div(
                      class_type = "info", 
                      panel_title = "Tutorial",
                      content = shiny::HTML(
                        "Vignette describing the different panels in this app and how to use them:
          <br>
          <a href='https://core-bioinformatics.github.io/ClustAssess/articles/stability-based-parameter-assessment.html'>
              https://core-bioinformatics.github.io/ClustAssess/articles/stability-based-parameter-assessment.html
          </a>"
                      )
                    )),
                    shiny::column(6, shinyLP::panel_div(
                      class_type = "primary", 
                      panel_title = "GitHub",
                      content = shiny::HTML(
                        "Latest stable version, bug reports and feature requests:
          <br>
          <a href='https://github.com/Core-Bioinformatics/ClustAssess'>
              https://github.com/Core-Bioinformatics/ClustAssess
          </a>")
                    ))
                  ),
                  shiny::fluidRow(
                    shiny::column(6, shinyLP::panel_div(
                      class_type = "primary", 
                      panel_title = "Manuscript",
                      content = shiny::HTML(
                        "ClustAssess: automated pipeline for assessing the robustness of single-cell clustering
          <br>
          <br>
          <a href='https://www.biorxiv.org/content/10.1101/2022.01.31.478592v1.full'>
              https://www.biorxiv.org/content/10.1101/2022.01.31.478592v1.full
          </a>
          <br>
          Authors: Arash Shahsavari, Andi Munteanu, Miguel Larraz, Liviu Ciortuz and Irina Mohorianu"
                      )
                    )),
                    shiny::column(6, shinyLP::panel_div(
                      class_type = "info", 
                      panel_title = "Our website",
                      content = HTML(
                        "Wellcome MRC - Cambridge Stem Cell Institute Core Bioinformatics Group:
          <br>
          <a href='https://www.corebioinf.stemcells.cam.ac.uk'>
              https://www.corebioinf.stemcells.cam.ac.uk
          </a>
          <br>
                 ClustAssess: <a hrer='https://www.corebioinf.stemcells.cam.ac.uk/pipelines-tools/tools/clustering-evaluation-clustassess'> https://www.corebioinf.stemcells.cam.ac.uk/pipelines-tools/tools/clustering-evaluation-clustassess </a>")
                    ))
                  ),
                  shiny::fluidRow(
                    shiny::column(6, shinyLP::panel_div(
                      class_type = "info", 
                      panel_title = "Created by",
                      content = shiny::HTML(
                        "Arash Shahsavari, Andi Munteanu, Miguel Larraz, Liviu Ciortuz and Irina Mohorianu
                 <br>"
                      )
                    )),
                    shiny::column(6, shinyLP::panel_div(
                      class_type = "primary", 
                      panel_title = "CRAN",
                      content = shiny::HTML(
                        "For more information, please visit our CRAN page:
          <br>
          <a href='https://cran.r-project.org/package=ClustAssess'>
              https://cran.r-project.org/package=ClustAssess
          </a>")
                    ))
                  )
                  
  )
  
}


#' Landing page - server side
#'
#' @description to be completed
#'
#' @export
server_landing_page <- function(id, height_ratio, dimension, parent_session) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
        print(paste(Sys.time(), "landing - loading"))
        ggplot2::theme_set(ggplot2::theme_classic())
        add_env_variable("feature_ordering", rhdf5::h5read("stability.h5", "feature_ordering"))
      
        genes <- rhdf5::h5read("expression.h5", "genes_of_interest")
        index <- seq_along(genes)
        names(index) <- genes
        add_env_variable("genes_of_interest", index)
        genes <- rhdf5::h5read("expression.h5", "genes_others")
        index <- seq_along(genes)
        names(index) <- genes
        add_env_variable("genes_others", index)
        rm(genes)
        rm(index)
        # gc()
        add_env_variable("cells", rhdf5::h5read("expression.h5", "cells"))
        add_env_variable("feature_types", names(pkg_env$feature_ordering$original))
        add_env_variable("height_ratio", height_ratio)

        add_env_variable("dimension", dimension)
        mdt <- readRDS("metadata.rds")

        add_env_variable("metadata", mdt$metadata)
        add_env_variable("metadata_colors", mdt$metadata_colors)
        add_env_variable("metadata_unique", mdt$metadata_unique)
        add_env_variable("enable_markers_button", shiny::reactiveVal(-1))
        add_env_variable("find_markers_button", shiny::reactiveVal(-1))
        rm(mdt)

        shiny::observe({
          shiny::req(pkg_env$dimension()[1] > 0,
                      pkg_env$dimension()[2] > 0)
                      
          
          shiny::showTab("tabset_id", "Dimensionality Reduction", select = FALSE, session = parent_session)
          shiny::showTab("tabset_id", "Sandbox", select = FALSE, session = parent_session)
        }) %>% shiny::bindEvent(pkg_env$dimension(), once = TRUE)
        print(paste(Sys.time(), "landing - finished loading"))
        gc()
    }
    
  )
}
