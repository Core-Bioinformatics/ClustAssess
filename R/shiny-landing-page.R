ui_landing_page <- function(id){
  ns <- shiny::NS(id)
  shiny::tabPanel("Home",
           icon = shiny::icon("home"),
           # dashboardBody(
           #   tags$img(
           #     src = "starry_night_blur.png",
           #     style = 'height:450px;width:100%;position: absolute'
           #   )),
           shiny::headerPanel(
             list(tags$head(tags$style(".container {
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
                          <img src="starry_night_4.png" alt="starry" style="width:100%"/>
                          <!--<div class="bottom-left">Automated pipeline for assessing the robustness of single-cell clustering</div>-->
                       </div>')
             )),
           shiny::fluidRow(
             shiny::column(6, shinyLP::panel_div(
               class_type = "danger", 
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
               class_type = "success", 
               panel_title = "Manuscript",
               content = shiny::HTML(
                 "For a more detailed overview of the Package and app functionality:
                 <br>
                 ClustAssess: automated pipeline for assessing the robustness of single-cell clustering
          <br>
          <a href='https://www.biorxiv.org/content/10.1101/2022.01.31.478592v1.full'>
              https://www.biorxiv.org/content/10.1101/2022.01.31.478592v1.full
          </a>"
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
          </a>")
             ))
           ),
           shiny::fluidRow(
             shiny::column(6, shinyLP::panel_div(
               class_type = "warning", 
               panel_title = "Created by",
               content = shiny::HTML(
                 "Arash Shahsavari, Andi Munteanu, Miguel Larraz, Liviu Ciortuz and Irina Mohorianu"
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