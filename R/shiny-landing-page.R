####### UI #######

#' UI - Landing page module
#'
#' @description Creates the UI interface for the landing page module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to identify the UI elements.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
ui_landing_page <- function(id) {
    ns <- shiny::NS(id)
    shiny::tabPanel("Home",
        icon = shiny::icon("home"),
        # titlePanel("ClustAssess"),
        shiny::headerPanel(
            list(
                shiny::tags$head(shiny::tags$style(".container {
                                        position: relative;
                                        text-align: center;
                                        color: white;
                                      }

                                      .bottom-left {
                                        position: absolute;
                                        bottom: 8px;
                                        left: 126px;
                                      }")),
                shiny::tags$script(
                    shiny::HTML("
                        var globalGenes = null;
                        Shiny.addCustomMessageHandler('genes', function(message) {
                            globalGenes = message;
                        });
                    ")
                ),
                shiny::HTML('<div class="container" style="width:100%; margin-top: 120px;">
                          <img src="https://raw.githubusercontent.com/Core-Bioinformatics/ClustAssess/main/man/figures/ClustAssess_starry_night.png" alt="starry" style="width:100%"/>
                          <!--<div class="bottom-left">Automated pipeline for assessing the robustness of single-cell clustering</div>-->
                       </div>')
            )
        ),
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
          </a>"
                )
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
          Authors: Andi Munteanu, Arash Shahsavari, Rafael Kollyfas, Miguel Larraz, Liviu Ciortuz and Irina Mohorianu"
                )
            )),
            shiny::column(6, shinyLP::panel_div(
                class_type = "info",
                panel_title = "Our website",
                content = shiny::HTML(
                    "Wellcome MRC - Cambridge Stem Cell Institute Core Bioinformatics Group:
          <br>
          <a href='https://www.corebioinf.stemcells.cam.ac.uk'>
              https://www.corebioinf.stemcells.cam.ac.uk
          </a>
          <br>
                 ClustAssess: <a hrer='https://www.corebioinf.stemcells.cam.ac.uk/pipelines-tools/tools/clustering-evaluation-clustassess'> https://www.corebioinf.stemcells.cam.ac.uk/pipelines-tools/tools/clustering-evaluation-clustassess </a>"
                )
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
          </a>"
                )
            ))
        )
    )
}

#' Server - Landing page module
#'
#' @description Creates the backend interface for the landing page module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to acess the UI elements.
#' @param height_ratio A reactive object that contains the height ratio of the
#' plots in the application (the height of the plot is calculated using the
#' height ratio and the height of the webpage).
#' @param dimension A reactive object that contains the dimensions of the
#' webpage.
#' @param parent_session The session of the parent module, used to control the
#' tabs of the application.
#' @param organism The organism of the dataset, which will be used in the
#' enrichment analysis.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
server_landing_page <- function(id, height_ratio, dimension, parent_session, organism = "hsapiens") {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            print(paste(Sys.time(), "landing - loading"))
            ggplot2::theme_set(ggplot2::theme_classic())
            add_env_variable("full_version", file.exists("stability.h5"))
            if (pkg_env$full_version) {
                add_env_variable("feature_ordering", rhdf5::h5read("stability.h5", "feature_ordering"))
            }

            if ("genes" %in% rhdf5::h5ls("expression.h5")$name) {
                genes <- rhdf5::h5read("expression.h5", "genes")
                index <- seq_along(genes)
                names(index) <- genes
                add_env_variable("genes", index)
            } else { # for backwards-compatibility purposes
                genes <- rhdf5::h5read("expression.h5", "genes_of_interest")
                index <- seq_along(genes)
                names(index) <- genes
                add_env_variable("genes_of_interest", index)
                genes <- rhdf5::h5read("expression.h5", "genes_others")
                index <- seq_along(genes)
                names(index) <- genes
                add_env_variable("genes_others", index)
            }

            shiny::observe({
                gene_dict <- as.list(names(genes))
                names(gene_dict) <- sapply(gene_dict, function(x) gene_name_transformation(x))
                session$sendCustomMessage("genes", gene_dict)
            })            

            rm(genes)
            rm(index)
            # gc()
            add_env_variable("cells", rhdf5::h5read("expression.h5", "cells"))
            if (pkg_env$full_version) {
                add_env_variable("feature_types", names(pkg_env$feature_ordering$original))
            }
            add_env_variable("height_ratio", height_ratio)

            add_env_variable("dimension", dimension)
            add_env_variable("plt_height", shiny::reactive(floor(pkg_env$height_ratio * pkg_env$dimension()[2])))
            add_env_variable("plt_width", shiny::reactive(pkg_env$dimension()[1]))
            add_env_variable("plt_width_half", shiny::reactive(pkg_env$dimension()[1] * 0.43))

            if (pkg_env$full_version) {
                discrete_colors <- rhdf5::h5read("stability.h5", "colors")
                bxplt_color <- rhdf5::h5read("stability.h5", "feature_stability/colours")
                discrete_colors[[as.character(length(bxplt_color))]] <- bxplt_color
            } else {
                discrete_colors <- list()
            }

            mtd <- readRDS("metadata.rds")
            for (color_options in mtd$metadata_colors) {
                ncolors <- as.character(length(color_options))
                discrete_colors[[ncolors]] <- color_options
            }

            if (utils::packageVersion("ClustAssess") <= "1.0.0") {
                for (discr_metadata in names(mtd$metadata_unique)) {
                    # transform NA values in factor
                    current_vals <- as.character(mtd$metadata[[discr_metadata]])
                    current_vals[is.na(current_vals)] <- "N/A"
                    current_vals <- factor(current_vals)

                    mtd$metadata[[discr_metadata]] <- current_vals
                    nlevs <- nlevels(current_vals)

                    if (nlevs == length(mtd$metadata_unique[[discr_metadata]])) {
                        next
                    }

                    mtd$metadata_unique[[discr_metadata]] <- levels(current_vals)

                    if (as.character(nlevs) %in% names(discrete_colors)) {
                        next
                    }

                    discrete_colors[[as.character(nlevs)]] <- generate_colours(nlevs, "pretty_dark")
                }
            }

            if (!pkg_env$full_version && length(intersect(c("UMAP_1", "UMAP_2"), colnames(mtd$metadata))) == 2) {
                umap_df <- mtd$metadata[, c("UMAP_1", "UMAP_2")]
                mtd$metadata <- mtd$metadata[, setdiff(colnames(mtd$metadata), c("UMAP_1", "UMAP_2"))]
                stab_obj <- list(umap = umap_df)
                add_env_variable("stab_obj", stab_obj)
            }

            add_env_variable("metadata", mtd$metadata)
            add_env_variable("metadata_temp", shiny::reactiveVal(mtd$metadata))
            # add_env_variable("metadata_colors", mtd$metadata_colors)
            add_env_variable("metadata_unique", mtd$metadata_unique)
            add_env_variable("metadata_unique_temp", shiny::reactiveVal(mtd$metadata_unique))

            add_env_variable("discrete_colors", discrete_colors)
            add_env_variable("organism", organism)

            add_env_variable("enable_markers_button", shiny::reactiveVal(-1))
            add_env_variable("find_markers_button", shiny::reactiveVal(-1))
            add_env_variable("annotation_button", shiny::reactiveVal(-1))
            add_env_variable("split_button", shiny::reactiveVal(-1))
            rm(mtd)

            shiny::observe({
                shiny::req(
                    pkg_env$dimension()[1] > 0,
                    pkg_env$dimension()[2] > 0
                )


                shiny::showTab("tabset_id", ifelse(pkg_env$full_version, "Dimensionality Reduction", "Comparison"), select = FALSE, session = parent_session)
                shiny::showTab("tabset_id", "Sandbox", select = FALSE, session = parent_session)
            }) %>% shiny::bindEvent(pkg_env$dimension(), once = TRUE)

            update_gears_width()
            print(paste(Sys.time(), "landing - finished loading"))
        }
    )
}
