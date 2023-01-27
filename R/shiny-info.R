
github <- shiny::a("https://github.com/Core-Bioinformatics/ClustAssess",
  href = "https://github.com/Core-Bioinformatics/ClustAssess",
  target = "_blank"
)

#### DIMENSIONAL REDUCTION ####

dr_title_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
      shiny::h3(shiny::strong("Assessing the stability of the dimensionality reduction")),
      shiny::br(),
      shiny::h5("Many single-cell experiments base their cell-type identification on the PhenoGraph pipeline by Levine et al. This pipeline can be sumarised into three main steps:"),
      shiny::h5("- dimensionality reduction"),
      shiny::h5("- graph construction"),
      shiny::h5("- graph clustering"),
      shiny::h5("Each of these steps contain a stochastic component that can alter the clustering results when the random seed is changed. The main goal of the ClustAssess package is to assess the stability of the clustering analysis and assist the user into choosing the most robust and reproducible configurations of parameters."),
      shiny::h5("\n"),
      shiny::h5("In the first step of the Dimensionality reduction, the main factor that affects the results is the feature set that is used for calculating the reduced space. The type of the feature size might be dependant upon the nature of the experiment."),
      shiny::h5("\n"),
      shiny::h5("In this page the user can interact and visualise the stability of different configurations of feature type - size of the feature set."),
      shiny::h1("\n"),
      shiny::h5("For more information please go to:"),
      shiny::tagList("", github),
      easyClose = TRUE
    ),
    session
  )
}

dr_individual_ecc_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
      shiny::h3(shiny::strong("Assessing the stability for each resolution parameter")),
      shiny::br(),
      shiny::h5("The following plot illustrates as a boxplot the distribution of the stability (measured using the ECC score) for different configurations."),
      shiny::h5("\n"),
      shiny::h5("Hovering over the boxplots will display more information about that specific distribution."),
      shiny::h5("Clicking a boxplot will trigger the display of the stability score distribution over an UMAP representation."),
      shiny::h5("\n"),
      shiny::h5("The results are displayed for each resolution parameter. When scrolling down the page the option to modify the resolution will appear. Changing the resolution value will affect the correspondent plots."),
      shiny::h1("\n"),
      shiny::h5("For more information please go to:"),
      shiny::tagList("", github),
      easyClose = TRUE
    ),
    session
  )
}

dr_individual_incremental_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
      shiny::h3(shiny::strong("Assessing the incremental stability for each resolution parameter")),
      shiny::br(),
      shiny::h5("One of the goals when choosing the feature set is to make sure we both have used a representative list that hasn't lead to introducing noise in the results."),
      shiny::h5("\n"),
      shiny::h5("The following plot compares the clustering obtained at two consecutive sizes of the set from the same feature type while varying the resolution parameter."),
      shiny::h5("\n"),
      shiny::h5("In our analysis, we consider a feature configuration to be stable if: (1) their consistency (ECC) across multiple seeds is high and (2) adding more features won't lead to significant changes in the clustering."),
      shiny::h1("\n"),
      shiny::h5("For more information please go to:"),
      shiny::tagList("", github),
      easyClose = TRUE
    ),
    session
  )
}

dr_overall_ecc_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
        shiny::h3(shiny::strong("Assessing the overall stability")),
        shiny::br(),
        shiny::h5("The previous plots help us infer the stability for each configuration while varying the resolution parameter."),
        shiny::h5("\n"),
        shiny::h5("The following plot summarizes the stability by picking the median value of the ECC score for each resolution value."),
        shiny::h1("\n"),
        shiny::h5("For more information please go to:"),
        shiny::tagList("", github),
        easyClose = TRUE
    ),
    session
  )
}

dr_overall_incremental_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
        shiny::h3(shiny::strong("Assessing the overall incremental stability")),
        shiny::br(),
        shiny::h5("The previous plots help us infer the incremental stability for each configuration while varying the resolution parameter."),
        shiny::h5("\n"),
        shiny::h5("The following plot summarizes the incremental stability by picking the median value of the ECC score for each resolution value."),
        shiny::h1("\n"),
        shiny::h5("For more information please go to:"),
        shiny::tagList("", github),
        easyClose = TRUE
    ),
    session
  )
}

dr_comparison_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
        shiny::h3(shiny::strong("Pairwise comparison between feature configurations")),
        shiny::br(),
        shiny::h5("The process of choosing the feature set can also be affected by the degree the metadata or the expression level distribution on the UMAP align with the biological interpretation."),
        shiny::h5("\n"),
        shiny::h5("Thus, we provide a section where we can select two different feature configurations and compare them by allowing the user to plot:"),
        shiny::h5("- the expression level of a gene on the UMAP produced by that feature set"),
        shiny::h5("- the cells that present a high expression on a list of genes above a threshold (voting scheme plot)"),
        shiny::h5("- the distribution of the metadata provided by the user on the UMAP"),
        shiny::h1("\n"),
        shiny::h5("For more information please go to:"),
        shiny::tagList("", github),
        easyClose = TRUE

    ),
    session
  )
}

dr_choice_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
        shiny::h3("Fixing the feature configuration"),
        shiny::br(),
        shiny::h5("Before going to the next steps of the pipeline, the user should decide upon one feature type and one size of the feature that will be fixed in the downstream analysis"),
        shiny::h5("\n"),
        shiny::h5("For each feature type, we will provide a recommendation in terms of the most stable number of features."),
        shiny::h1("\n"),
        shiny::h5("For more information please go to:"),
        shiny::tagList("", github),
        easyClose = TRUE
    ),
    session
  )
}


dr_individual_incremental_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
        shiny::h1("\n"),
        shiny::h5("For more information please go to:"),
        shiny::tagList("", github),
        easyClose = TRUE

    ),
    session
  )
}

dr_individual_incremental_info <- function(session) {
  shiny::showModal(
    shiny::modalDialog(
        shiny::h1("\n"),
        shiny::h5("For more information please go to:"),
        shiny::tagList("", github),
        easyClose = TRUE

    ),
    session
  )
}
# dimensionality_reduction_info <- function(ns) {
#     shinyWidgets::dropdownButton(
#         inputId = ns("info"),
#         label = "",
#             icon = shiny::icon("circle-info"),
#             status = "success",
#             size = "sm",
#         shiny::h3(shiny::strong("Assessing the stability of the dimensionality reduction")),
#         shiny::br(),
#         shiny::h5("lalala")
#         # placement = "top-end"
#     )
# }

# dimensionality_reduction_info <- function(ns) {
#     shiny::actionButton(
#         inputId = ns("info"),
#         label = "",
#         icon = shiny::icon("info"),
#         class = "btn-success",
#         width = "20px"
#     )
# }
