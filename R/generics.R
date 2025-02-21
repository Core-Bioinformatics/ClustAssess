#' Create the ClustAssess ShinyApp
#'
#' @description Creates the ClustAssess ShinyApp based on the output of the
#' automatic ClustAssess pipeline. In addition to that, the expression matrix
#' and the metadata dataframe are provided as input to the ShinyApp.

#' @param object A Seurat object or an expression matrix
#' @param clustassess_object The output of the ClustAssess automatic pipeline
#' @param project_folder The folder where the files will be written
#' @param metadata The metadata dataframe. This parameter will be ignored if
#' the object is a Seurat object.
#' @param assay_name The name of the assay to be used to extract the expression matrix
#' from the Seurat object. This parameter will be ignored if the object is not
#' a Seurat object.
#' @param compression_level The compression level for the h5 files (See `rhdf5::h5createFile`` for more details)
#' @param summary_function The function used for summarizing the stability values; the default is `median`
#' @param shiny_app_title The title of the shiny app
#' @param organism_enrichment The organism used for the enrichment analysis; the default is `hsapiens`
#' @param height_ratio The ratio of the height of the plot to the height of the browser; the default is `0.6`
#' @param qualpalr_colorspace The colorspace used for generating the colors; the default is `pretty`
#' @param prompt_feature_choice Should the user be prompted to choose if he wants to continue with the selection of features even if it is lower than median sequence depth; the default is `TRUE`
#' @rdname write_shiny_app
#' @export
write_shiny_app <- function(object,
                            metadata = NULL,
                            assay_name = NULL,
                            clustassess_object,
                            project_folder,
                            compression_level = 6,
                            summary_function = stats::median,
                            shiny_app_title = "",
                            organism_enrichment = "hsapiens",
                            height_ratio = 0.6,
                            qualpalr_colorspace = "pretty",
                            prompt_feature_choice = TRUE) {
    UseMethod(generic = "write_shiny_app", object = object)
}
