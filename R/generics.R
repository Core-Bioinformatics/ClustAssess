#' Writing the shiny app folder
#'
#' @description  to be completed
#'
#' @rdname write_shiny_app
#' @export
write_shiny_app <- function(object, ...) {
    UseMethod(generic = "write_shiny_app", object = object)
}
