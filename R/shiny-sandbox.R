####### UI #######
ui_colour_picker <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h1("Change discrete colour scheme"),
        shiny::splitLayout(
            cellWidths = c("20%", "80%"),
            colourpicker::colourInput(
                inputId = ns("colour_picker"),
                label = "Pick a colour to get the hexcode",
                value = "red",
                allowTransparent = TRUE,
                closeOnClick = TRUE
            ),
            shiny::tagList(
                shiny::h3("Qualpalr palette generator"),
                shiny::splitLayout(
                    cellWidths = c("120px", "170px", "220px", "220px", "220px"),
                    shiny::numericInput(
                        inputId = ns("n_colours"),
                        label = "Number of colours",
                        value = 10,
                        min = 1,
                        max = 99,
                        width = "100px",
                    ),
                    shiny::selectInput(
                        inputId = ns("preset"),
                        label = "Preset",
                        choices = c("pretty", "pretty_dark", "pastels", "pastels_dark"),
                        width = "150px"
                    ),
                    shiny::sliderInput(
                        inputId = ns("hue"),
                        label = "Hue",
                        min = -360, max = 360,
                        value = c(0, 360),
                        width = "200px"
                    ),
                    shiny::sliderInput(
                        inputId = ns("saturation"),
                        label = "Saturation",
                        min = 0, max = 1,
                        value = c(0.1, 0.3),
                        width = "200px"
                    ),
                    shiny::sliderInput(
                        inputId = ns("lightness"),
                        label = "Lightness",
                        min = 0, max = 1,
                        value = c(0.3, 0.6),
                        width = "200px"
                    )
                ),
                shiny::splitLayout(
                    cellWidths = c("15%", "85%"),
                        shiny::actionButton(
                            inputId = ns("generate_palette"),
                            label = "Generate palette",
                            style = "font-size:20px;",
                            class = "btn-success"
                        ),
                        shiny::verbatimTextOutput(ns("palette_output"))
                ),
                shiny::plotOutput(ns("palette_plot"))
            )
        ),
        shiny::splitLayout(
            shiny::tagList(
                shiny::h3("Update colour palette for a given number of groups"),
                shiny::splitLayout(
                    cellWidths = c("120px", "600px", "290px"),
                    shiny::numericInput(
                        inputId = ns("n_groups"),
                        label = "Number of groups",
                        value = 10,
                        min = 1,
                        max = 99,
                        width = "100px"
                    ),
                    shiny::textAreaInput(
                        inputId = ns("colours_input"),
                        label = "Enter hexcodes of colours, separated by commas",
                        value = "",
                        width = "100%"
                    )
                ),
                shiny::actionButton(
                    inputId = ns("update_colours"),
                    label = "Update colours for ngroups",
                    style = "font-size:20px;",
                    class = "btn-danger",
                    width = "270px"
                )
            ),
            shiny::tagList(
                shiny::h3("Save the current colour configurations"),
                shiny::downloadButton(ns("download_colours"), "Download current colours (JSON format)"),
                shiny::h3("Update colour configuration from an uploaded file"),
                shiny::p("The files should follow the JSON format, where the keys are the number of groups and the values are the hexcodes of the colours."),
                shiny::fileInput(ns("upload_colours"), "Upload a file with colour configurations (JSON format)", accept = ".json")
            )
        )
    )

}

ui_sandbox_config_choice <- function(id, draw_line) {
    ns <- shiny::NS(id)
    style <- ifelse(draw_line, "border-right:5px solid", "")

    shinyWidgets::panel(
        style = style,
        shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info_Sandoboxs"), status = "success", icon = shiny::icon("info"), size = "sm"),
            shiny::h3(shiny::strong("Compare any two configurations")),
            shiny::div(
                style = "white-space: pre-wrap; /* css-3 */
                                                                      white-space: -moz-pre-wrap; /* Mozilla, since 1999 */
                                                                      white-space: -pre-wrap; /* Opera 4-6 */
                                                                      white-space: -o-pre-wrap; /* Opera 7 */
                                                                      word-wrap: break-word; /* Internet Explorer 5.5+ */",
                shiny::h5("In this section you can compare any configuration to another. Here, you can also colour each configuration by the ECC, clusters, as well as any other metadata features that you have previously specified. You should choose the number of clusters that makes the most sense to you. Please note that this tab should not be used to select a final configuration. We encourage you to go throught the inital tabs to find the most suitable configuration for your study.")
            ),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href = "https://github.com/Core-Bioinformatics/ClustAssess", target = "_blank")),
            placement = "right",
            arrow = F,
            maxWidth = "700px"
        ),
        shiny::div(style = "width:90%;", shiny::verticalLayout(
            shiny::h1("Select a configuration"),
            shiny::splitLayout(
                cellWidths = c("50%", "50%"),
                shiny::verticalLayout(
                    shiny::uiOutput(ns("sandbox_sel_fset_render")),
                    shiny::uiOutput(ns("sandbox_sel_steps_render"))
                ),
                shiny::verticalLayout(
                    shiny::uiOutput(ns("sandbox_clustering_method_render")),
                    shinyWidgets::pickerInput(
                        inputId = ns("sandbox_select_n_clusters"),
                        label = "Select a number of clusters",
                        choices = NULL,
                        inline = FALSE,
                        options = list(
                            `actions-box` = TRUE,
                            title = "Select/deselect clusters",
                            size = 10,
                            width = "50%",
                            `selected-text-format` = "count > 3"
                        ),
                        multiple = TRUE
                    )
                )
            )
        )),
        shiny::actionButton(ns("fix_config"),
            "Fix Configuration",
            style = "font-size:20px;",
            class = "btn-danger"
        )
    )
}

ui_sandbox_metadata_panel <- function(id, draw_line) {
    ns <- shiny::NS(id)
    style <- ifelse(draw_line, "border-right:5px solid", "")

    shinyWidgets::panel(
        style = style,
        shiny::selectizeInput(
            inputId = ns("metadata"),
            label = "Metadata",
            choices = NULL
        ),
        shiny::verticalLayout(
            # shiny::column(6,
            shiny::splitLayout(
                cellWidths = c("40px", "40px"),
                gear_umaps(ns, "metadata"),
                gear_download(ns, "metadata", "metadata")
            ),
            shinyWidgets::pickerInput(
                inputId = ns("select_groups"),
                choices = NULL,
                inline = FALSE,
                # width = "100%",
                # width = "30%",
                options = list(
                    `actions-box` = TRUE,
                    title = "Select/deselect groups",
                    # actionsBox = TRUE,
                    size = 10,
                    width = "90%",
                    `selected-text-format` = "count > 3"
                ),
                multiple = TRUE
            )
            # ),
        ),
        # ),
        # shiny::uiOutput(ns("umap_metadata_generator"))
        shiny::plotOutput(ns("umap_metadata"), height = "auto")
    )
}

ui_sandbox_gene_panel <- function(id, draw_line) {
    ns <- shiny::NS(id)
    style <- ifelse(draw_line, "border-right:5px solid", "")

    shinyWidgets::panel(
        style = style,
        shiny::verticalLayout(
            shiny::selectizeInput(
                inputId = ns("gene_expr"),
                choices = NULL,
                label = "Gene name(s)",
                width = "95%",
                multiple = TRUE
            ),
            shiny::sliderInput(
                inputId = ns("expr_threshold"),
                label = "Gene expression threshold",
                min = 0, max = 10, value = 0,
                width = "95%"
            )
        ),
        shiny::splitLayout(
            cellWidths = c("40px", "40px"),
            gear_umaps(ns, "gene"),
            gear_download(ns, "gene", "gene"),
        ),
        shiny::plotOutput(ns("umap_gene"), height = "auto")
    )
}

ui_sandbox_jsi_panel <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h1("Jaccard Simmilarity Index (JSI)/Cells per cluster"),
        shinyWidgets::dropMenu(shinyWidgets::circleButton(ns("Info"), status = "success", icon = shiny::icon("info"), size = "sm"),
            shiny::h3(shiny::strong("Jaccard Simmilarity Index (JSI) between clusters")),
            shiny::br(),
            shiny::h5("This plot aims to showcase the behaviour of the individual clusters on the different partitions. JSI is calculated for the cell barcodes for every cluster, in both configurations, in a pair-wise manner."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", shiny::a("https://github.com/Core-Bioinformatics/ClustAssess", href = "https://github.com/Core-Bioinformatics/ClustAssess", target = "_blank")),
            placement = "right",
            arrow = FALSE,
            maxWidth = "700px"
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("download"),
            status = "success",
            size = "sm",
            shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
            shiny::textInput(ns("filename_heatmap"), "File name:", width = "80%"),
            shiny::numericInput(ns("width_heatmap"), "Width (in):", 7, 3, 100, 0.1),
            shiny::numericInput(ns("height_heatmap"), "Height (in):", 7, 3, 100, 0.1),
            shiny::selectInput(ns("heatmap_filetype"), "Filetype", choices = c("PDF", "PNG", "SVG"), selected = "PDF", width = "100%"),
            shiny::downloadButton(ns("download_heatmap"), label = "Download Plot")
        ),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("cog"),
            status = "success",
            size = "sm",
            shiny::radioButtons(ns("heatmap_type"), "Calculate similarity", choices = c("JSI", "Cells per cluster"), width = "100%")
        ),
        shiny::selectizeInput(
            inputId = ns("jsi_k_1"),
            label = "Select the number of clusters (k) for the first comparison",
            choices = NULL
        ),
        shiny::selectizeInput(
            inputId = ns("jsi_k_2"),
            label = "Select the number of clusters (k) for the second comparison",
            choices = NULL
        ),
        shiny::plotOutput(ns("barcode_heatmap"), height = "auto")
    )
}

#' UI - Sandbox module
#'
#' @description Creates the UI interface for the sandbox module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to identify the UI elements.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
ui_sandbox <- function(id) {
    ns <- shiny::NS(id)
    shiny::tabPanel(
        "Sandbox",
        ui_colour_picker(ns("colour_picker")),
        shiny::fluidRow(
            shiny::h1("Compare your current configuration", style = "margin-bottom:10px ")
        ),
        shiny::splitLayout(
            cellWidths = c("48%", "48%"),
            ui_sandbox_config_choice(ns("config_choice_left"), TRUE),
            ui_sandbox_config_choice(ns("config_choice_right"), FALSE)
        ),
        shiny::splitLayout(
            cellWidths = c("48%", "48%"),
            ui_sandbox_metadata_panel(ns("sbx_metadata_panel_left"), TRUE),
            ui_sandbox_metadata_panel(ns("sbx_metadata_panel_right"), FALSE)
        ),
        shiny::splitLayout(
            cellWidths = c("48%", "48%"),
            ui_sandbox_gene_panel(ns("sbx_gene_panel_left"), TRUE),
            ui_sandbox_gene_panel(ns("sbx_gene_panel_right"), FALSE)
        ),
        ui_sandbox_jsi_panel(ns("sbx_jsi")),
        style = "margin-left: 25px;margin-top:130px; margin-bottom:70px;"
    )
}
####### SERVER #######
server_colour_picker <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::observe({
                chosen_preset <- input$preset

                if (chosen_preset == "pretty") {
                    shiny::updateSliderInput(session, "hue", value = c(0, 360))
                    shiny::updateSliderInput(session, "saturation", value = c(0.1, 0.5))
                    shiny::updateSliderInput(session, "lightness", value = c(0.5, 0.85))
                    return()
                }

                if (chosen_preset == "pretty_dark") {
                    shiny::updateSliderInput(session, "hue", value = c(0, 360))
                    shiny::updateSliderInput(session, "saturation", value = c(0.1, 0.5))
                    shiny::updateSliderInput(session, "lightness", value = c(0.2, 0.4))
                    return()
                }

                if (chosen_preset == "pastels") {
                    shiny::updateSliderInput(session, "hue", value = c(0, 360))
                    shiny::updateSliderInput(session, "saturation", value = c(0.2, 0.4))
                    shiny::updateSliderInput(session, "lightness", value = c(0.8, 0.9))
                    return()
                }

                shiny::updateSliderInput(session, "hue", value = c(0, 360))
                shiny::updateSliderInput(session, "saturation", value = c(0.1, 0.3))
                shiny::updateSliderInput(session, "lightness", value = c(0.3, 0.6))
            }) %>% shiny::bindEvent(input$preset)

            shiny::observe({
                pallette <- qualpalr::qualpal(
                    n = input$n_colours,
                    colorspace = list(
                        h = input$hue,
                        s = input$saturation,
                        l = input$lightness
                    )
                )$hex

                output$palette_output <- shiny::renderPrint({
                    paste(pallette, collapse = ",")
                })

                output$palette_plot <- shiny::renderPlot({
                    graphics::par(mar = c(0, 0, 0, 0))
                    graphics::barplot(
                        rep(1, length(pallette)),
                        col = pallette,
                        axes = FALSE,
                        space = 0,
                        border = NA
                    )
                })

            }) %>% shiny::bindEvent(input$generate_palette)

            shiny::observe({
                shinyjs::disable("update_colours")
                ngroups <- input$n_groups
                colours_hexcode <- strsplit(input$colours_input, ",")[[1]]

                if (ngroups != length(colours_hexcode)) {
                    return()
                }

                shinyjs::enable("update_colours")
            })

            shiny::observe({
                ngroups <- input$n_groups
                colours_hexcode <- strsplit(input$colours_input, ",")[[1]]
                discrete_colors <- pkg_env$discrete_colors

                shiny::req(ngroups == length(colours_hexcode)) 

                discrete_colors[[as.character(ngroups)]] <- colours_hexcode
                add_env_variable("discrete_colors", discrete_colors)
            }) %>% shiny::bindEvent(input$update_colours)

            output$download_colours <- shiny::downloadHandler(
                filename = function() {
                    "clustassess_app_discrete_colours.json"
                },
                content = function(file) {
                    jsonlite::write_json(pkg_env$discrete_colors, file, pretty = TRUE)
                }
            )

            shiny::observe({
                shiny::req(input$upload_colours)
                shiny::req(input$upload_colours$datapath)

                discrete_colors <- pkg_env$discrete_colors
                new_discrete_colors <- jsonlite::read_json(input$upload_colours$datapath)
                for (discrete_col in names(new_discrete_colors)) {
                    new_discrete_colors[[discrete_col]] <- as.character(new_discrete_colors[[discrete_col]])
                }

                for (discrete_col in names(new_discrete_colors)) {
                    discrete_colors[[discrete_col]] <- new_discrete_colors[[discrete_col]]
                }
                add_env_variable("discrete_colors", discrete_colors)
            }) %>% shiny::bindEvent(input$upload_colours)
        }
    )
}

server_sandbox_config_choice <- function(id, side) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            output$sandbox_sel_fset_render <- shiny::renderUI({
                ns <- session$ns
                shiny::selectInput(
                    inputId = ns("sandbox_sel_fset"),
                    label = "Select feature - set:",
                    choices = names(pkg_env$fsets$fsets),
                    selected = names(pkg_env$fsets$fsets)[1],
                    multiple = FALSE
                )
            })


            output$sandbox_sel_steps_render <- shiny::renderUI({
                ns <- session$ns
                shiny::req(input$sandbox_sel_fset)
                shiny::selectInput(
                    inputId = ns("sandbox_sel_steps"),
                    label = "Select feature - size:",
                    choices = pkg_env$fsets$fsets[input$sandbox_sel_fset][[1]],
                    selected = pkg_env$fsets$fsets[input$sandbox_sel_fset][[1]][1],
                    multiple = FALSE
                )
            })

            shiny::observeEvent(input$sandbox_sel_steps, {
                shiny::req(input$sandbox_sel_fset)
                shiny::req(input$sandbox_sel_steps)
                add_env_variable("clustering_options", list(
                    clustering_options = rhdf5::h5read("stability.h5", paste(input$sandbox_sel_fset, input$sandbox_sel_steps, "clustering_stability/split_by_k/structure_list", sep = "/"))
                ))
            })

            output$sandbox_clustering_method_render <- shiny::renderUI({
                shiny::req(input$sandbox_sel_fset)
                shiny::req(input$sandbox_sel_steps)
                shiny::req(pkg_env$clustering_options)
                ns <- session$ns
                shiny::radioButtons(
                    inputId = ns("sandbox_clustering_method"),
                    label = "Select clustering Method:",
                    choices = names(pkg_env$clustering_options$clustering_options),
                    selected = names(pkg_env$clustering_options$clustering_options)[1],
                    width = "50%"
                )
            })

            toListen <- shiny::reactive(
                list(
                    input$sandbox_sel_fset,
                    input$sandbox_sel_steps,
                    input$sandbox_clustering_method
                )
            )

            shiny::observeEvent(toListen(), {
                shiny::req(input$sandbox_clustering_method)
                k_values <- pkg_env$clustering_options$clustering_options[input$sandbox_clustering_method][[1]]
                shinyWidgets::updatePickerInput(
                    session = session,
                    inputId = "sandbox_select_n_clusters",
                    choices = k_values,
                    selected = k_values[1]
                )
            })
            user_choice <- shiny::eventReactive(input$fix_config, {
                shiny::req(input$fix_config)
                shiny::req(input$sandbox_sel_fset)
                shiny::req(input$sandbox_sel_steps)
                shiny::req(input$sandbox_select_n_clusters)
                shiny::req(input$sandbox_clustering_method)
                user_choice <- list(
                    fset = input$sandbox_sel_fset,
                    fsize = input$sandbox_sel_steps,
                    k_vals = input$sandbox_select_n_clusters,
                    c_method = input$sandbox_clustering_method,
                    side = side
                )
                user_choice
            })
            shiny::observeEvent(input$fix_config, {
                shiny::req(user_choice())
                if (user_choice()$side == "left") {
                    add_env_variable("stab_obj_left", list(
                        mbs = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "clustering_stability", "split_by_k", "mbs", user_choice()$c_method, sep = "/"))[user_choice()$k_vals],
                        ecc = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "clustering_stability", "split_by_k", "ecc", sep = "/"))[paste(sprintf("%06d", as.integer(user_choice()$k_vals)), user_choice()$c_method, sep = ";")],
                        ecc_order = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "clustering_stability", "split_by_k", "ecc_order", sep = "/"))[paste(sprintf("%06d", as.integer(user_choice()$k_vals)), user_choice()$c_method, sep = ";")],
                        umap = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "umap", sep = "/"))
                    ))
                } else {
                    add_env_variable("stab_obj_right", list(
                        mbs = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "clustering_stability", "split_by_k", "mbs", user_choice()$c_method, sep = "/"))[user_choice()$k_vals],
                        ecc = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "clustering_stability", "split_by_k", "ecc", sep = "/"))[paste(sprintf("%06d", as.integer(user_choice()$k_vals)), user_choice()$c_method, sep = ";")],
                        ecc_order = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "clustering_stability", "split_by_k", "ecc_order", sep = "/"))[paste(sprintf("%06d", as.integer(user_choice()$k_vals)), user_choice()$c_method, sep = ";")],
                        umap = rhdf5::h5read("stability.h5", paste(user_choice()$fset, user_choice()$fsize, "umap", sep = "/"))
                    ))
                }
            })
            shiny::observeEvent(input$fix_config, {
                shiny::req(user_choice())
                if (user_choice()$side == "left") {
                    add_env_variable("selected_kvals_left", user_choice()$k_vals)
                } else {
                    add_env_variable("selected_kvals_right", user_choice()$k_vals)
                }
            })
        }
    )
}

server_sandbox_metadata_panel_left <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            plt_height <- shiny::reactive(
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.43))
            )
            shiny::observe({
                is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
                is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

                if (is_cluster || is_ecc) {
                    if (is_cluster) {
                        k_value <- as.numeric(strsplit(input$metadata, "_")[[1]][2])
                        shinyjs::show(id = "select_groups")
                        shinyWidgets::updatePickerInput(
                            session,
                            inputId = "select_groups",
                            choices = seq_len(k_value),
                            selected = seq_len(k_value)
                        )
                    } else {
                        shinyjs::hide(id = "select_groups")
                    }

                    return()
                }

                mtd_names <- pkg_env$metadata_unique[[input$metadata]]
                if (is.null(mtd_names)) {
                    shinyjs::hide(id = "select_groups")
                } else {
                    shinyjs::show(id = "select_groups")
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = "select_groups",
                        choices = mtd_names,
                        selected = mtd_names
                    )
                }
            }) %>% shiny::bindEvent(input$metadata)

            metadata_legend_height <- shiny::reactive({
                unique_values <- pkg_env$metadata_unique[[input$metadata]]
                # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
                grDevices::pdf(file = NULL, width = plt_height(), height = plt_height())
                if (is.null(unique_values)) {
                    graphics::par(mai = c(0.1, 0, 0.1, 0))
                    text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$metadata_text_size)
                    grDevices::dev.off()
                    return(text_height * ppi * 1.25)
                }

                graphics::par(mar = c(0, 0, 0, 0))
                predicted_width <- graphics::strwidth(c(" ", unique_values), units = "inches", cex = input$metadata_text_size) * ppi
                space_width <- predicted_width[1]
                predicted_width <- predicted_width[2:length(predicted_width)]

                number_columns <- min(
                    max(
                        plt_height() %/% (5 * space_width + max(predicted_width)),
                        1
                    ),
                    length(unique_values)
                )
                number_rows <- ceiling(length(unique_values) / number_columns)
                print(number_rows)

                text_height <- graphics::strheight(
                    paste(
                        rep("TEXT", number_rows + 1),
                        collapse = "\n"
                    ),
                    units = "inches",
                    cex = input$metadata_text_size
                )

                grDevices::dev.off()

                return(text_height * ppi * 1.25)
            })


            output$umap_metadata <- shiny::renderPlot(
                height = function() {
                    plt_height() + metadata_legend_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    shiny::req(input$metadata)
                    is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
                    is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

                    if (is_cluster || is_ecc) {
                        k <- strsplit(input$metadata, "_")[[1]][2]
                        if (is_ecc) {
                            cl_method <- strsplit(names(pkg_env$stab_obj_left$ecc)[1], ";")[[1]][2]
                            print(cl_method)
                            unique_values <- NULL
                            color_values <- NULL
                            color_info <- pkg_env$stab_obj_left$ecc[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]
                        } else {
                            color_values <- rhdf5::h5read("stability.h5", paste0("colors/", k))
                            unique_values <- seq_len(as.integer(k))
                            color_info <- pkg_env$stab_obj_left$mbs[[k]]
                        }
                    } else {
                        unique_values <- pkg_env$metadata_unique[[input$metadata]]
                        color_values <- pkg_env$metadata_colors[[input$metadata]]
                        color_info <- pkg_env$metadata[[input$metadata]]
                    }

                    color_plot2(
                        embedding = pkg_env$stab_obj_left$umap,
                        color_info = color_info,
                        color_values = color_values,
                        unique_values = unique_values,
                        plt_height = plt_height(),
                        plt_width = plt_height(),
                        display_legend = TRUE,
                        predicted_height = (metadata_legend_height() - 1) / ppi,
                        pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
                        pt_size = input$metadata_pt_size,
                        text_size = input$metadata_text_size,
                        axis_size = input$metadata_axis_size,
                        labels = input$metadata_labels,
                        groups_highlight = input$select_groups
                    )
                }
            )
            plot_data <- shiny::reactive({
                shiny::req(input$metadata)

                is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
                is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

                if (is_cluster || is_ecc) {
                    k <- strsplit(input$metadata, "_")[[1]][2]
                    if (is_ecc) {
                        cl_method <- strsplit(names(pkg_env$stab_obj_left$ecc)[1], ";")[[1]][2]
                        unique_values <- NULL
                        color_values <- NULL
                        color_info <- pkg_env$stab_obj_left$ecc[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]
                        color_info <- color_info[pkg_env$stab_obj_left$ecc_order[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]]
                    } else {
                        color_values <- rhdf5::h5read("stability.h5", paste0("colors/", k))
                        unique_values <- seq_len(as.integer(k))
                        color_info <- pkg_env$stab_obj_left$mbs[[k]]
                    }
                } else {
                    unique_values <- pkg_env$metadata_unique[[input$metadata]]
                    color_values <- pkg_env$metadata_colors[[input$metadata]]
                    color_info <- pkg_env$metadata[[input$metadata]]
                }

                list(
                    unique_values = unique_values,
                    color_values = color_values,
                    color_info = color_info
                )
            }) %>% shiny::bindEvent(input$metadata)

            output$download_metadata <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_metadata, ".", tolower(input$filetype_metadata))
                },
                content = function(file) {
                    shiny::req(input$metadata, input$width_metadata, input$height_metadata)
                    filetypes[[input$filetype_metadata]](file, width = input$width_metadata, height = input$height_metadata)
                    color_plot2(
                        embedding = pkg_env$stab_obj_left$umap,
                        color_info = plot_data()$color_info,
                        color_values = plot_data()$color_values,
                        unique_values = plot_data()$unique_values,
                        plt_height = input$height_metadata * ppi, # - metadata_legend_height(),
                        plt_width = input$width_metadata * ppi,
                        # predicted_height = (metadata_legend_height() - 1) / ppi,
                        pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
                        pt_size = input$metadata_pt_size,
                        text_size = input$metadata_text_size,
                        axis_size = input$metadata_axis_size,
                        legend_text_size = input$metadata_legend_size,
                        labels = input$metadata_labels,
                        groups_highlight = input$select_groups,
                        display_legend = TRUE
                    )
                    grDevices::dev.off()
                }
            )
        }
    )
}
server_sandbox_metadata_panel_right <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            plt_height <- shiny::reactive(
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.43))
            )
            shiny::observe({
                is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
                is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

                if (is_cluster || is_ecc) {
                    if (is_cluster) {
                        k_value <- as.numeric(strsplit(input$metadata, "_")[[1]][2])
                        shinyjs::show(id = "select_groups")
                        shinyWidgets::updatePickerInput(
                            session,
                            inputId = "select_groups",
                            choices = seq_len(k_value),
                            selected = seq_len(k_value)
                        )
                    } else {
                        shinyjs::hide(id = "select_groups")
                    }

                    return()
                }

                mtd_names <- pkg_env$metadata_unique[[input$metadata]]
                if (is.null(mtd_names)) {
                    shinyjs::hide(id = "select_groups")
                } else {
                    shinyjs::show(id = "select_groups")
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = "select_groups",
                        choices = mtd_names,
                        selected = mtd_names
                    )
                }
            }) %>% shiny::bindEvent(input$metadata)

            metadata_legend_height <- shiny::reactive({
                unique_values <- pkg_env$metadata_unique[[input$metadata]]
                # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
                grDevices::pdf(file = NULL, width = plt_height(), height = plt_height())
                if (is.null(unique_values)) {
                    graphics::par(mai = c(0.1, 0, 0.1, 0))
                    text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$metadata_text_size)
                    grDevices::dev.off()
                    return(text_height * ppi * 1.25)
                }

                graphics::par(mar = c(0, 0, 0, 0))
                predicted_width <- graphics::strwidth(c(" ", unique_values), units = "inches", cex = input$metadata_text_size) * ppi
                space_width <- predicted_width[1]
                predicted_width <- predicted_width[2:length(predicted_width)]

                number_columns <- min(
                    max(
                        plt_height() %/% (5 * space_width + max(predicted_width)),
                        1
                    ),
                    length(unique_values)
                )
                number_rows <- ceiling(length(unique_values) / number_columns)
                print(number_rows)

                text_height <- graphics::strheight(
                    paste(
                        rep("TEXT", number_rows + 1),
                        collapse = "\n"
                    ),
                    units = "inches",
                    cex = input$metadata_text_size
                )

                grDevices::dev.off()

                return(text_height * ppi * 1.25)
            })


            output$umap_metadata <- shiny::renderPlot(
                height = function() {
                    plt_height() + metadata_legend_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    shiny::req(input$metadata)
                    is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
                    is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

                    if (is_cluster || is_ecc) {
                        k <- strsplit(input$metadata, "_")[[1]][2]
                        if (is_ecc) {
                            cl_method <- strsplit(names(pkg_env$stab_obj_right$ecc)[1], ";")[[1]][2]
                            print(cl_method)
                            unique_values <- NULL
                            color_values <- NULL
                            color_info <- pkg_env$stab_obj_right$ecc[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]
                        } else {
                            color_values <- rhdf5::h5read("stability.h5", paste0("colors/", k))
                            unique_values <- seq_len(as.integer(k))
                            color_info <- pkg_env$stab_obj_right$mbs[[k]]
                        }
                    } else {
                        unique_values <- pkg_env$metadata_unique[[input$metadata]]
                        color_values <- pkg_env$metadata_colors[[input$metadata]]
                        color_info <- pkg_env$metadata[[input$metadata]]
                    }

                    color_plot2(
                        embedding = pkg_env$stab_obj_right$umap,
                        color_info = color_info,
                        color_values = color_values,
                        unique_values = unique_values,
                        plt_height = plt_height(),
                        plt_width = plt_height(),
                        display_legend = TRUE,
                        predicted_height = (metadata_legend_height() - 1) / ppi,
                        pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
                        pt_size = input$metadata_pt_size,
                        text_size = input$metadata_text_size,
                        axis_size = input$metadata_axis_size,
                        labels = input$metadata_labels,
                        groups_highlight = input$select_groups
                    )
                }
            )

            plot_data <- shiny::reactive({
                shiny::req(input$metadata)

                is_cluster <- stringr::str_detect(input$metadata, "stable_[0-9]+_clusters")
                is_ecc <- stringr::str_detect(input$metadata, "ecc_[0-9]+")

                if (is_cluster || is_ecc) {
                    k <- strsplit(input$metadata, "_")[[1]][2]
                    if (is_ecc) {
                        cl_method <- strsplit(names(pkg_env$stab_obj_right$ecc)[1], ";")[[1]][2]
                        unique_values <- NULL
                        color_values <- NULL
                        color_info <- pkg_env$stab_obj_right$ecc[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]
                        color_info <- color_info[pkg_env$stab_obj_rightt$ecc_order[[paste(sprintf("%06d", as.integer(k)), cl_method, sep = ";")]]]
                    } else {
                        color_values <- rhdf5::h5read("stability.h5", paste0("colors/", k))
                        unique_values <- seq_len(as.integer(k))
                        color_info <- pkg_env$stab_obj_right$mbs[[k]]
                    }
                } else {
                    unique_values <- pkg_env$metadata_unique[[input$metadata]]
                    color_values <- pkg_env$metadata_colors[[input$metadata]]
                    color_info <- pkg_env$metadata[[input$metadata]]
                }

                list(
                    unique_values = unique_values,
                    color_values = color_values,
                    color_info = color_info
                )
            }) %>% shiny::bindEvent(input$metadata)

            output$download_metadata <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_metadata, ".", tolower(input$filetype_metadata))
                },
                content = function(file) {
                    shiny::req(input$metadata, input$width_metadata, input$height_metadata)
                    filetypes[[input$filetype_metadata]](file, width = input$width_metadata, height = input$height_metadata)
                    color_plot2(
                        embedding = pkg_env$stab_obj_right$umap,
                        color_info = plot_data()$color_info,
                        color_values = plot_data()$color_values,
                        unique_values = plot_data()$unique_values,
                        plt_height = input$height_metadata * ppi, # - metadata_legend_height(),
                        plt_width = input$width_metadata * ppi,
                        # predicted_height = (metadata_legend_height() - 1) / ppi,
                        pch = ifelse(input$metadata_pt_type == "Pixel", ".", 19),
                        pt_size = input$metadata_pt_size,
                        text_size = input$metadata_text_size,
                        axis_size = input$metadata_axis_size,
                        legend_text_size = input$metadata_legend_size,
                        labels = input$metadata_labels,
                        groups_highlight = input$select_groups,
                        display_legend = TRUE
                    )
                    grDevices::dev.off()
                }
            )
        }
    )
}

server_sandbox_gene_panel_left <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            expr_matrix <- shiny::reactive({
                if ("genes" %in% names(pkg_env)) {
                    index_gene <- pkg_env$genes[input$gene_expr]
                    index_gene <- index_gene[!is.na(index_gene)] # not necesarry most probably

                    return(rhdf5::h5read("expression.h5", "expression_matrix", index = list(index_gene, NULL)))
                }

                # for backward-compatibility purposes
                index_interest <- pkg_env$genes_of_interest[input$gene_expr]
                index_interest <- index_interest[!is.na(index_interest)]

                index_others <- pkg_env$genes_others[input$gene_expr]
                index_others <- index_others[!is.na(index_others)]

                rbind(
                    rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(index_interest, NULL)),
                    rhdf5::h5read("expression.h5", "matrix_others", index = list(index_others, NULL))
                )
            }) %>% shiny::bindEvent(input$gene_expr)

            max_level_expr <- shiny::reactive(max(expr_matrix()))

            shiny::observe({
                shiny::updateSliderInput(session,
                    inputId = "expr_threshold",
                    max = round(max_level_expr(), 3),
                    step = round(max_level_expr() / 10, 3)
                )
            }) %>% shiny::bindEvent(input$gene_expr)

            plt_height <- shiny::reactive(
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.43))
            )

            gene_legend_height <- shiny::reactive({
                # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
                grDevices::pdf(NULL, width = plt_height(), height = plt_height())
                graphics::par(mai = c(0.1, 0, 0.1, 0))
                text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$gene_text_size)
                grDevices::dev.off()
                return((0.2 + text_height) * ppi)
            })

            output$umap_gene <- shiny::renderPlot(
                height = function() {
                    plt_height() + gene_legend_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    shiny::req(input$expr_threshold, input$gene_expr, expr_matrix())

                    unique_values <- NULL
                    used_matrix <- expr_matrix()
                    color_values <- function(n) {
                        grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n)
                    }
                    if (length(input$gene_expr) > 1) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- matrixStats::colSums2(used_matrix >= input$expr_threshold) == length(input$gene_expr)
                    } else if (input$expr_threshold > 0) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- used_matrix >= input$expr_threshold
                    }

                    color_plot2(
                        embedding = pkg_env$stab_obj_left$umap,
                        color_info = used_matrix,
                        plt_height = plt_height(),
                        plt_width = plt_height(),
                        predicted_height = (gene_legend_height() - 1) / ppi,
                        display_legend = TRUE,
                        unique_values = unique_values,
                        color_values = color_values,
                        pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
                        pt_size = input$gene_pt_size,
                        text_size = input$gene_text_size
                    )
                }
            )

            output$download_gene <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_gene, ".", tolower(input$filetype_gene))
                },
                content = function(file) {
                    shiny::req(input$expr_threshold, input$gene_expr, input$width_gene, input$height_gene, expr_matrix())
                    filetypes[[input$filetype_gene]](file, width = input$width_gene, height = input$height_gene)
                    unique_values <- NULL
                    used_matrix <- expr_matrix()
                    color_values <- function(n) {
                        grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n)
                    }
                    if (length(input$gene_expr) > 1) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- matrixStats::colSums2(used_matrix > input$expr_threshold) >= (length(input$gene_expr) - input$relaxation)
                    } else if (input$expr_threshold > 0) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- used_matrix > input$expr_threshold
                    }

                    color_plot2(
                        embedding = pkg_env$stab_obj_left$umap,
                        color_info = used_matrix,
                        plt_height = input$height_gene * ppi, #- gene_legend_height(),
                        plt_width = input$width_gene * ppi,
                        predicted_height = (gene_legend_height() - 1) / ppi,
                        # color_values = function(n) { paletteer::paletteer_c("grDevices::OrRd", n)},
                        unique_values = unique_values,
                        color_values = color_values,
                        pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
                        pt_size = input$gene_pt_size,
                        text_size = input$gene_text_size,
                        legend_text_size = input$gene_legend_size,
                        axis_size = input$gene_axis_size,
                        display_legend = TRUE
                    )
                    grDevices::dev.off()
                }
            )
        }
    )
}
server_sandbox_gene_panel_right <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            expr_matrix <- shiny::reactive({
                if ("genes" %in% names(pkg_env)) {
                    index_gene <- pkg_env$genes[input$gene_expr]
                    index_gene <- index_gene[!is.na(index_gene)] # not necesarry most probably

                    return(rhdf5::h5read("expression.h5", "expression_matrix", index = list(index_gene, NULL)))
                }

                # for backward-compatibility purposes
                index_interest <- pkg_env$genes_of_interest[input$gene_expr]
                index_interest <- index_interest[!is.na(index_interest)]

                index_others <- pkg_env$genes_others[input$gene_expr]
                index_others <- index_others[!is.na(index_others)]

                rbind(
                    rhdf5::h5read("expression.h5", "matrix_of_interest", index = list(index_interest, NULL)),
                    rhdf5::h5read("expression.h5", "matrix_others", index = list(index_others, NULL))
                )
            }) %>% shiny::bindEvent(input$gene_expr)

            max_level_expr <- shiny::reactive(max(expr_matrix()))

            shiny::observe({
                shiny::updateSliderInput(session,
                    inputId = "expr_threshold",
                    max = round(max_level_expr(), 3),
                    step = round(max_level_expr() / 10, 3)
                )
            }) %>% shiny::bindEvent(input$gene_expr)

            plt_height <- shiny::reactive(
                floor(min(pkg_env$height_ratio * pkg_env$dimension()[2], pkg_env$dimension()[1] * 0.43))
            )

            gene_legend_height <- shiny::reactive({
                # ragg::agg_png(res = ppi, width = plt_height(), height = plt_height())
                grDevices::pdf(NULL, width = plt_height(), height = plt_height())
                graphics::par(mai = c(0.1, 0, 0.1, 0))
                text_height <- graphics::strheight("TE\nXT\n", units = "inches", cex = input$gene_text_size)
                grDevices::dev.off()
                return((0.2 + text_height) * ppi)
            })

            output$umap_gene <- shiny::renderPlot(
                height = function() {
                    plt_height() + gene_legend_height()
                },
                width = function() {
                    plt_height()
                },
                {
                    shiny::req(input$expr_threshold, input$gene_expr, expr_matrix())

                    unique_values <- NULL
                    used_matrix <- expr_matrix()
                    color_values <- function(n) {
                        grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n)
                    }
                    if (length(input$gene_expr) > 1) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- matrixStats::colSums2(used_matrix >= input$expr_threshold) == length(input$gene_expr)
                    } else if (input$expr_threshold > 0) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- used_matrix >= input$expr_threshold
                    }

                    color_plot2(
                        embedding = pkg_env$stab_obj_right$umap,
                        color_info = used_matrix,
                        plt_height = plt_height(),
                        plt_width = plt_height(),
                        predicted_height = (gene_legend_height() - 1) / ppi,
                        display_legend = TRUE,
                        unique_values = unique_values,
                        color_values = color_values,
                        pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
                        pt_size = input$gene_pt_size,
                        text_size = input$gene_text_size
                    )
                }
            )

            output$download_gene <- shiny::downloadHandler(
                filename = function() {
                    paste0(input$filename_gene, ".", tolower(input$filetype_gene))
                },
                content = function(file) {
                    shiny::req(input$expr_threshold, input$gene_expr, input$width_gene, input$height_gene, expr_matrix())
                    filetypes[[input$filetype_gene]](file, width = input$width_gene, height = input$height_gene)
                    unique_values <- NULL
                    used_matrix <- expr_matrix()
                    color_values <- function(n) {
                        grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "OrRd"))(n)
                    }
                    if (length(input$gene_expr) > 1) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- matrixStats::colSums2(used_matrix > input$expr_threshold) >= (length(input$gene_expr) - input$relaxation)
                    } else if (input$expr_threshold > 0) {
                        unique_values <- c("other", "cells above threshold")
                        color_values <- c("lightgray", "red")
                        used_matrix <- used_matrix > input$expr_threshold
                    }

                    color_plot2(
                        embedding = pkg_env$stab_obj_right$umap,
                        color_info = used_matrix,
                        plt_height = input$height_gene * ppi, #- gene_legend_height(),
                        plt_width = input$width_gene * ppi,
                        predicted_height = (gene_legend_height() - 1) / ppi,
                        # color_values = function(n) { paletteer::paletteer_c("grDevices::OrRd", n)},
                        unique_values = unique_values,
                        color_values = color_values,
                        pch = ifelse(input$gene_pt_type == "Pixel", ".", 19),
                        pt_size = input$gene_pt_size,
                        text_size = input$gene_text_size,
                        legend_text_size = input$gene_legend_size,
                        axis_size = input$gene_axis_size,
                        display_legend = TRUE
                    )
                    grDevices::dev.off()
                }
            )
        }
    )
}

server_sandbox_jsi <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            plt_width <- shiny::reactive(
                pkg_env$dimension()[1]
            )

            plt_height <- shiny::reactive(
                floor(pkg_env$height_ratio * pkg_env$dimension()[2])
            )

            barcode_heatmap <- shiny::reactive({
                shiny::req(input$jsi_k_1, input$jsi_k_2)
                if (!is.na(as.numeric(input$jsi_k_1))) {
                    clustering_1 <- as.matrix(pkg_env$stab_obj$mbs[[as.character(input$jsi_k_1)]])
                    df_1 <- data.frame(clustering_1)
                } else {
                    meta_category <- pkg_env$metadata[, input$jsi_k_1]
                    df_1 <- data.frame(meta_category)
                }
                df_1$cell <- rownames(df_1)
                if (!is.na(as.numeric(input$jsi_k_2))) {
                    clustering_2 <- as.matrix(pkg_env$stab_obj$mbs[[as.character(input$jsi_k_2)]])
                    df_2 <- data.frame(clustering_2)
                } else {
                    meta_category <- pkg_env$metadata[, input$jsi_k_2]
                    df_2 <- data.frame(meta_category)
                }
                df_2$cell <- rownames(df_2)
                all_clusters_1 <- unique(df_1[, 1])
                all_clusters_2 <- unique(df_2[, 1])

                mat <- matrix(,
                    nrow = length(all_clusters_2),
                    ncol = length(all_clusters_1)
                )
                colnames(mat) <- sort(all_clusters_1)
                rownames(mat) <- sort(all_clusters_2)
                if (input$heatmap_type == "JSI") {
                    for (m in all_clusters_1) {
                        cluster_1 <- rownames(df_1[df_1[, 1] == m, ])
                        for (n in all_clusters_2) {
                            cluster_2 <- rownames(df_2[df_2[, 1] == n, ])
                            mat[as.character(n), as.character(m)] <- jaccard_index(cluster_1, cluster_2)
                            label <- "JSI"
                        }
                    }
                } else {
                    for (m in all_clusters_1) {
                        cluster_1 <- rownames(df_1[df_1[, 1] == m, ])
                        for (n in all_clusters_2) {
                            cluster_2 <- rownames(df_2[df_2[, 1] == n, ])
                            mat[as.character(n), as.character(m)] <- length(intersect(cluster_1, cluster_2))
                            label <- "Shared cells"
                        }
                    }
                }
                df_mat <- reshape2::melt(mat)

                ggplot2::ggplot(df_mat, ggplot2::aes(as.factor(.data$Var1), as.factor(.data$Var2))) +
                    ggplot2::geom_tile(ggplot2::aes(fill = .data$value)) +
                    ggplot2::geom_text(ggplot2::aes(label = round(.data$value, 2))) +
                    ggplot2::scale_fill_gradient2(
                        low = scales::muted("darkred"),
                        mid = "white",
                        high = scales::muted("green"),
                        midpoint = 0
                    ) +
                    ggplot2::theme(
                        panel.background = ggplot2::element_rect(fill = "white"),
                        axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1, size = 10, face = "bold"),
                        axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
                        axis.title = ggplot2::element_text(size = 14, face = "bold"),
                        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 20, l = 30)),
                        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 20, b = 30))
                    ) +
                    ggplot2::xlab("Clusters in Configuration 2") +
                    ggplot2::ylab("Clusters in Configuration 1") +
                    ggplot2::labs(fill = label)
            })

            output$barcode_heatmap <- shiny::renderPlot(
                {
                    if (!shiny::isTruthy(input$jsi_k_1) | !shiny::isTruthy(input$jsi_k_2)) {
                        return(ggplot2::ggplot() +
                            ggplot2::theme_void())
                    } else {
                        barcode_heatmap()
                    }
                },
                height = plt_height(),
                width = plt_width()
            )

            heatmap_filetype <- shiny::reactive({
                if (input$heatmap_filetype == "PDF") {
                    filename <- paste0(input$filename_heatmap, ".pdf")
                    return(filename)
                } else if (input$heatmap_filetype == "PNG") {
                    filename <- paste0(input$filename_heatmap, ".png")
                    return(filename)
                } else {
                    filename <- paste0(input$filename_heatmap, ".svg")
                    return(filename)
                }
            })

            output$download_heatmap <- shiny::downloadHandler(
                filename = function() {
                    heatmap_filetype()
                },
                content = function(file) {
                    ggplot2::ggsave(file, barcode_heatmap(),
                        width = input$width_heatmap,
                        height = input$height_heatmap,
                        units = "in",
                        limitsize = FALSE
                    )
                }
            )
        }
    )
}

#' Server - Sandbox module
#'
#' @description Creates the backend interface for the sandbox module inside
#' the ClustAssess Shiny application.
#'
#' @param id The id of the module, used to acess the UI elements.
#'
#' @note This function should not be called directly, but in the context of the
#' app that is created using the `write_shiny_app` function.
#'
#' @export
server_sandbox <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            server_colour_picker("colour_picker")

            add_env_variable("fsets", list(
                fsets = rhdf5::h5read("stability.h5", "feature_ordering/stable")
            ))
            genes <- rhdf5::h5read("expression.h5", "genes")
            index <- seq_along(genes)
            names(index) <- genes
            add_env_variable("genes", index)


            server_sandbox_config_choice("config_choice_left", "left")
            shiny::observeEvent(input$"config_choice_left-fix_config", {
                shiny::req(pkg_env$selected_kvals_left)
                shiny::req(pkg_env$clustering_options)
                shiny::updateSelectizeInput(
                    session,
                    inputId = "sbx_metadata_panel_left-metadata",
                    server = FALSE,
                    choices = c(colnames(pkg_env$metadata), paste0("stable_", pkg_env$selected_kvals_left, "_clusters"), paste0("ecc_", pkg_env$selected_kvals_left)),
                    selected = paste0("stable_", pkg_env$selected_kvals_left[1], "_clusters")
                )

                shiny::updateSelectizeInput(
                    session,
                    inputId = glue::glue("sbx_gene_panel_left-gene_expr"),
                    choices = c(names(pkg_env$genes)),
                    # selected = NULL,
                    selected = names(pkg_env$genes[1]),
                    server = TRUE,
                    options = list(
                        maxOptions = 7,
                        create = TRUE,
                        persist = TRUE
                    )
                )
            })

            server_sandbox_config_choice("config_choice_right", "right")
            shiny::observeEvent(input$"config_choice_right-fix_config", {
                shiny::req(pkg_env$selected_kvals_right)
                shiny::updateSelectizeInput(
                    session,
                    inputId = "sbx_metadata_panel_right-metadata",
                    server = FALSE,
                    choices = c(colnames(pkg_env$metadata), paste0("stable_", pkg_env$selected_kvals_right, "_clusters"), paste0("ecc_", pkg_env$selected_kvals_right)),
                    selected = paste0("stable_", pkg_env$selected_kvals_right[1], "_clusters")
                )

                shiny::updateSelectizeInput(
                    session,
                    inputId = glue::glue("sbx_gene_panel_right-gene_expr"),
                    choices = c(names(pkg_env$genes)),
                    # selected = NULL,
                    selected = names(pkg_env$genes[1]),
                    server = TRUE,
                    options = list(
                        maxOptions = 7,
                        create = TRUE,
                        persist = TRUE
                    )
                )
            })

            server_sandbox_metadata_panel_left("sbx_metadata_panel_left")
            server_sandbox_metadata_panel_right("sbx_metadata_panel_right")
            server_sandbox_gene_panel_left("sbx_gene_panel_left")
            server_sandbox_gene_panel_right("sbx_gene_panel_right")

            discrete <- c()
            for (category in colnames(pkg_env$metadata)) {
                if (length(unique(pkg_env$metadata[, category])) < 20) {
                    discrete <- append(discrete, category)
                }
            }

            # And the JSI
            shiny::observeEvent(input$"config_choice_left-fix_config", {
                shiny::req(pkg_env$selected_kvals_left)
                shiny::req(pkg_env$clustering_options)
                shiny::updateSelectizeInput(
                    session = session,
                    inputId = "sbx_jsi-jsi_k_1",
                    choices = append(pkg_env$selected_kvals_left, discrete),
                    selected = pkg_env$selected_kvals_left[1]
                )
            })
            shiny::observeEvent(input$"config_choice_right-fix_config", {
                shiny::req(pkg_env$selected_kvals_right)
                shiny::req(pkg_env$clustering_options)
                shiny::updateSelectizeInput(
                    session = session,
                    inputId = "sbx_jsi-jsi_k_2",
                    choices = append(pkg_env$selected_kvals_right, discrete),
                    selected = pkg_env$selected_kvals_right[1]
                )
            })
            server_sandbox_jsi("sbx_jsi")

            gc()
        }
    )
}
