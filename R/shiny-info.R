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
            shiny::h5("The recommendation is based on a rank-based selection that is done on the overall stability and overall incremental stability."),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

# dr_template_info <- function(session) {
#     shiny::showModal(
#         shiny::modalDialog(
#             shiny::h1("\n"),
#             shiny::h5("For more information please go to:"),
#             shiny::tagList("", github),
#             easyClose = TRUE
#         ),
#         session
#     )
# }

# dr_template_info <- function(session) {
#     shiny::showModal(
#         shiny::modalDialog(
#             shiny::h1("\n"),
#             shiny::h5("For more information please go to:"),
#             shiny::tagList("", github),
#             easyClose = TRUE
#         ),
#         session
#     )
# }

##### GRAPH CONSTRUCTION #####
##### GRAPH CLUSTERING #####
gclust_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Assessing the stability of the graph clustering")),
            shiny::br(),
            shiny::h5("In the third step, the Graph Clustering, the main factor that affects the results is the community detection method that is used, as well as the resolution parameter that is used for the clustering."),
            shiny::h5("\n"),
            shiny::h5("In this page the user can interact and visualise the stability of the mainly used community detection methods - Louvain, Louvain with multilevel refinement, SLM and Leiden, with respect to the number of clusters and resolution. The user can also assess the relationship between the resolution parameter and the number of clusters, as well the stability of the clustering results."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

gclust_distro_res_boxplot_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Asessing the stability of the clustering method for each resolution parameter")),
            shiny::br(),
            shiny::h5("The following plot illustrates as a boxplot the distribution of the stability (measured using the ECC score) for different clustering methods."),
            shiny::h5("\n"),
            shiny::h5("Given the range of resolution parameter used in the assessment, the user can visualise the stability when grouping the partitions by the number of clusters or the resolution."),
            shiny::h1("\n"),
            shiny::h5("Below this plot, the user can visualise the distribution of the obtained clusters and of the stability (as measured using ECC) on an UMAP representation."),
            shiny::br(),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

gclust_distro_overall_boxplot_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Assessing the overall stability of the clustering method")),
            shiny::br(),
            shiny::h5("The following plots summaries the stability for each clustering method by extracting the ECC median from each group. The group can be either the resolution parameter, as shown in the left panel, or the number of clusters, as shown in the right panel."),
            shiny::h5("\n"),
            shiny::h5("These plots should be used by the user to decide which clustering method to use further."),
            shiny::h5("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

gclust_k_corresp_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Finding the relationship between the resolution parameter and the number of clusters")),
            shiny::br(),
            shiny::h5("The choice of the resolution parameter can be done arbitrary, without a good intuition about an estimate number of resulting clusters."),
            shiny::h5("\n"),
            shiny::h5("The following plot showcases the relationship between the resolution parameter and the number of clusters obtained by the clustering method."),
            shiny::h5("\n"),
            shiny::h5("The shape of the point represents the clustering method used. The size of the point is proportional with the frequency of the pair resolution - number of clusters. The colour of the point represents the stability of the clustering method. Lighter colours (yellow) represent a higher stability."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

gclust_k_stab_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Assessing the stability of the number of clusters")),
            shiny::br(),
            shiny::h5("After selecting a clustering method, the final step is to assess the stability of the number of clusters obtained using that method and choose the stable options."),
            shiny::h5("\n"),
            shiny::h5("This plot illustrates the stability of the number of clusters using three metrics:"),
            shiny::h5(" y axis: the number of different partitions having the same number of clusters (k); high values can be a sign that, for this k, the partition's composition is very likely to be changed"),
            shiny::h5(" - point colour: ECC metric showing the stability of the number of clusters; lighter colours (yellow) represent a higher stability; this shows how similar are the different partitions having the same number of clusters"),
            shiny::h5("- point size: the frequency of the number of clusters; this shows how reliable, from a statistical perspective, are the assessments done in the other dimensions; the size is proportionate to the frequency"),
            shiny::h5("\n"),
            shiny::h5("The user can also filter, from the plot settings, which k values to show based on ECC and frequency thresholds."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

gclust_choice_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Choosing the clustering method and number of clusters")),
            shiny::br(),
            shiny::h5("For the downstream analysis to be robust and reproducible, the user should decide upon parameters and number of clusters that are as stable to the seed change as possible."),
            shiny::h5("\n"),
            shiny::h5("Based on the plots from this tab, the user should fix the clustering method and the number of clusters that will be used in the downstream analysis from the next tab."),
            shiny::h5("\n"),
            shiny::h5("For the number of clusters, we recommend using the ECC and the frequency as the most important metrics to use when considering a number of clusters or not. The frequency should be high enough (30, for example) for the assessment to be statistically reliable."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

##### CLUSTERING COMPARISON #####

compar_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Comparison - dowsntream analysis")),
            shiny::br(),
            shiny::h5("This tab has exploratory purposes. The user is invited to compare all the selected number of stable clusters by determining their identity with respect to the biological signal."),
            shiny::h5("\n"),
            shiny::h5("Besides functional features, such as identification of differentially expressed genes and enrichment analysis, the user can also generate paper-ready plots that describe the patterns of the clusters and set of genes."),
            shiny::h5("\n"),
            shiny::h5("All plots are available for download in multiple formats. We encourage the user to check the cog icon, where multiple plot settings can be adjusted."),
            shiny::h1("\n"),
            shiny::h5("For more information please go to:"),
            shiny::tagList("", github),
            easyClose = TRUE
        ),
        session
    )
}

compar_annotation_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Cluster annotation")),
            shiny::br(),
            shiny::h5("Here the user can create new metadata by annotating the clusters chosen from the previous tab."),
            shiny::h5("\n"),
            shiny::h5("This can be done by selecting any of the existing metadata columns and, based on the provided number of classes (the slider from right), the user can create a mapping between the new annotation and the groups of the selected metadata."),
            shiny::h5("\n"),
            shiny::h5("The annotation can be performed only when all the groups of the metadata are involved in the mapping. The user should also provide names for all the annotation groups, as well as for the name of the annotation."),
            easyClose = TRUE
        ),
        session
    )
}

compar_split_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Metadata splitting")),
            shiny::br(),
            shiny::h5("Here the user can create new metadata by combining the unique values from two discrete metadata columns."),
            shiny::h5("\n"),
            shiny::h5("This can be done by selecting any of the existing metadata columns. The two metadata columns must be different. By default, the name of the new column will be the concatenation of the two selected columns, but this behaviour can be changed by the user. The new column name shouldn't overlap with the name of an existing one."),
            shiny::h5("\n"),
            shiny::h5("Example: column A has values a, b, c; column B has values 1, 2, 3. After splitting, the new column A_B will have values a_1, b_2, c_3."),
            easyClose = TRUE
        ),
        session
    )
}

compar_distribution_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Gene and metadata distribution")),
            shiny::br(),
            shiny::h5("Similar with the 'Dimensionality reduction' tab, the user can interact with the metadata of the dataset and the expression level of the genes by visualising their distribution on the UMAP."),
            shiny::h5("\n"),
            shiny::h5("For gene distribution, the user can input one or multiple genes. In the latter case, a voting scheme is applied and the cells that express the genes above a specified threshold are highlighted with red. The threshold, as well as the number of genes where the threshold should be achieved simoultaneously, can be edited by the user"),
            shiny::h5("\n"),
            shiny::h5("The cells can also be included or discarded based on the existing metadata. Discarded cell groups will be coloured with gray."),
            easyClose = TRUE
        ),
        session
    )
}

compar_metadata_jsi_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Metadata JSI")),
            shiny::br(),
            shiny::h5("Here the user can visualise the relationship between the discrete metadata and the clustering results."),
            shiny::h5("\n"),
            shiny::h5("This can be useful in the context of identifying relationships between 'super-clusters' and more fine-grained metadata."),
            shiny::h5("\n"),
            shiny::h5("The displayed information can be either the number of common cells, or the Jaccard Similarity Index (JSI); the JSI is calculated row-wise."),
            easyClose = TRUE
        ),
        session
    )
}

compar_violin_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Violin distribution")),
            shiny::br(),
            shiny::h5("Here the user can distribution of any continuous metadata or of the expression of a gene as a violin plot. The user can change the visualisation to boxplots as well."),
            shiny::h5("\n"),
            shiny::h5("The information can be split using any discrete metadata available (including the chosen stable clusters). Also, the user can group the violins / boxplots by other metadata. Filtering on metadata groups is also available."),
            shiny::h5("\n"),
            shiny::h5("Below the violin plot, we provide a table with important stats related to the selected distribution."),
            easyClose = TRUE
        ),
        session
    )
}

compar_heatmap_gene_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Gene heatmap and bubbleplot")),
            shiny::br(),
            shiny::h5("The gene expression can be summarized at pseudobulk level using the discrete metadata. The summary is done using the average."),
            shiny::h5("\n"),
            shiny::h5("The resulting genes by metadata heatmap can be further scaled with respect to each gene. The values can be clipped to a threshold specified by the user."),
            shiny::h5("\n"),
            shiny::h5("The visualisation can be either a heatmap or a bubbleplot. For the bubbleplot, the size of the bubble is proportionate with the percentage of cells from the metadata group expressing the gene."),
            easyClose = TRUE
        ),
        session
    )
}

compar_markers_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Marker identification")),
            shiny::br(),
            shiny::h5("Using the Wilcoxon rank sum test, the user can identify differentially expressed genes by comparing two groups of cells. These group of cells are determined using metadata groups."),
            shiny::h5("\n"),
            shiny::h5("The candidate genes can be filtered using the average log2 fold change, the p-value, the average expression and the percentage of cells expressing the gene."),
            shiny::h5("\n"),
            shiny::h5("The analysis will output a table with the marker genes ordered by the average log2 fold change."),
            easyClose = TRUE
        ),
        session
    )
}

compar_enrichment_info <- function(session) {
    shiny::showModal(
        shiny::modalDialog(
            shiny::h3(shiny::strong("Enrichment analysis")),
            shiny::br(),
            shiny::h5("The identity of the clusters can be further investigated by performing an enrichment analysis on the marker genes."),
            shiny::h5("\n"),
            shiny::h5("The user can choose whether to use the positive or negative markers for the enrichment analysis. The analysis is compared by comparing the genes to the following databases: GO, KEGG, Reactome, TF and MIRNA."),
            shiny::h5("\n"),
            shiny::h5("The analysis will output a table with the enriched terms ordered by the adjusted p-value. The results can be visualised in an interactive dotplot."),
            easyClose = TRUE
        ),
        session
    )
}
