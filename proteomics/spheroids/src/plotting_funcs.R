#' Plotting utilities for spheroid proteomics analyses.
#'
#' This file contains helper functions for visualizing proteomics data from the
#' spheroid analysis workflow. The functions in this script cover PCA plotting,
#' missing-value and outlier summaries, differential-expression visualization,
#' and layout helpers for combining multiple figures into a single annotated
#' output. Most functions operate on `SummarizedExperiment` objects and use
#' metadata in `colData()` to label samples and groups.
#'
#' @docType package
#' @name plotting_funcs
#' @title Plotting helpers for spheroid proteomics workflows
#' @keywords internal
#'
#' @description
#' The plotting utilities in this file provide a consistent interface for
#' exploratory data checks and downstream analysis summaries. They are designed
#' to support PCA exploration, missingness diagnostics, sample outlier
#' detection, and differential-expression visualization for both global
#' proteomics and phosphoproteomics data.

library(tibble)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(rstatix)

color_vals_grey_first <- c(
  '#BBBBBB',
  '#EE7733',
  '#0077BB',
  '#33BBEE',
  '#EE3377',
  '#CC3311',
  '#009988'
)
color_vals_black_first <- c(
  '#000000',
  '#EE7733',
  '#0077BB',
  '#33BBEE',
  '#EE3377',
  '#CC3311',
  '#009988'
)
color_vals_default <- c(
  '#EE7733',
  '#0077BB',
  '#33BBEE',
  '#EE3377',
  '#CC3311',
  '#009988',
  '#BBBBBB'
)


#' Generate a PCA plot for a proteomics experiment
#'
#' Creates a principal component analysis plot from a `SummarizedExperiment`
#' object, optionally subsetting the data before plotting and mapping aesthetics
#' such as color, shape, and labels from column metadata.
#'
#' @param experiment A `SummarizedExperiment` object containing assay values and
#'   sample metadata in `colData()`.
#' @param subset_by An optional column name in `colData(experiment)` used to
#'   subset samples before PCA.
#' @param subset_values A vector of values to keep when `subset_by` is supplied.
#' @param ... Additional optional plotting arguments. Supported values include:
#'   \itemize{
#'     \item{`color_by`}{Column name in `colData(experiment)` used to color
#'       points in the PCA plot.}
#'     \item{`shape_by`}{Column name in `colData(experiment)` used to set point
#'       shape.}
#'     \item{`label_by`}{Column name in `colData(experiment)` used to label
#'       points with text or text-repulsion.}
#'     \item{`color_vals`}{Optional named or unnamed vector of colors supplied
#'       to `scale_color_manual()`. If omitted, default palette values are used.}
#'     \item{`title`}{Custom plot title to override the default
#'       `"Spheroid proteomics PCA"`.}
#'     \item{`subtitle`}{Optional subtitle displayed below the title.}
#'     \item{`legend_position`}{Position for the legend passed to
#'       `theme(legend.position = ...)`. Common values include `"bottom"`,
#'       `"right"`, and `"none"`.}
#'   }
#' @returns A `ggplot` object showing the first two principal components.
#'
#' @examples
#' \dontrun{
#' pca_plot <- plot_pca(
#'   experiment = se,
#'   subset_by = "agent",
#'   subset_values = c("DMSO", "Trab"),
#'   color_by = "agent",
#'   shape_by = "timepoint",
#'   label_by = "condition",
#'   title = "PCA by agent",
#'   subtitle = "Global proteomics"
#' )
#' }
plot_pca <- function(
  experiment,
  subset_by = NULL,
  subset_values = NULL,
  ...
) {
  # Get the additional arguments passed to the function and store them in a list.
  extra_args <- list(...)

  # Check if the user has defined whether to subset the experiment object (and if so to which value).
  # If so, subset the experiment object accordingly.
  # Otherwise, use the entire experiment object.
  if (!is.null(subset_by) & !is.null(subset_values)) {
    exp <- experiment[, colData(experiment)[[subset_by]] %in% subset_values]
  } else {
    exp <- experiment
  }

  # Run the PCA analysis on the transposed assay data (prcomp expects samples as rows) of the experiment object.
  # Then combine the results with the column data of the experiment object (i.e. the metadata).
  principal_components <- prcomp(t(assay(exp))) %$%
    cbind(x, colData(exp))

  # Setup the aesthetics for the PCA plot based on the user-defined parameters.
  aes_args <- list(x = rlang::sym("PC1"), y = rlang::sym("PC2"))
  if (!is.null(extra_args$color_by)) {
    aes_args$colour <- rlang::sym(extra_args$color_by)
  }
  if (!is.null(extra_args$shape_by)) {
    aes_args$shape <- rlang::sym(extra_args$shape_by)
  }
  if (!is.null(extra_args$label_by)) {
    aes_args$label <- rlang::sym(extra_args$label_by)
  }

  # Setting up the color values for the PCA plot based on the user-defined parameters.
  # If the user has defined a variable to color by check if the user has defined custom color values.
  # If not, create a named vector of color values based on the unique values of the variable to color by.
  # If the user has not defined a variable to color by, set the color values to NULL resulting in no color mapping in the PCA plot.
  if (!is.null(extra_args$color_by)) {
    if (!is.null(extra_args$color_vals)) {
      colors <- extra_args$color_vals
    } else {
      colors <- color_vals_default[
        1:length(unique(principal_components[[extra_args$color_by]]))
      ]
      setNames(
        colors,
        unique(principal_components[[extra_args$color_by]])
      )
    }
  } else {
    colors <- NULL
  }

  # Create the PCA plot using ggplot2 with the defined aesthetics.
  plot <- (ggplot(
    principal_components,
    do.call(aes, aes_args)
  ) +
    geom_point(size = 2) +
    theme_minimal() +
    theme(legend.position = "bottom"))

  # Add additional layers to the PCA plot based on the user-defined parameters.
  if (!is.null(extra_args$label_by)) {
    plot <- plot +
      geom_text_repel(max.overlaps = 10, size = 2, show.legend = FALSE)
  }

  # Add custom shape and color scales to the PCA plot based on the user-defined parameters.
  if (!is.null(extra_args$shape_by)) {
    plot <- plot +
      scale_shape_manual(
        values = 1:length(unique(principal_components[[extra_args$shape_by]]))
      )
  }
  if (!is.null(extra_args$color_by)) {
    plot <- plot + scale_color_manual(values = colors)
  }

  # Add title and subtitle to the PCA plot based on the user-defined parameters.
  # If the user has not defined a title or subtitle, use default values.
  if (!is.null(extra_args$title)) {
    plot_title <- extra_args$title
  } else {
    plot_title <- "Spheroid proteomics PCA"
  }
  if (!is.null(extra_args$subtitle)) {
    plot_subtitle <- extra_args$subtitle
  } else {
    plot_subtitle <- NULL
  }
  plot <- plot +
    labs(
      title = plot_title,
      subtitle = plot_subtitle
    )

  return(plot)
}

#' Plot PCA by agent (deprecated)
#'
#' This function is deprecated and should no longer be used. Prefer
#' `plot_pca()` for PCA visualizations because it supports more flexible
#' subsetting and metadata-driven aesthetics through a consistent `...` argument
#' interface.
#'
#' @description
#' Creates a PCA scatter plot for samples belonging to a specific agent, while
#' also highlighting the corresponding control condition and any associated plate
#' data. The plot is built from the first two principal components of the assay
#' matrix and annotates samples by metadata columns such as `agent`,
#' `timepoint`, and `condition`.
#'
#' @param experiment A `SummarizedExperiment` object containing assay values and
#'   sample metadata in `colData()`.
#' @param agent Character scalar naming the agent whose samples should be plotted.
#' @param type Character scalar specifying the assay type. Must be either
#'   `"global"` or `"phospho"`.
#' @param missingness Character string describing the missingness filter applied
#'   before plotting. This is included in the subtitle for recordkeeping.
#' @param imputation Character string describing the imputation strategy used
#'   before plotting. This is included in the subtitle for recordkeeping.
#'
#' @return A `ggplot` object showing the first two principal components for the
#'   selected agent, with points colored by agent and shaped by timepoint.
#'
#' @details This helper builds a subset of the experiment containing the target
#' agent and its control condition, then computes PCA on the assay matrix. It is
#' retained only for backward compatibility and is not recommended for new code.
#'
#' @examples
#' \dontrun{
#' pca_plot <- plot_pca_by_agent(
#'   experiment = se,
#'   agent = "Trab",
#'   type = "global",
#'   missingness = "50%",
#'   imputation = "SampMin"
#' )
#' }
plot_pca_by_agent <- function(
  experiment,
  agent,
  type = c('global', 'phospho'),
  missingness = 'None',
  imputation = "None"
) {
  prot_type <- match.arg(type)

  plates <- unique(colData(experiment[, experiment$agent == agent])$plate)
  plates <- append(plates, 10)
  ctrl <- ifelse(test = (agent == 'Trab'), yes = "Water", no = "DMSO")
  sel <- c(ctrl, agent)

  experiment_subset <- experiment[,
    (experiment$plate %in% plates) & experiment$agent %in% sel
  ]

  principal_components <- prcomp(t(assay(experiment_subset))) %$%
    cbind(x, colData(experiment_subset))

  plot_title <- base::paste("Spheroid", prot_type, "proteomics for:", agent)
  plot_subtitle <- base::paste(
    "Missingness filter:",
    missingness,
    "; Imputation:",
    imputation
  )
  color_vals <- c("#66c2a5", "#fc8d62", "#8da0cb")
  colors <- setNames(color_vals, c(agent, "DMSO", "Water"))

  plot <- (ggplot(
    principal_components,
    aes(
      x = PC1,
      y = PC2,
      # col = timepoint,
      # shape = agent,
      col = agent,
      shape = timepoint,
      label = condition
    )
  ) +
    geom_point(size = 2) +
    # geom_text_repel(max.overlaps = 10, size = 2) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = colors) +
    # + scale_color_manual(values = c("DMSO" = "#e66101", "Mirda" = "#5e3c99", "NA" = "grey"))
    labs(
      title = plot_title,
      subtitle = "50% missingness filter; SampMin imputation"
    ))
  return(plot)
}

#' Plot PCA by plate (deprecated)
#'
#' This function is deprecated and should no longer be used. Prefer
#' `plot_pca()` for PCA visualizations because it supports more flexible
#' metadata-based customizations and a unified `...` interface.
#'
#' @description
#' Creates a PCA scatter plot for a single plate by extracting the samples from
#' that plate and visualizing the first two principal components. Points are
#' colored by agent and labeled with the sample condition metadata, which makes
#' it useful for quick plate-level exploratory checks.
#'
#' @param experiment A `SummarizedExperiment` object containing assay values and
#'   sample metadata in `colData()`.
#' @param plate Character or numeric scalar identifying the plate to subset.
#' @param type Character scalar specifying the assay type. Must be either
#'   `"global"` or `"phospho"`.
#' @param missingness Character string describing the missingness filter applied
#'   before plotting. This is stored in the subtitle for provenance.
#' @param imputation Character string describing the imputation strategy used
#'   before plotting. This is included in the subtitle for provenance.
#' @param ... Additional optional plotting arguments. Supported values include
#'   `title`, `subtitle`, `colors`, and `label`.
#'
#' @return A `ggplot` object showing the first two principal components for the
#'   selected plate.
#'
#' @details This helper subsets the experiment to one plate, runs PCA on the
#' assay matrix, and then overlays agent-based coloring and condition labels.
#' It is retained only for backward compatibility and is not recommended for new
#' code.
#'
#' @examples
#' \dontrun{
#' pca_plot <- plot_pca_by_plate(
#'   experiment = se,
#'   plate = 10,
#'   type = "global",
#'   missingness = "50%",
#'   imputation = "SampMin",
#'   title = "PCA for Plate 10",
#'   label = "condition"
#' )
#' }
plot_pca_by_plate <- function(
  experiment,
  plate,
  type = c('global', 'phospho'),
  missingness = 'None',
  imputation = "None",
  ...
) {
  prot_type <- match.arg(type)

  experiment_subset <- experiment[, experiment$plate == plate]

  principal_components <- prcomp(t(assay(experiment_subset))) %$%
    cbind(x, colData(experiment_subset))

  extra_args <- list(...)
  if (!is.null(extra_args$title)) {
    plot_title <- extra_args$title
  } else {
    plot_title <- base::paste0(
      "Spheroid ",
      prot_type,
      " proteomics for Plate: ",
      plate
    )
  }

  if (!is.null(extra_args$subtitle)) {
    plot_subtitle <- extra_args$subtitle
  } else {
    plot_subtitle <- base::paste0(
      "Missingness filter: ",
      missingness,
      "; Imputation: ",
      imputation
    )
  }

  if (!is.null(extra_args$colors)) {
    colors <- extra_args$colors
  } else {
    agents <- unique(colData(experiment_subset)$agent)
    ctrl <- ifelse(test = ("DMSO" %in% agents), yes = "DMSO", no = "Water")
    agents_sel <- agents[agents != ctrl]
    colors <- setNames(color_vals_black_first, append(c(ctrl), agents_sel))
  }

  if (!is.null(extra_args$label)) {
    label_ <- extra_args$label
  } else {
    label_ <- 'condition'
  }

  plot <- (ggplot(
    principal_components,
    aes(
      x = PC1,
      y = PC2,
      # col = timepoint,
      shape = agent,
      col = agent,
      # shape = timepoint,
      label = get(label_)
    )
  ) +
    geom_point(size = 2.5) +
    # geom_text_repel(max.overlaps = 10, size = 3, show.legend = F) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = colors) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle
    ))
  return(plot)
}

#' Identify outlier values using the IQR rule
#'
#' @description
#' Detects values in a numeric vector that fall outside the usual interquartile
#' range bounds. This helper is useful for flagging obvious outliers in summary
#' plots or exploratory checks before downstream filtering or visualization.
#'
#' @param x A numeric vector to evaluate for outliers.
#'
#' @return A logical vector of the same length as `x`, where `TRUE` indicates
#'   that the value lies outside the lower or upper IQR threshold.
#'
#' @details The function uses the standard Tukey-style rule: any value below
#' `Q1 - 1.5 * IQR` or above `Q3 + 1.5 * IQR` is considered an outlier.
#'
#' @examples
#' \dontrun{
#' x <- c(1, 2, 3, 4, 100)
#' is_outlier(x)
#' }
is_outlier <- function(x) {
  return(
    x < quantile(x, 0.25) - 1.5 * IQR(x) |
      x > quantile(x, 0.75) + 1.5 * IQR(x)
  )
}


#' Plot missing-value overview by plate
#'
#' @description
#' Summarizes the number of missing values per sample and visualizes the
#' distribution of missingness across plates. The plot highlights samples that
#' appear to be outliers based on the IQR rule and is useful for identifying
#' low-quality or irregular samples before downstream analysis.
#'
#' @param experiment A `SummarizedExperiment` object containing assay values and
#'   sample metadata in `colData()`.
#' @param type Character scalar specifying the assay type. Must be either
#'   `"global"` or `"phospho"`.
#'
#' @return A `ggplot` object showing a boxplot of the number of missing values
#'   per plate, with any outlier samples labeled by condition.
#'
#' @details The function calculates the total number of missing values for each
#'   sample, joins this summary to the sample metadata, and generates a plate-
#'   level boxplot. Samples flagged by `is_outlier()` are labeled with their
#'   condition name.
#'
#' @examples
#' \dontrun{
#' na_plot <- plot_per_plate_na_overview(
#'   experiment = se,
#'   type = "global"
#' )
#' }
plot_per_plate_na_overview <- function(
  experiment,
  type = c('global', 'phospho')
) {
  prot_type <- match.arg(type)

  summarized_data <- as_tibble(assay(experiment)) %>%
    summarise(across(everything(), ~ sum(is.na(.)))) %>%
    pivot_longer(cols = everything(), names_to = "sample_id", values_to = "NAs")

  metadata <- colData(experiment) %>%
    as.data.frame(.) %>%
    rownames_to_column(., var = "sample_id") %>%
    as_tibble(.)

  data_plot <- inner_join(metadata, summarized_data, by = 'sample_id') %>%
    group_by(plate) %>%
    mutate(outlier = ifelse(is_outlier(NAs), condition, NA)) %>%
    mutate(
      .,
      agent = replace_when(agent, agent %in% c("DMSO", "Water") ~ "DMSO/Water")
    )

  plot_title <- base::paste("Outliers in spheroid", prot_type, "proteomics")
  plot_subtitle <- base::paste(
    "# of NAs per sample out of",
    dim(assay(experiment))[1],
    ifelse(prot_type == "global", "proteins", "phospho-sites")
  )

  plot <- (ggplot(
    data_plot,
    aes(
      x = plate,
      y = NAs,
      group = plate,
      xmin = 1,
      xmax = 10
    )
  ) +
    geom_boxplot() +
    geom_label_repel(
      aes(label = outlier),
      na.rm = TRUE,
      max.overlaps = 10,
      size = 2,
      min.segment.length = 0
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # + scale_x_continuous(breaks = seq(1, 10, 1))
    labs(
      title = plot_title,
      subtitle = plot_subtitle
    ))

  return(plot)
}

#' Plot differential-expression volcano plot for a single agent and timepoint.
#'
#' @description
#' Builds a volcano-style visualization for a pairwise comparison between an
#' experimental agent and its control group. The function computes a differential
#' expression table for the requested comparison and plots log fold-change against
#' negative log10-adjusted p-value, highlighting features that pass the chosen
#' thresholds.
#'
#' @param experiment A `SummarizedExperiment` object containing proteomics assay
#'   values and sample metadata in `colData()`.
#' @param agent Character scalar naming the treatment agent to compare against its
#'   control.
#' @param timepoint Character scalar naming the timepoint to compare.
#' @param log_fc Numeric scalar defining the absolute log2 fold-change cutoff for
#'   highlighting differential features. Defaults to `1`.
#' @param pval Numeric scalar specifying the adjusted p-value threshold for
#'   significance. Defaults to `0.05`.
#' @param prot_type Character scalar specifying the assay type. Must be either
#'   `"global"` or `"phospho"`.
#' @param spheroid Optional character scalar naming the spheroid label when a
#'   nested group structure is present.
#'
#' @return A `ggplot` object showing the differential-expression summary for the
#'   selected comparison, with vertical and horizontal threshold lines and color-
#'   coded expression changes.
#'
#' @details The helper derives the appropriate control condition from `agent`,
#' calls `calc_diff_ex()` for the requested comparison, and then formats the
#' result with `extract_diff_ex_features()`. The resulting plot emphasizes genes
#' or phosphosites whose log fold-change and adjusted p-value exceed the selected
#' thresholds.
#'
#' @examples
#' \dontrun{
#' p <- plot_diff_ex(
#'   experiment = se,
#'   agent = "Trab",
#'   timepoint = "24hr",
#'   log_fc = 1,
#'   pval = 0.05,
#'   prot_type = "global"
#' )
#' }
plot_diff_ex <- function(
  experiment,
  agent,
  timepoint,
  log_fc = 1,
  pval = 0.05,
  prot_type = c('global', 'phospho'),
  spheroid = NULL
) {
  prot_type <- match.arg(prot_type)
  align_by_group <- ifelse(
    test = (prot_type == 'global'),
    yes = "Protein.Group",
    no = "Protein.With.Phospho.Site"
  )

  ctrl <- ifelse(test = (agent == 'Trab'), yes = "Water", no = "DMSO")
  if (!is.null(spheroid)) {
    group_1 <- paste0(spheroid, "_", timepoint, "_", ctrl)
    group_2 <- paste0(spheroid, "_", timepoint, "_", agent)
  } else {
    group_1 <- paste0(timepoint, "_", ctrl)
    group_2 <- paste0(timepoint, "_", agent)
  }

  diff_ex <- calc_diff_ex(
    experiment,
    group_1,
    group_2,
    align_by = align_by_group
  )

  if (prot_type == "global") {
    dfPlot <- extract_diff_ex_features(diff_ex, "Genes", log_fc, pval)
  } else if (prot_type == "phospho") {
    dfPlot <- extract_diff_ex_features(
      diff_ex,
      "Gene.With.Phospho.Site",
      log_fc,
      pval
    )
  } else {
    stop(paste0("Invalid prot_type selected: ", prot_type))
  }

  # plot_title = paste('Spheroid', prot_type ,'proteomics:', agent, "vs.", ctrl, "at", timepoint)
  plot_subtitle <- base::paste0(
    timepoint,
    "; ",
    "#down: ",
    dim(dfPlot[dfPlot$expression_change == 'down', ])[1],
    " / ",
    "#up: ",
    dim(dfPlot[dfPlot$expression_change == 'up', ])[1]
  )

  plot <- (ggplot(
    data = dfPlot,
    aes(
      x = logFC,
      y = -log10(adj.P.Val),
      col = expression_change,
      label = labels
    )
  ) +
    geom_vline(
      xintercept = c(-log_fc, log_fc),
      col = "gray",
      linetype = 'dashed'
    ) +
    geom_hline(yintercept = -log10(pval), col = "gray", linetype = 'dashed') +
    geom_point(size = .75) +
    # geom_text_repel(max.overlaps = 5, size = 2) +
    scale_color_manual(
      values = c(
        'up' = '#74add1',
        'down' = '#f46d43',
        'not significant' = 'grey'
      )
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      # title = plot_title,
      subtitle = plot_subtitle,
      col = 'expression',
      y = expression("-log"[10] * "p-Value")
    ))

  plot
}

#' Arrange multiple differential-expression plots in a grid
#'
#' @description
#' Combines a list of ggplot objects into a single multi-panel layout and adds an
#' optional figure title above the grid. Useful to combine multiple
#' differential-expression volcano plots for different timepoints or agents
#' into a single figure.
#'
#' @param plots A list of `ggplot` objects to arrange into a common grid.
#' @param plot_title Optional character string used as the title for the final
#'   combined figure. If `NULL`, no title is added.
#'
#' @return A `ggplot`/`ggpubr` figure object created by `ggpubr::annotate_figure()`
#'   containing the arranged plots and an optional title.
#'
#' @details The function first lays out the plots with `ggpubr::ggarrange()`
#' using three columns and a shared legend. It then adds a top title using
#' `ggpubr::text_grob()`, returning the final combined figure object for display
#' or export.
#'
#' @examples
#' \dontrun{
#' p1 <- plot_diff_ex(experiment = se, agent = "Trab", timepoint = "8hr")
#' p2 <- plot_diff_ex(experiment = se, agent = "Trab", timepoint = "24hr")
#' p3 <- plot_diff_ex(experiment = se, agent = "Trab", timepoint = "48hr")
#'
#' combined_plot <- plot_diff_ex_grid(
#'   plots = list(p1, p2, p3),
#'   plot_title = "Differential expression summary"
#' )
#' }
plot_diff_ex_grid <- function(plots, plot_title = NULL) {
  plot_grid <- ggpubr::ggarrange(
    plotlist = plots,
    ncol = 3,
    common.legend = TRUE,
    legend = "bottom"
  )
  plot_annotated <- ggpubr::annotate_figure(
    plot_grid,
    top = ggpubr::text_grob(plot_title, size = 14)
  )

  plot_annotated
}


plot_features_across_timepoints <- function(
  experiment,
  features,
  agent,
  timepoints = c('8hr', '24hr', '48hr'),
  ftype = 'Gene.Names',
  fmapping = 'Gene.With.Phospho.Site'
) {
  fexp <- experiment[rowData(experiment)[, ftype] %in% features, ]

  if (nrow(fexp) == 0) {
    print(paste(paste(features, collapse = ','), 'not found in this data'))
    return(NULL)
  }

  rd <- rowData(experiment) |>
    as.data.frame() |>
    tibble::rownames_to_column('feature') |>
    dplyr::rename(feature_name = fmapping) |>
    dplyr::select(feature, feature_name)

  ##now we can get those features across timepoints and agents
  full_dat <- do.call(
    rbind,
    lapply(timepoints, function(tp) {
      do.call(
        rbind,
        lapply(agent, function(ag) {
          ##annoying but we have to do this for every time point
          te <- exp_subset(fexp, ag, tp)

          if (nrow(te) == 0) {
            return(data.frame())
          }
          te <- te |>
            impute_proteomics() #|>
          #log_transform()

          df <- assay(te) |>
            as.data.frame() |>
            tibble::rownames_to_column('feature') |>
            left_join(rd) |>
            tidyr::pivot_longer(
              -c(feature, feature_name),
              names_to = 'sample',
              values_to = 'expression'
            ) |>
            left_join(tibble::rownames_to_column(
              as.data.frame(colData(te)),
              'sample'
            ))
        })
      )
    })
  )

  ##now calculate p-values
  #first filter for complete values

  complete <- full_dat |>
    drop_na() |>
    group_by(timepoint, feature_name) |>
    summarize(ngroups = n_distinct(agent))

  if (length(agent) == 1) {
    complete <- complete |>
      subset(ngroups == 2) |>
      ungroup()

    pdat <- full_dat |>
      inner_join(complete) |> ###only filters for those that are complete
      group_by(timepoint, feature_name) |>
      summarize(pval = t.test(expression ~ agent)$p.value) |>
      mutate(pval = sprintf("p = %.3f", pval))
  } else {
    pdat <- full_dat |>
      dplyr::select('agent', 'timepoint', 'feature_name') |>
      distinct() |>
      mutate(pval = "")
  }

  ninf <- which(is.finite(full_dat$expression))

  full_dat <- full_dat[ninf, ]

  #now do the plot
  full_dat |>
    left_join(pdat) |>
    rowwise() |>
    mutate(newgroup = paste(feature_name, agent, timepoint, sep = '_')) |>
    ggplot(aes(
      x = feature_name,
      y = expression,
      color = agent,
      fill = agent,
      label = pval
    )) +
    geom_boxplot() +
    #    geom_jitter()+
    geom_text(
      data = pdat,
      color = 'black',
      aes(
        label = pval,
        x = feature_name,
        y = max(full_dat$expression, na.rm = T) * .9
      ),
      size = 2
    ) +
    facet_grid(timepoint ~ .) +
    scale_fill_manual(
      values = c(
        "#74ADDc",
        "#BFC3A4",
        "#F8BC8D",
        '#f46d43',
        '#a35667',
        '#a4ba5a',
        '#eee110',
        '#0088cf',
        '#aa55bb',
        '#a3aabb',
        '#a4ee55',
        '#88aabb'
      )
    ) +
    scale_color_manual(
      values = c(
        "#74ADDc",
        "#BFC3A4",
        "#F8BC8D",
        '#f46d43',
        '#a35667',
        '#a4ba5a',
        '#eee110',
        '#0088Cf',
        '#aa55bb',
        '#a3aabb',
        '#a4ee55',
        '#88aabb'
      )
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_enrich_result <- function(
  all.kin,
  expdatalist, #list of differential expression results
  enrichlist, #list of enrichment
  agent,
  fmapping = 'Gene.With.Phospho.Site'
) {
  tps <- all.kin |>
    arrange(rev(zscore))

  forder <- unique(tps$feature)
  #get max z for color
  maxz <- max(abs(all.kin$zscore))

  p1 <- ggplot(
    all.kin,
    aes(y = feature, x = timepoint, color = zscore, size = -log10(pval))
  ) +
    geom_point() +
    ggplot2::scale_y_discrete(limits = forder) +
    ggplot2::scale_x_discrete(limits = unique(all.kin$timepoint)) +
    scale_color_gradient2(
      low = "#2266bb",
      high = "#cc3311",
      mid = "lightgrey",
      limits = c(-maxz, maxz)
    ) +
    theme_minimal(base_size = 12) +
    geom_point(
      data = subset(all.kin, sig),
      col = "black",
      stroke = 1,
      shape = 21
    ) +
    ggtitle(paste(agent, 'Enrichment'))

  ##which kinases show up at all three timepoints?

  tops <- tps |>
    subset(sig) |>
    group_by(feature) |>
    summarize(tps = n()) |>
    subset(tps == length(timepoints))

  print(paste('found', nrow(tops), 'that are enriched at all time points'))
  names(expdatalist) <- unique(all.kin$timepoint)
  names(enrichlist) <- unique(all.kin$timepoint)
  ##now we plot the substrates in distribution compared to other features
  sub.tab <- do.call(
    rbind,
    lapply(timepoints, function(tp) {
      exp <- expdatalist[[tp]]

      ##get all features
      all.feats <- rowData(exp)[, fmapping]
      ##now we have, for each kinase/pathway, the features that drive it
      el <- enrichlist[[tp]] |>
        subset(feature %in% tops$feature) |>
        dplyr::select(feature, ingroupnames) |>
        tidyr::separate_longer_delim(ingroupnames, delim = ', ') |>
        subset(ingroupnames %in% all.feats) |>
        mutate(type = 'in group')

      ##now get off targets
      nt <- do.call(
        rbind,
        lapply(tops$feature, function(path) {
          targs <- el |>
            subset(feature == path)
          offtargs <- setdiff(all.feats, targs$ingroupnames)
          data.frame(
            feature = path,
            ingroupnames = offtargs,
            type = 'out group'
          )
        })
      )

      ##now we have full table with annotations
      full_df <- rbind(el, nt) |>
        mutate(timepoint = tp) |>
        left_join(dplyr::select(
          as.data.frame(rowData(exp)),
          ingroupnames = fmapping,
          'logFC'
        ))

      return(full_df)
    })
  )
  ##now we can plot the substrates

  if (nrow(sub.tab) > 0) {
    p2 <- ggplot(sub.tab, aes(x = logFC, fill = type, alpha = 0.5)) +
      geom_density() +
      facet_grid(feature ~ timepoint) +
      scale_fill_manual(values = c("#2266bb", "#cc3311")) +
      theme_minimal(base_size = 12) +
      ggtitle(paste(agent, 'Enrichment drivers'))
  } else {
    p2 <- NULL
  }
  return(c(p1, p2))
}
