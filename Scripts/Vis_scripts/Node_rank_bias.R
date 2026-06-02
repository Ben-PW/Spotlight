################################################################################

# This is the script for visualising the various node level outcome metrics which
# distinguish between spotlit and non spotlit nodes

make_topN_plot_df <- function(df,
                              metric_choice,
                              outcome_col,
                              spotlight_pct_choice = 0.05) {
  
  df %>%
    filter(
      metric == metric_choice,
      spotlight_pct == spotlight_pct_choice
    ) %>%
    mutate(
      gt_cent = case_when(
        str_detect(dataset, "_c1$") ~ "Low",
        str_detect(dataset, "_c3$") ~ "Med",
        str_detect(dataset, "_c5$") ~ "High",
        TRUE ~ "WARNING"
      )
    ) %>%
    filter(gt_cent != "WARNING") %>%
    group_by(
      gt_cent,
      alpha,
      b,
      miss_level
    ) %>%
    summarise(
      mean_value = mean(.data[[outcome_col]], na.rm = TRUE),
      q25 = quantile(.data[[outcome_col]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[outcome_col]], 0.75, na.rm = TRUE),
      n_graphs = sum(!is.na(.data[[outcome_col]])),
      .groups = "drop"
    ) %>%
    mutate(
      gt_cent = factor(gt_cent, levels = c("High", "Med", "Low")),
      alpha = factor(alpha),
      b = factor(b)
    )
}

plot_topN_outcome <- function(plot_df,
                              metric_choice,
                              outcome_label,
                              y_limits = NULL,
                              reference_line = NULL) {
  
  p <- ggplot(
    plot_df,
    aes(
      x = miss_level,
      y = mean_value,
      colour = b,
      fill = b,
      group = b
    )
  )
  
  if (!is.null(reference_line)) {
    p <- p +
      geom_hline(
        yintercept = reference_line,
        linetype = "dashed",
        linewidth = 0.4
      )
  }
  
  p <- p +
    geom_ribbon(
      aes(ymin = q25, ymax = q75),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(
      gt_cent ~ alpha,
      labeller = label_both
    ) +
    labs(
      x = "Missingness level",
      y = outcome_label,
      colour = "Non-spotlit\ntie weight (b)",
      fill = "Non-spotlit\ntie weight (b)",
      title = paste(metric_choice, "-", outcome_label),
      subtitle = "Top 10% node-set outcomes, summarised across graph-condition replicates"
    ) +
    theme_minimal()
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  p
}

topN_outcomes <- c(
  "precision",
  "recall",
  "jaccard_overlap",
  "spotlight_lift_obs_top",
  "spotlight_lift_gt_top",
  "excess_spotlight_lift"
)

topN_labels <- c(
  precision = "Top-10 precision",
  recall = "Top-10 recall",
  jaccard_overlap = "Top-10 Jaccard overlap",
  spotlight_lift_obs_top = "Observed Top-10 spotlight lift",
  spotlight_lift_gt_top = "GT Top-10 spotlight lift",
  excess_spotlight_lift = "Excess observed spotlight lift"
)

topN_metrics <- c(
  "Degree",
  "Betweenness",
  "Closeness",
  "Eigenvector"
)

topN_plots <- list()

for (m in topN_metrics) {
  for (outcome in topN_outcomes) {
    
    plot_df <- make_topN_plot_df(
      df = node_rank_df,
      metric_choice = m,
      outcome_col = outcome,
      spotlight_pct_choice = 0.05
    )
    
    # Sensible y-axis and reference lines depending on outcome type
    y_limits <- NULL
    reference_line <- NULL
    
    if (outcome %in% c("precision", "recall", "jaccard_overlap")) {
      y_limits <- c(0, 1)
    }
    
    if (outcome %in% c("spotlight_lift_obs_top", "spotlight_lift_gt_top")) {
      reference_line <- 1
    }
    
    if (outcome == "excess_spotlight_lift") {
      reference_line <- 0
    }
    
    plot_name <- paste(m, outcome, sep = "_")
    
    topN_plots[[plot_name]] <- plot_topN_outcome(
      plot_df = plot_df,
      metric_choice = m,
      outcome_label = topN_labels[[outcome]],
      y_limits = y_limits,
      reference_line = reference_line
    )
  }
}

topN_plots$Degree_excess_spotlight_lift
topN_plots$Betweenness_excess_spotlight_lift
topN_plots$Closeness_excess_spotlight_lift
topN_plots$Eigenvector_excess_spotlight_lift

topN_plots$Degree_jaccard_overlap
topN_plots$Betweenness_jaccard_overlap
topN_plots$Closeness_jaccard_overlap
topN_plots$Eigenvector_jaccard_overlap

