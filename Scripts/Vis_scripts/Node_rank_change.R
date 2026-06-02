################################################################################

# This is the script for visualising node level rank change due to spotlight 
# effects, without replying on common outcome metrics such as TopN

################################################################################

make_rank_change_plot_df <- function(df,
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

plot_rank_change_outcome <- function(plot_df,
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
      subtitle = "Rank-change outcomes, summarised across graph-condition replicates"
    ) +
    theme_minimal()
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  p
}

rank_change_outcomes <- c(
  "mean_norm_rank_lift",
  "mean_abs_norm_rank_change",
  "mean_norm_rank_lift_spotlit",
  "mean_norm_rank_lift_nonspotlit",
  "mean_abs_norm_rank_change_spotlit",
  "mean_abs_norm_rank_change_nonspotlit",
  "spotlight_rank_lift_gap",
  "spotlight_abs_rank_change_gap"
)

rank_change_labels <- c(
  mean_norm_rank_lift = "Mean normalised rank lift",
  mean_abs_norm_rank_change = "Mean absolute normalised rank change",
  mean_norm_rank_lift_spotlit = "Mean normalised rank lift: spotlit nodes",
  mean_norm_rank_lift_nonspotlit = "Mean normalised rank lift: non-spotlit nodes",
  mean_abs_norm_rank_change_spotlit = "Mean absolute normalised rank change: spotlit nodes",
  mean_abs_norm_rank_change_nonspotlit = "Mean absolute normalised rank change: non-spotlit nodes",
  spotlight_rank_lift_gap = "Spotlight rank lift gap",
  spotlight_abs_rank_change_gap = "Spotlight absolute rank-change gap"
)

rank_change_metrics <- c(
  "Degree",
  "Betweenness",
  "Closeness",
  "Eigenvector"
)

rank_change_plots <- list()

for (m in rank_change_metrics) {
  for (outcome in rank_change_outcomes) {
    
    plot_df <- make_rank_change_plot_df(
      df = rank_lift_df,
      metric_choice = m,
      outcome_col = outcome,
      spotlight_pct_choice = 0.05
    )
    
    y_limits <- NULL
    reference_line <- NULL
    
    # Directional rank lift and gap metrics are centred around 0
    if (outcome %in% c(
      "mean_norm_rank_lift",
      "mean_norm_rank_lift_spotlit",
      "mean_norm_rank_lift_nonspotlit",
      "spotlight_rank_lift_gap"
    )) {
      reference_line <- 0
    }
    
    # Absolute rank-change metrics cannot go below 0
    if (outcome %in% c(
      "mean_abs_norm_rank_change",
      "mean_abs_norm_rank_change_spotlit",
      "mean_abs_norm_rank_change_nonspotlit"
    )) {
      y_limits <- c(0, NA)
    }
    
    # Absolute gap can be positive or negative, so centre at 0
    if (outcome == "spotlight_abs_rank_change_gap") {
      reference_line <- 0
    }
    
    plot_name <- paste(m, outcome, sep = "_")
    
    rank_change_plots[[plot_name]] <- plot_rank_change_outcome(
      plot_df = plot_df,
      metric_choice = m,
      outcome_label = rank_change_labels[[outcome]],
      y_limits = y_limits,
      reference_line = reference_line
    )
  }
}

rank_change_plots$Degree_spotlight_rank_lift_gap
rank_change_plots$Betweenness_spotlight_rank_lift_gap
rank_change_plots$Closeness_spotlight_rank_lift_gap
rank_change_plots$Eigenvector_spotlight_rank_lift_gap

rank_change_plots$Degree_mean_norm_rank_lift_spotlit
rank_change_plots$Betweenness_mean_norm_rank_lift_spotlit
rank_change_plots$Closeness_mean_norm_rank_lift_spotlit
rank_change_plots$Eigenvector_mean_norm_rank_lift_spotlit
