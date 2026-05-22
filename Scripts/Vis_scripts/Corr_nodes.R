################################################################################

# This script is for visualisations of AGGREGATED node level centrality correlations, both
# rank and pearson. Any spotlight specific visualisations should go in the dedicated 
# script which doesn't exist yet but fingers crossed will by the end of the day.

################################################################################

make_node_corr_plot_df <- function(df,
                                   corr_col,
                                   spotlight_pct_choice = 0.10) {
  
  df %>%
    filter(spotlight_pct == spotlight_pct_choice) %>%
    mutate(
      gt_cent = case_when(
        str_detect(dataset, "_c1$") ~ "Low",
        str_detect(dataset, "_c3$") ~ "Med",
        str_detect(dataset, "_c5$") ~ "High"
      )
    ) %>%
    group_by(
      gt_cent,
      alpha,
      b,
      miss_level
    ) %>%
    summarise(
      mean_corr = mean(.data[[corr_col]], na.rm = TRUE),
      q25 = quantile(.data[[corr_col]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[corr_col]], 0.75, na.rm = TRUE),
      n_graphs = n(),
      .groups = "drop"
    ) %>%
    mutate(
      gt_cent = factor(gt_cent, levels = c("High", "Med", "Low")),
      alpha = factor(alpha),
      b = factor(b)
    )
}

plot_node_corr <- function(plot_df, corr_label) {
  
  ggplot(
    plot_df,
    aes(
      x = miss_level,
      y = mean_corr,
      colour = b,
      fill = b,
      group = b
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
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
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      x = "Missingness level",
      y = "Mean correlation",
      colour = "Non-spotlit\ntie weight (b)",
      fill = "Non-spotlit\ntie weight (b)",
      title = corr_label,
      subtitle = "Correlations calculated within graph, then summarised across graphs"
    ) +
    theme_minimal()
}

degree_rank_plot_df <- make_node_corr_plot_df(
  df = node_corr_df,
  corr_col = "degree_rank_corr",
  spotlight_pct_choice = 0.05
)

plot_node_corr(
  plot_df = degree_rank_plot_df,
  corr_label = "Degree centrality rank correlation"
)

corr_cols <- c(
  "degree_rank_corr",
  "closeness_rank_corr",
  "betweenness_rank_corr",
  "eigenvector_rank_corr",
  "degree_corr",
  "closeness_corr",
  "betweenness_corr",
  "eigenvector_corr"
)

corr_labels <- c(
  degree_rank_corr = "Degree centrality rank correlation",
  closeness_rank_corr = "Closeness centrality rank correlation",
  betweenness_rank_corr = "Betweenness centrality rank correlation",
  eigenvector_rank_corr = "Eigenvector centrality rank correlation",
  degree_corr = "Degree centrality value correlation",
  closeness_corr = "Closeness centrality value correlation",
  betweenness_corr = "Betweenness centrality value correlation",
  eigenvector_corr = "Eigenvector centrality value correlation"
)

node_corr_plots <- lapply(corr_cols, function(corr_col) {
  
  plot_df <- make_node_corr_plot_df(
    df = node_corr_df,
    corr_col = corr_col,
    spotlight_pct_choice = 0.10
  )
  
  plot_node_corr(
    plot_df = plot_df,
    corr_label = corr_labels[[corr_col]]
  )
})

names(node_corr_plots) <- corr_cols

node_corr_plots$degree_rank_corr
node_corr_plots$closeness_rank_corr
node_corr_plots$betweenness_rank_corr
node_corr_plots$eigenvector_rank_corr
node_corr_plots$degree_corr
node_corr_plots$closeness_corr
node_corr_plots$betweenness_corr
node_corr_plots$eigenvector_corr
