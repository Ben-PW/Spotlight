################################################################################

# This script is for visualisatino of node level mean absolute and relative bias,
# however these values are aggregated over all nodes, distinctions are not made 
# between spotlit and non-spotlit nodes

################################################################################

node_bias_agg1 <- node_bias_agg1 %>%
  mutate(
    baseline_centralisation = factor(
      baseline_centralisation,
      levels = c("High", "Med", "Low")
    ),
    b = factor(b),
    alpha = factor(alpha)
  )

plot_df <- node_bias_summary %>%
  filter(metric == metric_choice) %>%
  group_by(
    metric,
    baseline_centralisation,
    alpha,
    b,
    miss_level
  ) %>%
  summarise(
    mean_relative_bias = mean(mean_relative_bias, na.rm = TRUE),
    rb_q25 = quantile(mean_relative_bias, 0.25, na.rm = TRUE),
    rb_q75 = quantile(mean_relative_bias, 0.75, na.rm = TRUE),
    
    mean_abs_relative_bias = mean(mean_abs_relative_bias, na.rm = TRUE),
    abs_rb_q25 = quantile(mean_abs_relative_bias, 0.25, na.rm = TRUE),
    abs_rb_q75 = quantile(mean_abs_relative_bias, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

################################################################################

############ Plotting function ##############

plot_node_mean_relative_bias <- function(df, metric_choice) {
  
  plot_df <- df %>%
    filter(metric == metric_choice)
  
  ggplot(
    plot_df,
    aes(
      x = miss_level,
      y = mean_relative_bias,
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
      aes(
        ymin = rb_q25,
        ymax = rb_q75
      ),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(
      baseline_centralisation ~ alpha,
      labeller = label_both
    ) +
    labs(
      x = "Missingness level",
      y = "Mean relative bias",
      colour = "Non-spotlit\ntie weight (b)",
      fill = "Non-spotlit\ntie weight (b)",
      title = paste("Mean relative bias in", metric_choice),
      subtitle = "Node-level centrality bias, averaged within graphs then across conditions"
    ) +
    theme_minimal()
}

################# Plots ################

plot_node_mean_relative_bias(node_bias_agg1, "Degree_raw")

##################################### Thoughts ####################################

# The above didn't work, the query didn't aggregate over spotlight pct, meaning there were
# large numbers of lines. The data will need to be re-queried and then post-processed
# outside of SQL, which is fine as long as it's aggregated

node_bias_agg2 %>%
  count(metric, baseline_centralisation, alpha, spotlight_pct, b, miss_level) %>%
  summarise(
    min_graphs = min(n),
    max_graphs = max(n)
  )

node_bias_agg2 %>%
  group_by(metric) %>%
  summarise(
    mean_nodes_total = mean(n_nodes_total, na.rm = TRUE),
    mean_nodes_used = mean(n_nodes_used, na.rm = TRUE),
    mean_nodes_dropped = mean(n_nodes_dropped, na.rm = TRUE),
    .groups = "drop"
  )

##################################### Attempt 2 #####################################

# Current thoughts here is that spotligh tpct isn't doing a huge amount here, at 
# at least for degree

metric_choice <- "Closeness_raw"
spotlight_pct_choice <- 0.05


plot1df <- node_bias_agg2 %>%
  filter(
    metric == metric_choice,
    spotlight_pct == spotlight_pct_choice
  ) %>%
  group_by(
    metric,
    baseline_centralisation,
    alpha,
    b,
    miss_level
  ) %>%
  summarise(
    mean_relative_bias = mean(graph_mean_relative_bias, na.rm = TRUE),
    rb_q25 = quantile(graph_mean_relative_bias, 0.25, na.rm = TRUE),
    rb_q75 = quantile(graph_mean_relative_bias, 0.75, na.rm = TRUE),
    
    mean_abs_relative_bias = mean(graph_mean_abs_relative_bias, na.rm = TRUE),
    abs_rb_q25 = quantile(graph_mean_abs_relative_bias, 0.25, na.rm = TRUE),
    abs_rb_q75 = quantile(graph_mean_abs_relative_bias, 0.75, na.rm = TRUE),
    
    mean_nodes_dropped = mean(n_nodes_dropped, na.rm = TRUE),
    n_graphs = n(),
    .groups = "drop"
  ) %>%
  mutate(
    baseline_centralisation = factor(
      baseline_centralisation,
      levels = c("High", "Med", "Low")
    ),
    alpha = factor(alpha),
    b = factor(b)
  )

# Plot mean relative binas

mean_rel <- ggplot(
  plot1df,
  aes(
    x = miss_level,
    y = mean_relative_bias,
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
    aes(ymin = rb_q25, ymax = rb_q75),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(
    baseline_centralisation ~ alpha,
    labeller = label_both
  ) +
  labs(
    x = "Missingness level",
    y = "Mean relative bias",
    colour = "Non-spotlit\ntie weight (b)",
    fill = "Non-spotlit\ntie weight (b)",
    title = paste("Mean relative bias in", metric_choice),
    subtitle = "Node bias averaged within graph, then summarised across graphs"
  ) +
  theme_minimal()

mean_rel

mean_abs <- ggplot(
  plot1df,
  aes(
    x = miss_level,
    y = mean_abs_relative_bias,
    colour = b,
    fill = b,
    group = b
  )
) +
  geom_ribbon(
    aes(ymin = abs_rb_q25, ymax = abs_rb_q75),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(
    baseline_centralisation ~ alpha,
    labeller = label_both
  ) +
  labs(
    x = "Missingness level",
    y = "Mean absolute relative bias",
    colour = "Non-spotlit\ntie weight (b)",
    fill = "Non-spotlit\ntie weight (b)",
    title = paste("Mean absolute relative bias in", metric_choice),
    subtitle = "Node bias averaged within graph, then summarised across graphs"
  ) +
  theme_minimal()

mean_abs
