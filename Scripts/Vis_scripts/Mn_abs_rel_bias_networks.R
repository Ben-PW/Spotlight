################################################################################

# This script is for mean absolute relative bias plots for network level metrics
# This outcome metric was chosen to visualise the absoliute level of bias, which 
# could otherwise be masked if networks were biased in different directions.

################################################################################

# Input dataframe: mn_abs_rel_bias_nets

############################ Degree centralisation #############################

metric_choice <- "dcent"

plot_df_1 <- mn_abs_rel_bias_nets %>%
  filter(metric == metric_choice) %>%
  filter(baseline_centralisation != "WARNING") %>%
  rename(gt_cent = baseline_centralisation) %>%
  mutate(
    gt_cent = factor(
      gt_cent,
      levels = c(
        "Low",
        "Med",
        "High"
      )
    )
  ) %>%
  group_by(
    miss_level,
    b,
    alpha,
    gt_cent
  ) %>%
  summarise(
    mean_abs_rb = mean(abs_relative_bias, na.rm = TRUE),
    q25 = quantile(abs_relative_bias, 0.25, na.rm = TRUE),
    q75 = quantile(abs_relative_bias, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  plot_df_1,
  aes(
    x = miss_level,
    y = mean_abs_rb,
    colour = factor(b),
    fill = factor(b),
    group = factor(b)
  )
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
  labs(
    x = "Missingness level",
    y = "Mean absolute relative bias",
    colour = "Non-spotlit\ntie weight (b)",
    fill = "Non-spotlit\ntie weight (b)",
    title = paste("Absolute relative bias in", metric_choice),
    subtitle = "Faceted by baseline centralisation and spotlight selection bias"
  ) +
  theme_minimal()

########################################## APL #################################

metric_choice <- "APL"

plot_df_1 <- network_bias_long %>%
  filter(metric == metric_choice) %>%
  filter(baseline_centralisation != "WARNING") %>%
  rename(gt_cent = baseline_centralisation) %>%
  mutate(
    gt_cent = factor(
      gt_cent,
      levels = c(
        "Low",
        "Med",
        "High"
      )
    )
  ) %>%
  group_by(
    miss_level,
    b,
    alpha,
    gt_cent
  ) %>%
  summarise(
    mean_abs_rb = mean(abs_relative_bias, na.rm = TRUE),
    q25 = quantile(abs_relative_bias, 0.25, na.rm = TRUE),
    q75 = quantile(abs_relative_bias, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  plot_df_1,
  aes(
    x = miss_level,
    y = mean_abs_rb,
    colour = factor(b),
    fill = factor(b),
    group = factor(b)
  )
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
  labs(
    x = "Missingness level",
    y = "Mean absolute relative bias",
    colour = "Non-spotlit\ntie weight (b)",
    fill = "Non-spotlit\ntie weight (b)",
    title = paste("Absolute relative bias in", metric_choice),
    subtitle = "Faceted by baseline centralisation and spotlight selection bias"
  ) +
  theme_minimal()

####################################### Clustering ##############################

metric_choice <- "clustering"

plot_df_1 <- network_bias_long %>%
  filter(metric == metric_choice) %>%
  filter(baseline_centralisation != "WARNING") %>%
  rename(gt_cent = baseline_centralisation) %>%
  mutate(
    gt_cent = factor(
      gt_cent,
      levels = c(
        "Low",
        "Med",
        "High"
      )
    )
  ) %>%
  group_by(
    miss_level,
    b,
    alpha,
    gt_cent
  ) %>%
  summarise(
    mean_abs_rb = mean(abs_relative_bias, na.rm = TRUE),
    q25 = quantile(abs_relative_bias, 0.25, na.rm = TRUE),
    q75 = quantile(abs_relative_bias, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  plot_df_1,
  aes(
    x = miss_level,
    y = mean_abs_rb,
    colour = factor(b),
    fill = factor(b),
    group = factor(b)
  )
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
  labs(
    x = "Missingness level",
    y = "Mean absolute relative bias",
    colour = "Non-spotlit\ntie weight (b)",
    fill = "Non-spotlit\ntie weight (b)",
    title = paste("Absolute relative bias in", metric_choice),
    subtitle = "Faceted by baseline centralisation and spotlight selection bias"
  ) +
  theme_minimal()
