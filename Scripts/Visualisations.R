########################### Bias slopes test ####################################

library(ggplot2)

# This needs to be faceted by dataset as well
ggplot(
  bias_curve_test_1,
  aes(
    x = miss_level,
    y = mean_dcent_rel_bias,
    colour = b,
    group = b
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(dataset + alpha ~ spotlight_pct) +
  labs(
    x = "Missingness level",
    y = "Mean degree centralisation bias",
    colour = "Non-spotlit deletion weight (b)",
    title = "Bias in degree centralisation under spotlighted missingness"
  ) +
  theme_minimal()

ggplot(
  bias_curve_test_1,
  aes(
    x = miss_level,
    y = mean_dcent_rel_bias,
    colour = factor(b),
    group = factor(b)
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(dataset + alpha ~ spotlight_pct) +
  labs(
    x = "Missingness level",
    y = "Mean degree centralisation bias",
    colour = "b"
  ) +
  theme_minimal()

################################ Simple network bias plots #########################

# Density (obvious due to being a variable)
ggplot(network_bias_df, aes(x = miss_level_obs, y = density_ARB, colour = factor(b_obs))) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.15, aes(fill = factor(b_obs)), colour = NA) +
  facet_wrap(~ alpha_obs)

# Degree centralistion (this one is decent)
ggplot(network_bias_df, aes(x = miss_level_obs, y = dcent_ARB, colour = factor(b_obs))) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.15, aes(fill = factor(b_obs)), colour = NA) +
  facet_wrap(~ alpha_obs)

# Clustering
ggplot(network_bias_df, aes(x = miss_level_obs, y = clustering_ARB, colour = factor(b_obs))) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.15, aes(fill = factor(b_obs)), colour = NA) +
  facet_wrap(~ alpha_obs)

################################# Heat map attempt (traumatic) ##########################

ggplot(network_bias_df,
       aes(x = factor(spotlight_pct_obs), y = factor(b_obs), fill = mean_dcent_ARB)) +
  geom_tile() +
  facet_grid(miss_level_obs ~ alpha_obs) +
  labs(
    x = "Spotlight proportion",
    y = "Non-spotlit tie deletion weight",
    fill = "Mean relative bias"
  )

network_bias_df %>%
  dplyr::mutate(
    dcent_gt_bin = cut(dcent_gt, breaks = 2)
  ) %>%
  ggplot(aes(x = dcent_gt_bin, y = dcent_ARB)) +
  geom_boxplot() +
  facet_grid(miss_level_obs ~ b_obs) +
  labs(
    x = "Ground-truth degree centralisation",
    y = "Degree centralisation relative bias"
  )
ggplot(network_bias_df, aes(x = dcent_gt, y = dcent_ARB)) +
  geom_point(alpha = 0.2) +
  #geom_smooth(se = TRUE) +
  facet_grid(miss_level_obs ~ b_obs)

############################### Node lift plots ##################################

library(ggplot2)
library(dplyr)

lift_df <- lift_df %>%
  mutate(
    b = factor(b),
    alpha = factor(alpha),
    spotlight_pct = factor(spotlight_pct)
  )

# Degree lift

ggplot(
  lift_df,
  aes(
    x = b,
    y = degree_spotlight_lift,
    group = interaction(miss_level, dataset),
    colour = factor(miss_level)
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_grid(dataset+ alpha ~ spotlight_pct) +
  labs(
    x = "Spotlight strength (b)",
    y = "Degree spotlight lift",
    colour = "Missingness",
    title = "Spotlight-induced inflation of node degree"
  ) +
  theme_minimal()

# Betweenness lift
# Normalised betweenness will be required for the full run, this plot is 
# illogical

ggplot(
  lift_df,
  aes(
    x = b,
    y = betweenness_spotlight_lift,
    group = interaction(miss_level, dataset),
    colour = factor(miss_level)
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_grid(dataset+ alpha ~ spotlight_pct) +
  labs(
    x = "Spotlight strength (b)",
    y = "Betweenness spotlight lift",
    colour = "Missingness",
    title = "Spotlight-induced inflation of node betweenness"
  ) +
  theme_minimal()

ggplot(
  lift_df,
  aes(
    x = b,
    y = eigenvector_spotlight_lift,
    group = interaction(miss_level, dataset),
    colour = factor(miss_level)
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_grid(dataset+ alpha ~ spotlight_pct) +
  labs(
    x = "Spotlight strength (b)",
    y = "Betweenness spotlight lift",
    colour = "Missingness",
    title = "Spotlight-induced inflation of node betweenness"
  ) +
  theme_minimal()
