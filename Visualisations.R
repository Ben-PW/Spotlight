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
