library(ggplot2)
library(dplyr)

here::here()

############################## Network level metrics #############################

# Mean absolute relative bais in network level metrics. Variables visualised are:
# Target metrics: 
# - clustering
# - degree centralisation
# - APL
# Design variables:
# - alpha
# - missingness level
# - spotlight strength
# - ground truth centralisation

# Required dataframe: mn_abs_rel_bias_nets
source(here::here("Scripts", "Vis_scripts", "Mn_abs_rel_bias_networks.R"))

# Average relative bias in network level metrics. Variables visualised are:
# Target metrics: 
# - clustering
# - degree centralisation
# - APL
# Design variables:
# - alpha
# - missingness level
# - spotlight strength
# - ground truth centralisation

# Required dataframe: network_bias_df
source(here::here("Scripts", "Vis_scripts", "Mn_rel_bias_networks.R"))

############################### Node level metrics ##############################

############## Aggregating over all nodes (no spitlight distinction) ############

############## IMPORTANT
# Below will not work, relative bias requires dividing by the ground truth value,
# which means any value of 0 will result in a comparison being discarded
# Node level plots will have to rely on correlations and the other planned metrics
# Mean absolute relative bias and mean relative bias in node centrality metrics
# induced by spotlight effects. Variables visualised are:
# Target metrics:
# - Degree centrality
# - Betweenness centrality
# - Closeness centrality
# - Eigenvector ccentrality
# Design variables:
# - alpha
# - missingness level
# - spotlight strength
# - ground truth centralisation

# Required dataframe: node_bias_agg1
#source(here::here("Scripts", "Vis_scripts", "Mn_bias_nodes_agg.R"))

# Correlation plots between node level values 
# Target metrics:
# - Degree centrality
# - Betweenness centrality
# - Closeness centrality
# - Eigenvector centrality
# Design variables
# - alpha
# - missingness level
# - spotlight strength
# - ground truth centralisation

# Required dataframe: node_corr_df
library(stringr)
source(here::here("Scripts", "Vis_scripts", "Corr_nodes.R"))

############################# Non aggregated node metrics ########################

# Plots of node level outcomes such as TopN, OverlapN etc
# Target metrics:
# - Degree centrality
# - Betweenness centrality
# - Closeness centrality
# - Eigenvector centrality
# Design variables
# - alpha
# - missingness level
# - spotlight strength
# - ground truth centralisation

# Required dataframe: node_rank_df
source(here::here("Scripts", "Vis_scripts", "Node_rank_bias.R"))

# Plots of node level rank change statistics
# Target metrics:
# - Degree centrality rank change
# - Betweenness centrality rank change
# - Closeness centrality rank change
# - Eigenvector centrality rank change
# Design variables
# - alpha
# - missingness level
# - spotlight strength
# - ground truth centralisation

# Required dataframe: rank_lift_df
source(here::here("Scripts", "Vis_scripts", "Node_rank_change.R"))


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

# APL
ggplot(network_bias_df, aes(x = miss_level_obs, y = APL_ARB, colour = factor(b_obs))) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.15, aes(fill = factor(b_obs)), colour = NA) +
  facet_wrap(~ alpha_obs)

# Above plots but binned by centralisation

# Thoughts about these plots
# Design variables displayed:
# - Spotlight degree bias
# - Missingness level
# - Centralisation
# - Spotlight sampling bias
# Design variables omitted:
# - Density
# - Size
# - Spotlight sampling fraction

# APL
network_bias_df %>%
  dplyr::mutate(
    dcent_gt_bin = dplyr::ntile(dcent_gt, 3),
    dcent_gt_bin = factor(
      dcent_gt_bin,
      labels = c("Low GT centralisation",
                 "Mid GT centralisation",
                 "High GT centralisation")
    )
  ) %>%
  ggplot(aes(
    x = miss_level_obs,
    y = dcent_ARB,
    colour = factor(b_obs),
    fill = factor(b_obs)
  )) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(
    fun.data = mean_cl_boot,
    geom = "ribbon",
    alpha = 0.15,
    colour = NA
  ) +
  facet_grid(dcent_gt_bin ~ alpha_obs) +
  labs(
    x = "Missingness level",
    y = "Average Relative Bias in dcent",
    colour = "b",
    fill = "b"
  )


# APL
network_bias_df %>%
  dplyr::mutate(
    dcent_gt_bin = dplyr::ntile(dcent_gt, 3),
    dcent_gt_bin = factor(
      dcent_gt_bin,
      labels = c("Low GT centralisation",
                 "Mid GT centralisation",
                 "High GT centralisation")
    )
  ) %>%
  ggplot(aes(
    x = miss_level_obs,
    y = dcent_ARB,
    colour = factor(b_obs),
    fill = factor(b_obs)
  )) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(
    fun.data = mean_cl_boot,
    geom = "ribbon",
    alpha = 0.15,
    colour = NA
  ) +
  facet_grid(dcent_gt_bin ~ alpha_obs) +
  labs(
    x = "Missingness level",
    y = "Average relative bias in degree centralisation",
    colour = "b",
    fill = "b"
  )

# Clustering

network_bias_df %>%
  dplyr::mutate(
    dcent_gt_bin = dplyr::ntile(dcent_gt, 3),
    dcent_gt_bin = factor(
      dcent_gt_bin,
      labels = c("Low GT centralisation",
                 "Mid GT centralisation",
                 "High GT centralisation")
    )
  ) %>%
  ggplot(aes(
    x = miss_level_obs,
    y = clustering_ARB,
    colour = factor(b_obs),
    fill = factor(b_obs)
  )) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(
    fun.data = mean_cl_boot,
    geom = "ribbon",
    alpha = 0.15,
    colour = NA
  ) +
  facet_grid(dcent_gt_bin ~ alpha_obs) +
  labs(
    x = "Missingness level",
    y = "Average relative bias in global clustering coefficient",
    colour = "b",
    fill = "b"
  )

################################# Heat map attempt (traumatic) ##########################

ggplot(network_bias_df,
       aes(x = factor(spotlight_pct_obs), y = factor(b_obs), fill = dcent_ARB)) +
  geom_tile() +
  facet_grid(miss_level_obs ~ alpha_obs) +
  labs(
    x = "Spotlight proportion",
    y = "Non-spotlit tie deletion weight",
    fill = "Mean relative bias"
  )

library(viridis)

library(viridis)

ggplot(
  network_bias_df,
  aes(
    x = factor(spotlight_pct_obs),
    y = factor(b_obs),
    fill = dcent_ARB
  )
) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(-1, 0.5),
    oob = scales::squish,
    name = "Mean relative bias"
  ) +
  facet_grid(miss_level_obs ~ alpha_obs) +
  labs(
    x = "Spotlight proportion",
    y = "Non-spotlit tie deletion weight"
  ) +
  theme_minimal()

ggplot(network_bias_df,
       aes(x = factor(spotlight_pct_obs), y = factor(b_obs), fill = clustering_ARB)) +
  geom_tile() +
  facet_grid(miss_level_obs ~ alpha_obs) +
  labs(
    x = "Spotlight proportion",
    y = "Non-spotlit tie deletion weight",
    fill = "Mean relative bias"
  )

network_bias_df %>%
  dplyr::mutate(
    dcent_gt_bin = cut(dcent_gt, breaks = 3)
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
  geom_smooth(se = TRUE) + 
  facet_grid(miss_level_obs ~ b_obs) +
  labs(
    x = "Ground-truth degree centralisation",
    y = "Degree centralisation relative bias"
  )

ggplot(network_bias_df, aes(x = dcent_gt, y = dcent_ARB)) +
  geom_point(alpha = 0.2) +
  geom_smooth(se = TRUE) +
  facet_grid(miss_level_obs ~ b_obs)


