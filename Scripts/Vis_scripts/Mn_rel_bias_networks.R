################################################################################

# This script is for visualisations of mean relative bias in network level
# metrics, which provides an indication of the direction of bias, even though 
# the scale of the bias can be slightly masked

################################################################################

############################ Degree centralisation #############################

network_bias_df <- network_bias_df %>%
  dplyr::mutate(
    dcent_gt_bin = dplyr::ntile(dcent_gt, 3),
    dcent_gt_bin = factor(
      dcent_gt_bin,
      labels = c("Low dcent",
                 "Med dcent",
                 "High dcent")
    )
  ) 

network_bias_df %>%
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

network_bias_df %>%
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
    y = "Average Relative Bias in clustering",
    colour = "b",
    fill = "b"
  )

network_bias_df %>%
  ggplot(aes(
    x = miss_level_obs,
    y = APL_ARB,
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
    y = "Average Relative Bias in Av Path Length",
    colour = "b",
    fill = "b"
  )
