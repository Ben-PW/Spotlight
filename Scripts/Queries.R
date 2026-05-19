library(DBI)
library(duckdb)
library(here)
library(dplyr)
library(dbplyr)

con <- DBI::dbConnect(
  duckdb::duckdb(),
  dbdir = here::here("Results", "spotlight_results.duckdb")
)

on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

dbListTables(con)

############# Query top of databases to give a visual on structure ##############

node_results_gt <- DBI::dbGetQuery(con, "
                                SELECT *
                                FROM node_results_gt
                                ORDER BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level, NodeID
                                LIMIT 50;
                                ")

node_results <- DBI::dbGetQuery(con, "
                                SELECT *
                                FROM node_results
                                ORDER BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level, NodeID
                                LIMIT 50;
                                ")

network_results_gt <- DBI::dbGetQuery(con, "
                                      SELECT *
                                      FROM network_results_gt
                                      LIMIT 50;
                                      ")

network_results <- DBI::dbGetQuery(con, "
                                   SELECT *
                                   FROM network_results
                                   LIMIT 50;
                                   ")

################################ Add rank order columns ############################

# The original simulation did not calculate rank order of nodes by centrality scores,
# This will be useful for outcome metrics and is efficient to do using SQL

# Create table for rank order in ground truth networks

DBI::dbExecute(con, "
CREATE OR REPLACE TABLE node_results_GT_ranked AS
SELECT
  *,
  
  RANK() OVER (
    PARTITION BY dataset, replicate_id
    ORDER BY Degree_raw DESC
  ) AS Degree_raw_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id
    ORDER BY Degree_norm DESC
  ) AS Degree_norm_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id
    ORDER BY Betweenness_raw DESC
  ) AS Betweenness_raw_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id
    ORDER BY Betweenness_norm DESC
  ) AS Betweenness_norm_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id
    ORDER BY Closeness_raw DESC
  ) AS Closeness_raw_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id
    ORDER BY Closeness_norm DESC
  ) AS Closeness_norm_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id
    ORDER BY Eigenvector DESC
  ) AS Eigenvector_rank

FROM node_results_GT
")

# Create table for rank order in obsereved networks

DBI::dbExecute(con, "
CREATE OR REPLACE TABLE node_results_ranked AS
SELECT
  *,

  RANK() OVER (
    PARTITION BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level
    ORDER BY Degree_raw DESC
  ) AS Degree_raw_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level
    ORDER BY Degree_norm DESC
  ) AS Degree_norm_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level
    ORDER BY Betweenness_raw DESC
  ) AS Betweenness_raw_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level
    ORDER BY Betweenness_norm DESC
  ) AS Betweenness_norm_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level
    ORDER BY Closeness_raw DESC NULLS LAST
  ) AS Closeness_raw_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level
    ORDER BY Closeness_norm DESC NULLS LAST
  ) AS Closeness_norm_rank,

  RANK() OVER (
    PARTITION BY dataset, replicate_id, alpha, spotlight_pct, b, miss_level
    ORDER BY Eigenvector DESC
  ) AS Eigenvector_rank

FROM node_results
")

######################## Query to get network level differences ######################

# Because the tables for network level results are comparatively small, they can
# just be queried straight from the database. For larger runs, this would not be 
# feasible

network_bias_df <- tbl(con, "network_results") %>%
  inner_join(
    tbl(con, "network_results_gt"),
    by = c("dataset", "replicate_id"),
    suffix = c("_obs", "_gt")
  ) %>%
  mutate(
    density_ARB = (density_obs - density_gt)/density_gt,
    dcent_ARB = (dcent_obs - dcent_gt)/dcent_gt,
    clustering_ARB = (clustering_obs - clustering_gt)/clustering_gt,
    APL_ARB = (APL_obs - APL_gt)/APL_gt
  ) %>%
  collect()

# Query to check the spotlight simulation is properly biased towards degree

test <- DBI::dbGetQuery(con, "
  WITH spotlight_assignments AS (
    SELECT DISTINCT
      dataset,
      replicate_id,
      alpha,
      spotlight_pct,
      NodeID,
      Spotlight
    FROM node_results
  ),

  joined AS (
    SELECT
      s.dataset,
      s.replicate_id,
      s.alpha,
      s.spotlight_pct,
      s.NodeID,
      s.Spotlight,
      gt.Degree AS true_degree
    FROM spotlight_assignments s
    LEFT JOIN node_results_gt gt
      ON s.dataset = gt.dataset
     AND s.replicate_id = gt.replicate_id
     AND s.NodeID = gt.NodeID
  ),

  per_network AS (
    SELECT
      dataset,
      replicate_id,
      alpha,
      spotlight_pct,
      AVG(CASE WHEN Spotlight = 1 THEN true_degree END) AS mean_spotlit_degree,
      MAX(true_degree) AS max_degree_network
    FROM joined
    GROUP BY dataset, replicate_id, alpha, spotlight_pct
  )

  SELECT
    dataset,
    alpha,
    spotlight_pct,
    AVG(mean_spotlit_degree) AS mean_spotlit_degree,
    AVG(max_degree_network) AS mean_max_degree,
    AVG(mean_spotlit_degree / max_degree_network) AS relative_to_max
  FROM per_network
  GROUP BY dataset, alpha, spotlight_pct
  ORDER BY dataset, alpha, spotlight_pct;
")

############################################### BIAS CURVES #######################################

bias_curve_test_1 <- DBI::dbGetQuery(con, "
    SELECT
    obs.dataset,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,

    COUNT(*) AS n_graphs,

    AVG(obs.density - gt.density) AS mean_density_bias,
    AVG((obs.density - gt.density) / NULLIF(gt.density, 0)) AS mean_density_rel_bias,

    AVG(obs.dcent - gt.dcent) AS mean_dcent_bias,
    AVG((obs.dcent - gt.dcent) / NULLIF(gt.dcent, 0)) AS mean_dcent_rel_bias,

    AVG(obs.clustering - gt.clustering) AS mean_clustering_bias,
    AVG((obs.clustering - gt.clustering) / NULLIF(gt.clustering, 0)) AS mean_clustering_rel_bias,

    AVG(obs.APL - gt.APL) AS mean_APL_bias,
    AVG((obs.APL - gt.APL) / NULLIF(gt.APL, 0)) AS mean_APL_rel_bias

FROM network_results AS obs
JOIN network_results_GT AS gt
  ON obs.dataset = gt.dataset
 AND obs.replicate_id = gt.replicate_id

WHERE obs.source = 'observed'
  AND gt.source = 'true'

GROUP BY
    obs.dataset,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level

ORDER BY
    obs.dataset,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level;
                ")

################################## Node bias plots ##############################

# NB this query removes any nodes with infinite or NaN values, meaning these values 
# only hold for nodes connected in some way to a component

node_bias_df <- DBI::dbGetQuery(con, "
    SELECT
    obs.dataset,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,

    COUNT(*) AS n_nodes,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN (obs.Degree_raw - gt.Degree_raw) / gt.Degree_raw 
        END) 
        AS mean_degree_bias_spotlit,
        
    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN (obs.Degree_raw - gt.Degree_raw) / gt.Degree_raw 
        END) 
        AS mean_degree_bias_nonspotlit,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN (obs.Degree_raw - gt.Degree_raw) / gt.Degree_raw 
        END)
    -
    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN (obs.Degree_raw - gt.Degree_raw) / gt.Degree_raw 
        END) 
        AS degree_spotlight_lift,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN (obs.Betweenness_raw - gt.Betweenness_raw) / NULLIF(gt.Betweenness_raw, 0) 
        END)
    -
    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN (obs.Betweenness_raw - gt.Betweenness_raw) / NULLIF(gt.Betweenness_raw, 0) 
        END) 
        AS Betweenness_spotlight_lift,
        
    AVG(CASE 
          WHEN obs.Spotlight = 1
          AND NOT isnan(obs.Closeness_raw)
          AND NOT isnan(gt.Closeness_raw)
          AND gt.Closeness_raw != 0
          THEN (obs.Closeness_raw - gt.Closeness_raw) / gt.Closeness_raw
        END)
    -
    AVG(CASE 
          WHEN obs.Spotlight = 0
          AND NOT isnan(obs.Closeness_raw)
          AND NOT isnan(gt.Closeness_raw)
          AND gt.Closeness_raw != 0
          THEN (obs.Closeness_raw - gt.Closeness_raw) / gt.Closeness_raw
        END) 
        AS closeness_spotlight_lift,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN (obs.Eigenvector - gt.Eigenvector) / NULLIF(gt.Eigenvector, 0) 
        END)
    -
    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN (obs.Eigenvector - gt.Eigenvector) / NULLIF(gt.Eigenvector, 0) 
        END) 
        AS eigenvector_spotlight_lift

FROM node_results AS obs
JOIN node_results_GT AS gt
  ON obs.dataset = gt.dataset
 AND obs.replicate_id = gt.replicate_id
 AND obs.NodeID = gt.NodeID

WHERE obs.source = 'observed'
  AND gt.source = 'true'

GROUP BY
    obs.dataset,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level

ORDER BY
    obs.dataset,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level
")

# Problems with node bias plots
# turned out problem was some of the variables were normalised, some werent

node_bias_check <- DBI::dbGetQuery(con, "
                                   WITH node_bias AS (
  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    obs.Spotlight,

    (obs.Degree - gt.Degree) / NULLIF(gt.Degree, 0) AS degree_rb

  FROM node_results AS obs
  JOIN node_results_GT AS gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID

  WHERE obs.source = 'observed'
    AND gt.source = 'true'
),

graph_gap AS (
  SELECT
    dataset,
    replicate_id,
    alpha,
    spotlight_pct,
    b,
    miss_level,

    AVG(CASE WHEN Spotlight = 1 THEN degree_rb END) AS degree_rb_spotlit,
    AVG(CASE WHEN Spotlight = 0 THEN degree_rb END) AS degree_rb_nonspotlit,

    AVG(CASE WHEN Spotlight = 1 THEN degree_rb END)
    -
    AVG(CASE WHEN Spotlight = 0 THEN degree_rb END) AS degree_bias_gap

  FROM node_bias
  GROUP BY
    dataset,
    replicate_id,
    alpha,
    spotlight_pct,
    b,
    miss_level
)

SELECT *
FROM graph_gap
ORDER BY ABS(degree_bias_gap) ASC
LIMIT 50;")

############ Query for data for network level mean absolute bias plots #########

# These should be presented alongside average relative bias plots, to show the 
# difference between the magnitude of bias and the direction of bias

network_bias_long <- DBI::dbGetQuery(con, "
WITH gt_aug AS (
  SELECT
    *,
    CASE 
      WHEN dcent BETWEEN 0.05 AND 0.15 THEN 'Low'
      WHEN dcent BETWEEN 0.25 AND 0.35 THEN 'Med'
      WHEN dcent BETWEEN 0.45 AND 0.55 THEN 'High'
      ELSE 'WARNING'
    END AS 'baseline_centralisation'
  FROM network_results_gt
),

bias_long AS (

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,


    'density' AS metric,
    obs.density AS observed_value,
    gt.density AS gt_value,
    (obs.density - gt.density) / NULLIF(gt.density, 0) AS relative_bias,
    ABS((obs.density - gt.density) / NULLIF(gt.density, 0)) AS abs_relative_bias

  FROM network_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,


    'dcent' AS metric,
    obs.dcent AS observed_value,
    gt.dcent AS gt_value,
    (obs.dcent - gt.dcent) / NULLIF(gt.dcent, 0) AS relative_bias,
    ABS((obs.dcent - gt.dcent) / NULLIF(gt.dcent, 0)) AS abs_relative_bias

  FROM network_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,

    'clustering' AS metric,
    obs.clustering AS observed_value,
    gt.clustering AS gt_value,
    (obs.clustering - gt.clustering) / NULLIF(gt.clustering, 0) AS relative_bias,
    ABS((obs.clustering - gt.clustering) / NULLIF(gt.clustering, 0)) AS abs_relative_bias

  FROM network_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,


    'APL' AS metric,
    obs.APL AS observed_value,
    gt.APL AS gt_value,
    (obs.APL - gt.APL) / NULLIF(gt.APL, 0) AS relative_bias,
    ABS((obs.APL - gt.APL) / NULLIF(gt.APL, 0)) AS abs_relative_bias

  FROM network_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
)

SELECT *
FROM bias_long
WHERE abs_relative_bias IS NOT NULL
")

library(dplyr)
library(ggplot2)

metric_choice <- "dcent"

plot_df_1 <- network_bias_long %>%
  filter(metric == metric_choice) %>%
  filter(baseline_centralisation != "WARNING") %>%
  rename(gt_cent = baseline_centralisation) %>%
  mutate(
    gt_cent = factor(
      gt_cent,
      levels = c(
        "High",
        "Med",
        "Low"
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

################################# Trialing a model ###################################

model_df <- DBI::dbGetQuery(con, "
SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,

    net.density AS gt_density,
    net.dcent AS gt_centralisation,
    net.clustering AS gt_clustering,
    net.APL AS gt_APL,
    net.size AS gt_size,

    COUNT(*) AS n_nodes,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN obs.Degree - gt.Degree 
        END) AS mean_degree_bias_spotlit,

    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN obs.Degree - gt.Degree 
        END) AS mean_degree_bias_nonspotlit,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN obs.Degree - gt.Degree 
        END)
    -
    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN obs.Degree - gt.Degree 
        END) AS degree_spotlight_lift,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN obs.Betweenness - gt.Betweenness 
        END)
    -
    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN obs.Betweenness - gt.Betweenness 
        END) AS betweenness_spotlight_lift,

    AVG(CASE 
          WHEN obs.Spotlight = 1 
          THEN obs.Eigenvector - gt.Eigenvector 
        END)
    -
    AVG(CASE 
          WHEN obs.Spotlight = 0 
          THEN obs.Eigenvector - gt.Eigenvector 
        END) AS eigenvector_spotlight_lift

FROM node_results AS obs

JOIN node_results_GT AS gt
  ON obs.dataset = gt.dataset
 AND obs.replicate_id = gt.replicate_id
 AND obs.NodeID = gt.NodeID

LEFT JOIN network_results_GT AS net
  ON obs.dataset = net.dataset
 AND obs.replicate_id = net.replicate_id

WHERE obs.source = 'observed'
  AND gt.source = 'true'
  AND net.source = 'true'

GROUP BY
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    net.density,
    net.dcent,
    net.clustering,
    net.APL,
    net.size

ORDER BY
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level
")

m1 <- lm(
  degree_spotlight_lift ~ 
    miss_level * b +
    alpha +
    spotlight_pct +
    gt_centralisation +
    gt_density,
  data = model_df
)

m1
summary(m1)

m2 <- lm(
  degree_spotlight_lift ~ 
    miss_level * b +
    alpha * spotlight_pct +
    gt_centralisation * b +
    gt_density * miss_level,
  data = model_df
)

summary(m2)

model_df$base_graph_id <- interaction(
  model_df$dataset,
  model_df$replicate_id,
  drop = TRUE
)

m3 <- lme4::lmer(
  degree_spotlight_lift ~ 
    miss_level * b +
    alpha * spotlight_pct +
    gt_centralisation +
    gt_density +
    (1 | base_graph_id),
  data = model_df
)

summary(m3)

model_df$b_c <- model_df$b - 1
model_df$miss_c <- model_df$miss_level - mean(model_df$miss_level)
model_df$alpha_c <- model_df$alpha - mean(model_df$alpha)
model_df$spotlight_pct_c <- model_df$spotlight_pct - mean(model_df$spotlight_pct)

m3c <- lme4::lmer(
  degree_spotlight_lift ~ 
    miss_c * b_c +
    alpha * spotlight_pct +
    gt_centralisation +
    gt_density +
    (1 | base_graph_id),
  data = model_df
)

summary(m3c)
plot(m3c)
qqnorm(resid(m3c)); qqline(resid(m3c))

model_df$.resid <- resid(m3c)

model_df |>
  dplyr::arrange(abs(.resid) |> desc()) |>
  dplyr::select(
    dataset, replicate_id, alpha, spotlight_pct, b, miss_level,
    gt_centralisation, gt_density, degree_spotlight_lift, .resid
  ) |>
  head(20)

library(brms)

library(brms)

model_df <- model_df |>
  dplyr::mutate(
    base_graph_id = interaction(dataset, replicate_id, drop = TRUE),
    alpha_c = alpha - mean(alpha),
    spotlight_pct_c = spotlight_pct - mean(spotlight_pct),
    miss_c = miss_level - mean(miss_level),
    b_c = b - 1
  )

model_df$b_f <- factor(model_df$b)

model_df$b_f <- factor(model_df$b)

m_absurd <- brm(
  bf(
    degree_spotlight_lift ~
      miss_c * b_c * alpha_c * spotlight_pct_c *
      gt_centralisation * gt_density +
      
      s(gt_centralisation, k = 3) +
      s(gt_density, k = 3) +
      s(miss_level, by = b_f, k = 3) +
      
      (1 | dataset) +
      (1 + miss_c * b_c | base_graph_id),
    
    sigma ~ dataset + miss_level + spotlight_pct
  ),
  
  data = model_df,
  family = student(),
  chains = 4,
  cores = 3,
  iter = 4000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)
