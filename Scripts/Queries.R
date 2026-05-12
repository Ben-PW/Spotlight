library(DBI)
library(duckdb)
library(here)

con <- dbConnect(
  duckdb(),
  dbdir = here("Results", "spotlight_results.duckdb")
)

on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

dbListTables(con)

# Query to check the spotlight simulation is properly biased towards degree

DBI::dbGetQuery(con, "
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

################################## Node lift plots ##############################

lift_df <- DBI::dbGetQuery(con, "
    SELECT
    obs.dataset,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,

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
