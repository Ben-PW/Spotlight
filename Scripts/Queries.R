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

node_results_ranked <- DBI::dbGetQuery(con, "
                                       SELECT *
                                       FROM node_results_ranked
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


#################################### NETWORK BIAS ###################################

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


################################## Node bias plots ##############################

################################# Correlation plots #############################

# This dataframe is to calculate the pearson and spearman correlations between the
# node level centrality values of the data. After testing it turns out relative 
# bias metrics aren't suitable

# Bit of a gross query but I find creating the temporary table first easier to wrap
# my head around

node_corr_df <- DBI::dbGetQuery(con, "

-- Creating a temporary table temp of all the necessary variables

    WITH temp AS (
      SELECT
      
      -- listing obs_rank vars out explicitly for clarity instead of using SELECT *
      
        obs_rank.dataset,
        obs_rank.replicate_id,
        -- obs_rank.graph_id, removed as probably useless
        obs_rank.alpha,
        obs_rank.spotlight_pct,
        obs_rank.b,
        obs_rank.miss_level,
        
        -- observed rank centralisation
        
        obs_rank.Degree_raw_rank AS obs_degree_rank,
        obs_rank.Betweenness_raw_rank AS obs_betweenness_rank,
        obs_rank.Closeness_raw_rank AS obs_closeness_rank,
        obs_rank.Eigenvector_rank AS obs_eigenvector_rank,
      
        -- gt rank centralisation
        
        gt_rank.Degree_raw_rank AS gt_degree_rank,
        gt_rank.Betweenness_raw_rank AS gt_betweenness_rank,
        gt_rank.Closeness_raw_rank AS gt_closeness_rank,
        gt_rank.Eigenvector_rank AS gt_eigenvector_rank,
        
        -- obsereved centralisation
        
        obs_rank.Degree_raw AS obs_degree,
        obs_rank.Closeness_raw AS obs_closeness,
        obs_rank.Betweenness_raw AS obs_betweenness,
        obs_rank.Eigenvector AS obs_eigenvector,
        
        -- gt centralisation
        
        gt_rank.Degree_raw AS gt_degree,
        gt_rank.Closeness_raw AS gt_closeness,
        gt_rank.Betweenness_raw AS gt_betweenness,
        gt_rank.Eigenvector AS gt_eigenvector,
      
        -- get info on data simulation params for testing later
        
        net.density * (net.size - 1) AS net_av_degree,
        net.dcent AS net_centralisation,
        net.size AS net_size
      
      -- join the required tables to supply selected variables
      
      FROM node_results_ranked AS obs_rank
      
      JOIN node_results_GT_ranked AS gt_rank
        ON obs_rank.dataset = gt_rank.dataset
        AND obs_rank.replicate_id = gt_rank.replicate_id
        AND obs_rank.NodeID = gt_rank.NodeID
        
      JOIN network_results_gt AS net
        ON obs_rank.dataset = net.dataset
        AND obs_rank.replicate_id = net.replicate_id
      
    )
    
    -- Pull all the required variables from temp
    
    SELECT
    
      dataset,
      replicate_id,
      -- graph_id, removed, as above
      alpha,
      spotlight_pct,
      b,
      miss_level,
      
      net_size,
      net_av_degree,
      net_centralisation,
      
      corr(obs_degree_rank, gt_degree_rank) AS degree_rank_corr,
      corr(obs_closeness_rank, gt_closeness_rank) AS closeness_rank_corr,
      corr(obs_betweenness_rank, gt_betweenness_rank) AS betweenness_rank_corr,
      corr(obs_eigenvector_rank, gt_eigenvector_rank) AS eigenvector_rank_corr,
      
      corr(obs_degree, gt_degree) AS degree_corr,
      
      -- Closeness is apparently problematic
      
      corr(
      CASE 
        WHEN NOT isnan(obs_closeness)
          AND NOT isnan(gt_closeness)
        THEN obs_closeness
      END,
      CASE
        WHEN NOT isnan(obs_closeness)
          AND NOT isnan(gt_closeness)
        THEN gt_closeness
      END
      ) AS closeness_corr,
      
      corr(obs_betweenness, gt_betweenness) AS betweenness_corr,
      corr(obs_eigenvector, gt_eigenvector) AS eigenvector_corr,
      
      COUNT(*) AS n_nodes
      
    FROM temp
    
    -- Need to aggregate by every other variable
      
    GROUP BY
      
      dataset,
      replicate_id,
      -- graph_id,
      alpha,
      spotlight_pct,
      b,
      miss_level,
      net_size,
      net_av_degree,
      net_centralisation
      
    ORDER BY
      dataset,
      replicate_id,
      alpha,
      spotlight_pct,
      b,
      miss_level;
")

# Diagnostics

node_corr_df %>%
  count(net_size, n_nodes)

node_corr_df %>%
  summarise(
    degree_rank_na = sum(is.na(degree_rank_corr)),
    closeness_rank_na = sum(is.na(closeness_rank_corr)),
    betweenness_rank_na = sum(is.na(betweenness_rank_corr)),
    eigenvector_rank_na = sum(is.na(eigenvector_rank_corr)),
    
    degree_na = sum(is.na(degree_corr)),
    closeness_na = sum(is.na(closeness_corr)),
    betweenness_na = sum(is.na(betweenness_corr)),
    eigenvector_na = sum(is.na(eigenvector_corr))
  )

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


############ Query for data for network level mean absolute bias plots #########

# These should be presented alongside average relative bias plots, to show the 
# difference between the magnitude of bias and the direction of bias

mn_abs_rel_bias_nets <- DBI::dbGetQuery(con, "
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

################# Query for aggregate node-level spotlight effects ##################

# regex at start is to get target centralisation conditions from dataset names
# instead of relying on realised centralisation bands 

node_bias_agg1 <- DBI::dbGetQuery(con, "
WITH gt_aug AS (
  SELECT
    *,
    CASE regexp_extract(dataset, 'c([0-9]+)', 1)
      WHEN '1' THEN 'Low'
      WHEN '3' THEN 'Med'
      WHEN '5' THEN 'High'
      ELSE 'WARNING'
    END AS baseline_centralisation
  FROM node_results_gt
),

node_bias_long AS (

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.graph_id,
    obs.NodeID,
    
    'Degree_raw' AS metric,
    obs.Degree_raw AS observed_value,
    gt.Degree_raw AS gt_value,
    (obs.Degree_raw - gt.Degree_raw) / NULLIF(gt.Degree_raw, 0) AS relative_bias,
    ABS((obs.Degree_raw - gt.Degree_raw) / NULLIF(gt.Degree_raw, 0)) AS abs_relative_bias

  FROM node_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.graph_id,

    obs.NodeID,
    'Betweenness_raw' AS metric,
    obs.Betweenness_raw AS observed_value,
    gt.Betweenness_raw AS gt_value,
    (obs.Betweenness_raw - gt.Betweenness_raw) / NULLIF(gt.Betweenness_raw, 0) AS relative_bias,
    ABS((obs.Betweenness_raw - gt.Betweenness_raw) / NULLIF(gt.Betweenness_raw, 0)) AS abs_relative_bias

  FROM node_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.graph_id,

    obs.NodeID,
    'Closeness_raw' AS metric,
    obs.Closeness_raw AS observed_value,
    gt.Closeness_raw AS gt_value,
    (obs.Closeness_raw - gt.Closeness_raw) / NULLIF(gt.Closeness_raw, 0) AS relative_bias,
    ABS((obs.Closeness_raw - gt.Closeness_raw) / NULLIF(gt.Closeness_raw, 0)) AS abs_relative_bias

  FROM node_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.graph_id,

    obs.NodeID,
    'Eigenvector' AS metric,
    obs.Eigenvector AS observed_value,
    gt.Eigenvector AS gt_value,
    (obs.Eigenvector - gt.Eigenvector) / NULLIF(gt.Eigenvector, 0) AS relative_bias,
    ABS((obs.Eigenvector - gt.Eigenvector) / NULLIF(gt.Eigenvector, 0)) AS abs_relative_bias

  FROM node_results obs
  INNER JOIN gt_aug gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID
),

graph_level_bias AS (
  SELECT
    dataset,
    replicate_id,
    graph_id,
    alpha,
    spotlight_pct,
    b,
    miss_level,
    baseline_centralisation,
    metric,

    AVG(relative_bias) AS graph_mean_relative_bias,
    AVG(abs_relative_bias) AS graph_mean_abs_relative_bias

  FROM node_bias_long
  WHERE relative_bias IS NOT NULL
    AND abs_relative_bias IS NOT NULL

  GROUP BY
    dataset,
    replicate_id,
    graph_id,
    alpha,
    spotlight_pct,
    b,
    miss_level,
    baseline_centralisation,
    metric
),

condition_level_bias AS (
  SELECT
    alpha,
    spotlight_pct,
    b,
    miss_level,
    baseline_centralisation,
    metric,

    AVG(graph_mean_relative_bias) AS mean_relative_bias,
    quantile_cont(graph_mean_relative_bias, 0.25) AS rb_q25,
    quantile_cont(graph_mean_relative_bias, 0.75) AS rb_q75,

    AVG(graph_mean_abs_relative_bias) AS mean_abs_relative_bias,
    quantile_cont(graph_mean_abs_relative_bias, 0.25) AS abs_rb_q25,
    quantile_cont(graph_mean_abs_relative_bias, 0.75) AS abs_rb_q75,

    COUNT(*) AS n_graphs

  FROM graph_level_bias

  GROUP BY
    alpha,
    spotlight_pct,
    b,
    miss_level,
    baseline_centralisation,
    metric
)

SELECT *
FROM condition_level_bias 
ORDER BY
  metric,
  baseline_centralisation,
  alpha,
  spotlight_pct,
  b,
  miss_level;
")

node_bias_agg2 <- node_graph_bias <- DBI::dbGetQuery(con, "
WITH gt_aug AS (
  SELECT
    *,
    CASE regexp_extract(dataset, '_c([0-9]+)$', 1)
      WHEN '1' THEN 'Low'
      WHEN '3' THEN 'Med'
      WHEN '5' THEN 'High'
      ELSE 'WARNING'
    END AS baseline_centralisation
  FROM node_results_gt
),

node_bias_long AS (

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.graph_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.NodeID,

    'Degree_raw' AS metric,

    CASE
      WHEN gt.Degree_raw != 0
      THEN (obs.Degree_raw - gt.Degree_raw) / gt.Degree_raw
      ELSE NULL
    END AS relative_bias,

    CASE
      WHEN gt.Degree_raw != 0
      THEN ABS((obs.Degree_raw - gt.Degree_raw) / gt.Degree_raw)
      ELSE NULL
    END AS abs_relative_bias

  FROM node_results AS obs
  INNER JOIN gt_aug AS gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.graph_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.NodeID,

    'Betweenness_raw' AS metric,

    CASE
      WHEN gt.Betweenness_raw != 0
      THEN (obs.Betweenness_raw - gt.Betweenness_raw) / gt.Betweenness_raw
      ELSE NULL
    END AS relative_bias,

    CASE
      WHEN gt.Betweenness_raw != 0
      THEN ABS((obs.Betweenness_raw - gt.Betweenness_raw) / gt.Betweenness_raw)
      ELSE NULL
    END AS abs_relative_bias

  FROM node_results AS obs
  INNER JOIN gt_aug AS gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.graph_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.NodeID,

    'Closeness_raw' AS metric,

    CASE
      WHEN isfinite(obs.Closeness_raw)
       AND isfinite(gt.Closeness_raw)
       AND gt.Closeness_raw != 0
      THEN (obs.Closeness_raw - gt.Closeness_raw) / gt.Closeness_raw
      ELSE NULL
    END AS relative_bias,

    CASE
      WHEN isfinite(obs.Closeness_raw)
       AND isfinite(gt.Closeness_raw)
       AND gt.Closeness_raw != 0
      THEN ABS((obs.Closeness_raw - gt.Closeness_raw) / gt.Closeness_raw)
      ELSE NULL
    END AS abs_relative_bias

  FROM node_results AS obs
  INNER JOIN gt_aug AS gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID

  UNION ALL

  SELECT
    obs.dataset,
    obs.replicate_id,
    obs.graph_id,
    obs.alpha,
    obs.spotlight_pct,
    obs.b,
    obs.miss_level,
    gt.baseline_centralisation,
    obs.NodeID,

    'Eigenvector' AS metric,

    CASE
      WHEN gt.Eigenvector != 0
      THEN (obs.Eigenvector - gt.Eigenvector) / gt.Eigenvector
      ELSE NULL
    END AS relative_bias,

    CASE
      WHEN gt.Eigenvector != 0
      THEN ABS((obs.Eigenvector - gt.Eigenvector) / gt.Eigenvector)
      ELSE NULL
    END AS abs_relative_bias

  FROM node_results AS obs
  INNER JOIN gt_aug AS gt
    ON obs.dataset = gt.dataset
   AND obs.replicate_id = gt.replicate_id
   AND obs.NodeID = gt.NodeID
),

graph_level_bias AS (
  SELECT
    dataset,
    replicate_id,
    graph_id,
    alpha,
    spotlight_pct,
    b,
    miss_level,
    baseline_centralisation,
    metric,

    AVG(relative_bias) AS graph_mean_relative_bias,
    AVG(abs_relative_bias) AS graph_mean_abs_relative_bias,

    COUNT(*) AS n_nodes_total,
    COUNT(relative_bias) AS n_nodes_used,
    COUNT(*) - COUNT(relative_bias) AS n_nodes_dropped

  FROM node_bias_long

  WHERE baseline_centralisation != 'WARNING'

  GROUP BY
    dataset,
    replicate_id,
    graph_id,
    alpha,
    spotlight_pct,
    b,
    miss_level,
    baseline_centralisation,
    metric
)

SELECT *
FROM graph_level_bias
ORDER BY
  metric,
  baseline_centralisation,
  alpha,
  spotlight_pct,
  b,
  miss_level,
  dataset,
  replicate_id;
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
