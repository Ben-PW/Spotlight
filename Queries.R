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
