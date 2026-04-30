library(DBI)
library(duckdb)
library(here)

con <- dbConnect(
  duckdb(),
  dbdir = here("Results", "spotlight_results.duckdb")
)

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
      AVG(CASE WHEN Spotlight = 1 THEN true_degree END) AS mean_spotlit,
      AVG(true_degree) AS mean_network
    FROM joined
    GROUP BY dataset, replicate_id, alpha, spotlight_pct
  )

  SELECT
    dataset,
    alpha,
    spotlight_pct,
    AVG(mean_spotlit) AS mean_degree_spotlit,
    AVG(mean_network) AS mean_degree_network,
    AVG(mean_spotlit / mean_network) AS spotlight_degree_lift
  FROM per_network
  GROUP BY dataset, alpha, spotlight_pct
  ORDER BY dataset, alpha, spotlight_pct;
")
