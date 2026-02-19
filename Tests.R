
library(testthat)

test_that("assignSpotlight assigns expected vertices and edges given seed", {
  library(igraph)
  
  g  <- make_star(8, mode = "undirected")  # clear degree differences
  gl <- list(g)
  
  spotlight_pct <- 0.25
  alpha <- 2
  
  # expected
  withr::with_seed(123, {
    n <- vcount(g)
    k <- max(1L, round(n * spotlight_pct))
    deg <- degree(g)
    w <- (deg + 1)^alpha
    
    idx <- sample.int(n, size = k, replace = FALSE, prob = w)
    Spotlight_expected <- integer(n)
    Spotlight_expected[idx] <- 1L
  })
  
  # actual 
  out <- withr::with_seed(123, {
    assignSpotlight(gl, spotlight_pct = spotlight_pct, alpha = alpha)
  })
  g2 <- out[[1]]
  
  # node Spotlight matches expected
  expect_equal(V(g2)$Spotlight, Spotlight_expected)
  
  # edge Spotlight matches OR of endpoints
  ends_mat <- ends(g2, E(g2), names = FALSE)
  edge_expected <- as.integer(
    Spotlight_expected[ends_mat[, 1]] | Spotlight_expected[ends_mat[, 2]]
  )
  expect_equal(E(g2)$Spotlight, edge_expected)
})

test_that("assignSpotlight always assigns at least 1 spotlit node", {
  library(igraph)
  
  g2 <- assignSpotlight(list(make_ring(10)), spotlight_pct = 0, alpha = 0)[[1]]
  expect_equal(sum(V(g2)$Spotlight), 1L)
})

test_that("assignSpotlight works when graph has 0 edges", {
  library(igraph)
  
  g <- make_empty_graph(n = 6, directed = FALSE)
  g2 <- assignSpotlight(list(g), spotlight_pct = 0.5, alpha = 0)[[1]]
  
  expect_true(is.null(E(g2)$Spotlight) == FALSE)  # attribute exists
  expect_equal(length(E(g2)$Spotlight), 0L)
})

test_that("sampleSpotlight deletes expected edges given seed", {
  library(igraph)
  
  g <- make_ring(12)
  
  # assign Spotlight deterministically, then sample spotlight deterministically
  g <- withr::with_seed(1, assignSpotlight(list(g), spotlight_pct = 0.25, alpha = 0)[[1]])
  
  miss_level <- 0.3
  b <- 3
  
  m <- ecount(g)
  k <- round(m * miss_level)
  sp <- E(g)$Spotlight
  w  <- ifelse(sp == 0, b, 1)
  
  # expected dropped indices
  idx_drop <- withr::with_seed(999, sample.int(m, size = k, replace = FALSE, prob = w))
  expected <- delete_edges(g, E(g)[idx_drop])
  
  # actual
  actual <- withr::with_seed(999, sampleSpotlight(list(g), miss_level = miss_level, b = b)[[1]])
  
  expect_true(isomorphic(expected, actual))
})

test_that("sampleSpotlight miss_level=0 returns unchanged", {
  library(igraph)
  
  g <- make_ring(10)
  g <- assignSpotlight(list(g), spotlight_pct = 0.2, alpha = 0)[[1]]
  out <- sampleSpotlight(list(g), miss_level = 0, b = 2)[[1]]
  
  expect_true(isomorphic(g, out))
  expect_equal(ecount(out), ecount(g))
})

test_that("sampleSpotlight miss_level=1 removes all edges (when m>0)", {
  library(igraph)
  
  g <- make_ring(10)
  g <- assignSpotlight(list(g), spotlight_pct = 0.2, alpha = 0)[[1]]
  out <- sampleSpotlight(list(g), miss_level = 1, b = 2)[[1]]
  
  expect_equal(ecount(out), 0L)
})

test_that("sampleSpotlight errors if E(g)$Spotlight is missing", {
  library(igraph)
  
  g <- make_ring(10)
  expect_error(sampleSpotlight(list(g), miss_level = 0.2, b = 2))
})

