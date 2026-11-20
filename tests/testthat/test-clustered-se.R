# The standard error clustered by collection period decides which differences
# the article marks as credible.

test_that("with one observation per cluster it is the usual standard error", {
  set.seed(1)
  d <- rnorm(50)
  expect_equal(clustered_se(d, seq_along(d)), sqrt(length(d) * var(d)))
})

test_that("it matches an example worked by hand", {
  d <- c(1, 2, 3, 4, 5, 6)
  cluster <- c("a", "a", "b", "b", "c", "c")
  # Mean 3.5. Centred values: -2.5 -1.5 | -0.5 0.5 | 1.5 2.5.
  # Cluster sums: -4, 0, 4. Sum of squares 32. Three clusters: sqrt(3/2 * 32).
  expect_equal(clustered_se(d, cluster), sqrt(48))
})

test_that("it grows when the observations of a cluster move together", {
  d <- c(1, 1, -1, -1, 2, 2, -2, -2)
  together <- c(1, 1, 2, 2, 3, 3, 4, 4)
  apart <- c(1, 2, 1, 2, 3, 4, 3, 4)
  expect_gt(clustered_se(d, together), sqrt(length(d) * var(d)))
  expect_lt(clustered_se(d, apart), sqrt(length(d) * var(d)))
})

test_that("it does not depend on the order of the observations or on the labels", {
  set.seed(2)
  d <- rnorm(30)
  cluster <- rep(1:10, each = 3)
  shuffle <- sample(30)
  expect_equal(clustered_se(d[shuffle], cluster[shuffle]), clustered_se(d, cluster))
  expect_equal(clustered_se(d, letters[cluster]), clustered_se(d, cluster))
})

test_that("it needs two clusters and one label per value", {
  expect_true(is.na(clustered_se(c(1, 2, 3), c("a", "a", "a"))))
  expect_error(clustered_se(1:3, 1:2))
})

# Stand-ins for leave-one-out results: all that paired_elpd() reads is the
# pointwise ELPD column and the collection period of each observation.
fake_loo <- function(elpd, period) {
  l <- list(pointwise = cbind(elpd_loo = elpd))
  attr(l, "period") <- period
  l
}

test_that("paired_elpd() gives the summed difference and both standard errors", {
  period <- c("p1", "p1", "p2", "p3")
  a <- fake_loo(c(-1, -2, -3, -4), period)
  b <- fake_loo(c(-1.5, -2.5, -2, -5), period)
  d <- c(0.5, 0.5, -1, 1)
  p <- paired_elpd(a, b)
  expect_equal(p$delta, 1)
  expect_equal(p$se, sqrt(4 * var(d)))
  expect_equal(p$se_cl, clustered_se(d, period))
  expect_equal(p$n, 4)
  expect_equal(paired_elpd(b, a)$delta, -1)
})

test_that("paired_elpd() refuses models fitted to different observations", {
  a <- fake_loo(c(-1, -2, -3), c("p1", "p2", "p3"))
  expect_error(paired_elpd(a, fake_loo(c(-1, -2, -3), c("p1", "p2", "p4"))), "same observations")
  expect_error(paired_elpd(a, fake_loo(c(-1, -2), c("p1", "p2"))), "same observations")
  expect_error(paired_elpd(a, list(pointwise = cbind(elpd_loo = c(-1, -2, -3)))), "same observations")
})

test_that("loo_period_ids() labels exactly the rows a model would use", {
  data <- data.frame(
    y = c(3, NA, 5, 2, 1),
    x = c(0.1, 0.2, NA, 0.4, 0.5),
    unused = c(NA, 1, 1, 1, NA),
    trap_name = c("A", "A", "B", "B", "C"),
    start_date = as.Date("2019-06-01") + c(0, 7, 0, 7, 0),
    end_date = as.Date("2019-06-08") + c(0, 7, 0, 7, 0)
  )
  model <- list(formula = list(formula = y ~ x + (1 | trap_name)))
  # Rows 2 (no outcome) and 3 (no covariate) are dropped; a gap in a column the
  # model does not use changes nothing. Traps served on the same dates share a
  # label.
  expect_identical(loo_period_ids(model, data),
                   c("2019-06-01 2019-06-08", "2019-06-08 2019-06-15", "2019-06-01 2019-06-08"))
})
