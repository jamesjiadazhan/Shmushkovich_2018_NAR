test_that("AUC matches a direct pairwise computation", {
  set.seed(7)
  s <- rnorm(60); l <- s + rnorm(60) > 0
  pos <- s[l]; neg <- s[!l]
  brute <- mean(outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b)))
  expect_equal(sdrna_auc(s, l), brute)
})

test_that("a perfect and an inverted ranking give 1 and 0", {
  s <- c(5, 4, 3, 2, 1); l <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  expect_equal(sdrna_auc(s, l), 1)
  expect_equal(sdrna_auc(-s, l), 0)
})

test_that("AUC is NA when a class is empty", {
  expect_true(is.na(sdrna_auc(c(1, 2, 3), c(FALSE, FALSE, FALSE))))
})

test_that("PPP and sensitivity are consistent at the extremes", {
  wm <- sdrna_weight_matrix(seed = 1)
  s <- sdrna_score(sdrna356$sequence, wm)
  curve <- sdrna_ppp_sensitivity(s, sdrna356$expression, 24)
  expect_equal(nrow(curve), 356L)
  # Taking every compound recovers the base rate and full sensitivity.
  expect_equal(curve$sensitivity[356], 100)
  expect_equal(curve$ppp[356], 100 * 91 / 356)
  expect_true(all(diff(curve$n_true) >= 0))
})

test_that("sdrna_ppp_at returns the first qualifying threshold", {
  wm <- sdrna_weight_matrix(seed = 1)
  s <- sdrna_score(sdrna356$sequence, wm)
  curve <- sdrna_ppp_sensitivity(s, sdrna356$expression, 24)
  pt <- sdrna_ppp_at(curve, 40)
  expect_equal(nrow(pt), 1L)
  expect_gte(pt$sensitivity, 40)
})
