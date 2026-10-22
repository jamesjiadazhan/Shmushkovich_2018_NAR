test_that("weight matrix has the expected shape and zeroing", {
  wm <- sdrna_weight_matrix(seed = 1)
  expect_s3_class(wm, "sdrna_wm")
  expect_equal(dim(wm), c(4L, 50L))
  expect_equal(rownames(wm), c("A", "C", "G", "U"))
  # Positions 1 and 50 are excluded from the primary algorithm.
  expect_true(all(wm[, 1] == 0))
  expect_true(all(wm[, 50] == 0))
  # Non-significant cells are zero.
  sig <- attr(wm, "significant")
  expect_true(all(wm[!sig] == 0))
})

test_that("subset sizes are carried through", {
  wm <- sdrna_weight_matrix(functional_cutoff = 24, seed = 1)
  expect_equal(attr(wm, "n_functional"), 91L)
  expect_equal(attr(wm, "n_nonfunctional"), 157L)
})

test_that("weights are the functional minus non-functional frequency gap", {
  wm <- sdrna_weight_matrix(seed = 1)
  m <- sdrnapred:::sdrna_split(sdrna356$sequence)
  f <- sdrna356$expression < 24
  nf <- sdrna356$expression > 44
  cf <- sdrnapred:::sdrna_counts(m[f, ])
  cnf <- sdrnapred:::sdrna_counts(m[nf, ])
  expected <- (cf / sum(f) - cnf / sum(nf)) * 100
  sig <- attr(wm, "significant")
  expect_equal(as.vector(wm[sig]), as.vector(expected[sig]))
})

test_that("restricting to the targeting region works", {
  wm <- sdrna_weight_matrix(positions = 16:35, seed = 1)
  expect_true(all(wm[, 1:15] == 0))
  expect_true(all(wm[, 36:50] == 0))
  expect_equal(attr(wm, "positions"), 16:35)
})

test_that("the simulation is reproducible under a seed", {
  a <- sdrna_weight_matrix(seed = 42)
  b <- sdrna_weight_matrix(seed = 42)
  expect_equal(unclass(a), unclass(b))
})

test_that("invalid input is rejected", {
  expect_error(sdrna_weight_matrix(functional_cutoff = 50), "below")
  expect_error(sdrna_weight_matrix(positions = 1:99), "within")
})
