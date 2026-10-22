test_that("cross-validation returns held-out scores for every compound", {
  cv <- sdrna_cv(k = 5, seed = 1)
  expect_length(cv$scores, nrow(sdrna356))
  expect_false(anyNA(cv$scores))
  expect_equal(sort(unique(cv$fold)), 1:5)
})

test_that("out-of-sample AUC is below the resubstitution AUC", {
  wm <- sdrna_weight_matrix(functional_cutoff = 24, seed = 1)
  insample <- sdrna_auc(sdrna_score(sdrna356$sequence, wm),
                        sdrna356$expression < 24)
  cv <- sdrna_cv(functional_cutoff = 24, k = 10, seed = 1)
  expect_lt(cv$auc, insample)
  expect_gt(cv$auc, 0.5)
})

test_that("a fold with too few compounds in a subset is rejected cleanly", {
  tiny <- sdrna356[c(which(sdrna356$expression < 24)[1:3],
                     which(sdrna356$expression > 44)[1:3]), ]
  expect_error(sdrna_cv(data = tiny, k = 3, seed = 1), "too few")
})

test_that("k must be at least 2", {
  expect_error(sdrna_cv(k = 1), "k >= 2")
})
