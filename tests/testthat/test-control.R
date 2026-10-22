test_that("a reference composition replaces the non-functional side", {
  L <- 50L
  ref <- matrix(0.25, 4L, L, dimnames = list(c("A","C","G","U"), seq_len(L)))
  wm <- sdrna_weight_matrix(reference = ref, seed = 1)
  m <- sdrnapred:::sdrna_split(sdrna356$sequence)
  cf <- sdrnapred:::sdrna_counts(m[sdrna356$expression < 24, ])
  expected <- (cf / 91 - 0.25) * 100
  sig <- attr(wm, "significant")
  expect_equal(as.vector(wm[sig]), as.vector(expected[sig]))
})

test_that("the control matrix differs from the main matrix", {
  main <- sdrna_weight_matrix(seed = 1)
  ctl <- sdrna_control_matrix(seed = 1)
  expect_false(isTRUE(all.equal(unclass(main), unclass(ctl))))
})
