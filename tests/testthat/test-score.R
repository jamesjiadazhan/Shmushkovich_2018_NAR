test_that("scoring sums the weights of the bases present", {
  wm <- sdrna_weight_matrix(seed = 1)
  s <- sdrna_score(sdrna356$sequence[1:5], wm)
  expect_length(s, 5L)
  manual <- vapply(sdrna356$sequence[1:5], function(seq) {
    ch <- strsplit(seq, "")[[1]]
    sum(vapply(attr(wm, "positions"),
               function(p) wm[ch[p], p], numeric(1)))
  }, numeric(1), USE.NAMES = FALSE)
  expect_equal(unname(s), manual)
})

test_that("T is accepted as U", {
  wm <- sdrna_weight_matrix(seed = 1)
  seq_u <- sdrna356$sequence[1]
  seq_t <- gsub("U", "T", seq_u)
  expect_equal(sdrna_score(seq_u, wm), sdrna_score(seq_t, wm))
})

test_that("wrong-length sequences are rejected", {
  wm <- sdrna_weight_matrix(seed = 1)
  expect_error(sdrna_score(strrep("A", 30), wm), "nucleotides")
})

test_that("a non-matrix argument is rejected", {
  expect_error(sdrna_score(sdrna356$sequence[1], matrix(0, 4, 50)), "sdrna_wm")
})
