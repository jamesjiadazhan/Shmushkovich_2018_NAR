test_that("bundled training set matches the published table", {
  expect_equal(nrow(sdrna356), 356L)
  expect_equal(sdrna356$id, 1000:1355)
  expect_true(all(nchar(sdrna356$sequence) == 50L))
  expect_false(any(grepl("[^ACGU]", sdrna356$sequence)))
})

test_that("subset sizes reproduce the paper's stated bins", {
  expect_equal(sum(sdrna356$expression < 17), 41L)
  expect_equal(sum(sdrna356$expression < 24), 91L)
  expect_equal(sum(sdrna356$expression < 35), 138L)
  expect_equal(sum(sdrna356$expression > 44), 157L)
})

test_that("documented irregular values are present as printed", {
  expect_equal(sdrna356$expression[sdrna356$id == 1095], 25.405)
  expect_equal(sdrna356$expression[sdrna356$id == 1107], 29.2)
})
