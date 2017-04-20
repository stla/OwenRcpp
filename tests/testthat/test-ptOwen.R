context("ptOwen")

test_that("ptOwen = pt", {
  ncp <- 1
  expect_equal(ptOwen(2, nu=3, delta=ncp), pt(2, df=3, ncp=ncp), tolerance=1e-6)
  expect_equal(ptOwen(2, nu=4, delta=ncp), pt(2, df=4, ncp=ncp))
})