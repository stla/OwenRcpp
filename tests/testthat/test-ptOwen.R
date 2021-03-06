context("ptOwen")

test_that("ptOwen = pt", {
  ncp <- 1
  expect_equal(ptOwen(2, nu=3, delta=ncp), pt(2, df=3, ncp=ncp), tolerance=1e-13)
  expect_equal(ptOwen(2, nu=4, delta=ncp), pt(2, df=4, ncp=ncp), tolerance=1e-12)
  ncp <- 10
  expect_equal(ptOwen(2, nu=3, delta=ncp), pt(2, df=3, ncp=ncp), tolerance=1e-12)
  expect_equal(ptOwen(2, nu=4, delta=ncp), pt(2, df=4, ncp=ncp), tolerance=1e-13)
})