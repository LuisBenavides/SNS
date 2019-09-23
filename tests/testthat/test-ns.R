library("SNS")
context("test-ns")


test_that("rank and normal score", {
  Y = NULL
  X = c(3.0, 4.5, 8.6, 2.3, 2.8, 1.7,6.6)
  R_output = c(4, 5, 7, 2, 3, 1, 6)
  theta <- NULL
  Ftheta <- NULL
  # Test conditional
  ns <- SNS::NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
  expect_equal(ns$Z, c(-0.52440051,-0.38532047,0.08964235), tolerance=1e-7)
  expect_equal(ns$R, c(3.5,4,1), tolerance=1e-10)
})

test_that("normal score and rank unconditional", {
  theta <- NULL
  Ftheta <- NULL
  # Test unconditional
  ns <- SNS::NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
  expect_equal(ns$Z, c(-0.6045853,-0.4727891,-0.2298841), tolerance=1e-7)
  expect_equal(ns$R, c(3.5,4,5), tolerance=1e-10)

})


test_that("normal score and rank unconditional", {
  theta <- NULL
  Ftheta <- NULL
  # Test unconditional
  ns <- SNS::NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
  expect_equal(ns$Z, c(-0.6045853,-0.4727891,-0.2298841), tolerance=1e-7)
  expect_equal(ns$R, c(3.5,4,5), tolerance=1e-10)

})

testthat("Message wrong assignment of Ftheta and theta",{
  theta <- NULL
  Ftheta <- 0.5
  expect_message(SNS::NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta), "ERROR, theta or Ftheta missing")

  theta <- NULL
  Ftheta <- 0.5
  expect_message(SNS::NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta), "ERROR, theta or Ftheta missing")


})
