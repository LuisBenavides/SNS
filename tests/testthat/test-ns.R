library("SNS")
context("test-sns")

test_that("normal score and rank", {
  Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  X <- c(30, 35, 45)
  theta <- 40
  Ftheta <- 0.5
  # Test conditional
  ns <- NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
  #expect_equal(ns$Z, c(-0.52440051,-0.31863936,0.08964235), tolerance=1e-7)
  #expect_equal(ns$R, c(3.5,4,1), tolerance=1e-10)

  theta <- NULL
  Ftheta <- NULL
  # Test unconditional
  ns <- NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
  #expect_equal(ns$Z, c(-0.6045853,-0.4727891,-0.2298841), tolerance=1e-7)
  #expect_equal(ns$R, c(3.5,4,5), tolerance=1e-10)

})
