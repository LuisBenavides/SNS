library("SNS")

context("test-sns")

test_that("sequential normal scores", {
  Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
  X <- c(30, 35, 45)
  theta <- 40
  Ftheta <- 0.5
  sample.id <- c("a", "b", "c")
  # Test CONDITIONAL WITH REFERENCE SAMPLE
  output <- SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
  expectedOutput <- c(-0.52440051,-0.31863936,0.08964235)
  expect_equal(output$Z, expectedOutput, tolerance=1e-7)


  # Test UNCONDITIONAL WITH REFERENCE SAMPLE
  theta <- NULL
  Ftheta <- NULL
  output <- SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
  expectedOutput <- c(-0.6045853,-0.3186394,0.0000000)
  expect_equal(output$Z, expectedOutput, tolerance=1e-7)

  # EXAMPLE CONDITIONAL WITHOUT REFERENCE SAMPLE
  Y <- NULL
  theta <- 40
  Ftheta <- 0.5
  output <- SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
  expectedOutput <- c(-0.6744898,-0.2104284,0.6744898)
  expect_equal(output$Z, expectedOutput, tolerance=1e-7)

  # EXAMPLE UNCONDITIONAL WITHOUT REFERENCE SAMPLE
  theta <- NULL
  Ftheta <- NULL
  output <- SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
  expectedOutput <- c(0.0000000,0.9674216,1.1503494)
  expect_equal(output$Z, expectedOutput, tolerance=1e-7)
})
