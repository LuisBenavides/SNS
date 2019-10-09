library("SNS")

context("test-sns")

test_that("Sequential rank in batches (without ties)", {
  R_test = c(3, 1, 2, 3, 3, 4, 7, 2, 2, 2, 9, 10)
  Z_test = c(0.97, -0.97, 0, 0.32, 0.32, 1.15, 1.47, -0.79, -0.79, -1.04, 1.04, 1.64)

  n = 3
  X = c(6.1, 0.6, 3.6, 4.1, 5.6, 6.7, 8.6, 3.5, 2.5, 1.6, 6.9, 9.8)

  X.id = rep(1:(length(X)/n),each=n)
  output = SNS(X=X,X.id=X.id,snsRaw = TRUE)

  R = output$Rraw
  Z = round(output$Zraw,2)

  expect_equal(R,R_test)
  expect_equal(Z,Z_test)
})

test_that("Sequential rank and sequential normal score in batches of single observations (with ties)", {
  R_test = c(1, 1, 2.5, 2, 4, 1)
  Z_test = c(0, -0.67, 0.43, -0.32, 0.52, -1.38)

  n = 1
  X = c(4.99, 2.14, 4.99, 2.52, 4.99, 1.31)

  X.id = rep(1:(length(X)/n),each=n)
  output = SNS(X=X,X.id=X.id,snsRaw = TRUE)


  R = output$Rraw
  Z = round(output$Zraw,2)


  expect_equal(R, R_test)
  expect_equal(Z, Z_test)
})
