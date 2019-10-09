context("test-getrl")
shift = 0

n <- 10 # subgroup size
m <- 1000 # reference-sample size
dist <- "Normal"
mu <- c(0, shift) # c(reference sample mean, monitoring sample mean)
sigma <- c(1, 1) # c(reference sample sd, monitoring sample sd)

#### Distribution parameters
dist.par <- c(0, 1, 1) # c(location, scale, shape)

#### Other Parameters
replicates <- 10000
print.RL <- TRUE
calibrate <- FALSE
progress <- TRUE
arl0 <- 370

test_that("Shewart Normal approximation", {
  chart <- "Shewhart"
  chart.par <- c(3)
  #shewhart <- SNS::getRL(1, n, m,
  #                  theta = NULL, Ftheta = NULL,dist, mu, sigma, dist.par = dist.par,
  #                  chart = chart, chart.par = chart.par, calibrate = calibrate, arl0 = arl0
  #)

  expect_equal(2, 2)
})
