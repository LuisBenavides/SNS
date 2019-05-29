
<!-- README.md is generated from README.Rmd. Please edit that file -->
Sequential Normal Scores
========================

The methods discussed in this package are new nonparametric methods based on *sequential normal scores* (SNS), designed for sequences of observations, usually time series data, which may occur singly or in batches, and may be univariate or multivariate. These methods are designed to detect changes in the process, which may occur as changes in *location* (mean or median), changes in *scale* (standard deviation, or variance), or other changes of interest in the distribution of the observations, over the time observed. They usually apply to large data sets, so computations need to be simple enough to be done in a reasonable time on a computer, and easily updated as each new observation (or batch of observations) becomes available.

Installation
------------

You can install the released version of SNS from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SNS")
```

or install from the package hosted in [github](https://github.com/LuisBenavides/sns).

``` r
install_github("LuisBenavides/SNS")
```

**Note**: To use `install_github` it is needed the library devtools. \#\# Univariate Analysis

### Using SNS function

The reference sample `Y` and the monitoring sample `X`.

``` r
Y = c(10,20,30,40,50,60,70,80,90,100)
X = c(30, 35, 45)
```

#### Example of conditionsl SNS with a reference sample `Y`

``` r
Y = c(10,20,30,40,50,60,70,80,90,100)
X = c(30, 35, 45)
theta = 40
Ftheta = 0.5
sample.id = c("a", "b", "c")
SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
```

Output

    #> $coefficients
    #> $coefficients$n
    #> [1] 1
    #> 
    #> $coefficients$chart
    #> [1] "Shewhart"
    #> 
    #> $coefficients$chart.par
    #> [1] 3
    #> 
    #> 
    #> $Z
    #> [1] -0.52440051 -0.31863936  0.08964235
    #> 
    #> $X.id
    #> [1] "a" "b" "c"
    #> 
    #> $UCL
    #> [1] 3 3 3
    #> 
    #> $LCL
    #> [1] -3 -3 -3
    #> 
    #> $scoring
    #> [1] "Z"
    #> 
    #> attr(,"class")
    #> [1] "SNS"

#### Example of unconditionsl SNS with a reference sample `Y`

``` r
Y = c(10,20,30,40,50,60,70,80,90,100)
X = c(30, 35, 45)
theta = NULL
Ftheta = NULL
sample.id = c("a", "b", "c")
SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
```

Output

    #> $coefficients
    #> $coefficients$n
    #> [1] 1
    #> 
    #> $coefficients$chart
    #> [1] "Shewhart"
    #> 
    #> $coefficients$chart.par
    #> [1] 3
    #> 
    #> 
    #> $Z
    #> [1] -0.6045853 -0.3186394  0.0000000
    #> 
    #> $X.id
    #> [1] "a" "b" "c"
    #> 
    #> $UCL
    #> [1] 3 3 3
    #> 
    #> $LCL
    #> [1] -3 -3 -3
    #> 
    #> $scoring
    #> [1] "Z"
    #> 
    #> attr(,"class")
    #> [1] "SNS"

#### Example of conditional SNS without a reference sample `Y`

``` r
Y = NULL
X = c(30, 35, 45)
theta = 40
Ftheta = 0.5
sample.id = c("a", "b", "c")
SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
```

Output

    #> $coefficients
    #> $coefficients$n
    #> [1] 1
    #> 
    #> $coefficients$chart
    #> [1] "Shewhart"
    #> 
    #> $coefficients$chart.par
    #> [1] 3
    #> 
    #> 
    #> $Z
    #> [1] -0.6744898 -0.2104284  0.6744898
    #> 
    #> $X.id
    #> [1] "a" "b" "c"
    #> 
    #> $UCL
    #> [1] 3 3 3
    #> 
    #> $LCL
    #> [1] -3 -3 -3
    #> 
    #> $scoring
    #> [1] "Z"
    #> 
    #> attr(,"class")
    #> [1] "SNS"

#### Example of unconditional SNS without a reference sample `Y`

``` r
Y = NULL
X = c(30, 35, 45)
theta = NULL
Ftheta = NULL
sample.id = c("a", "b", "c")
SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
```

Output

    #> $coefficients
    #> $coefficients$n
    #> [1] 1
    #> 
    #> $coefficients$chart
    #> [1] "Shewhart"
    #> 
    #> $coefficients$chart.par
    #> [1] 3
    #> 
    #> 
    #> $Z
    #> [1] 0.0000000 0.9674216 1.1503494
    #> 
    #> $X.id
    #> [1] "a" "b" "c"
    #> 
    #> $UCL
    #> [1] 3 3 3
    #> 
    #> $LCL
    #> [1] -3 -3 -3
    #> 
    #> $scoring
    #> [1] "Z"
    #> 
    #> attr(,"class")
    #> [1] "SNS"
