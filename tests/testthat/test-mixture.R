test_that("Single-cluster results match the reference", {
  X =
    rbind(
      matrix( rnorm( 1e3 ), ncol = 10 )
    )
  knockoffs = computeGaussianMixtureKnockoffs(
    X,
    mus = list(10),
    sigmas = list(diag(10)),
    posterior_probs = matrix(1, ncol = 1, nrow = 100)
  )
  knockoffs_reference = computeGaussianKnockoffs(
    X,
    mu = 0,
    sigma = diag(10)
  )
  testthat::expect_equal(knockoffs, knockoffs_reference, tol = 1e-8)
})

test_that("Error control works", {
  X =
    rbind(
      matrix( rnorm( 1e3 ), ncol = 10 ) + 10,
      matrix( rnorm( 1e3 ), ncol = 10 ) - 10
    )
  knockoffs = computeGaussianMixtureKnockoffs(
    X,
    mus = list(10, -10),
    sigmas = list(diag(10), diag(10)),
    posterior_probs =
      data.frame(cluster = rep(c("a", "b"), each = 100)) %>%
      model.matrix(~cluster + 0, .)
  )
  calibration = simulateY(X, knockoffs)
})

