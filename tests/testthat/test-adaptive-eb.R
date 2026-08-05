test_that("adaptive EB weights are bounded and continuous", {
  fit <- combatrefql:::fit_adaptive_eb(c(-0.05, 0, 0.05, 0.1, 2),
    rep(0.04, 5))
  expect_true(all(fit$adaptive_weight >= 0 & fit$adaptive_weight <= 1))
  expect_true(all(fit$adaptive_weight >= fit$weight))
  expect_lt(abs(fit$posterior[2] - fit$standard_posterior[2]), 0.01)
  expect_lt(abs(fit$posterior[5] - 2),
    abs(fit$standard_posterior[5] - 2))
})

test_that("tau2 zero and invalid contrasts remain finite and explicit", {
  pooled <- combatrefql:::fit_adaptive_eb(rep(1, 4), rep(0.1, 4))
  expect_identical(pooled$variance, 0)
  expect_true(all(is.finite(pooled$posterior)))
  invalid <- combatrefql:::fit_adaptive_eb(c(1, NA_real_), c(0.1, NA_real_))
  expect_identical(invalid$status, "insufficient_genes")
  expect_true(is.na(invalid$posterior[2]))
})

test_that("production uses adaptive EB with stable count contracts", {
  z <- combatrefql_test_fixtures()[[1L]]
  a <- combat_ref_ql(z$counts, z$batch, z$group, reference = "1",
    verbosity = "quiet")
  b <- combat_ref_ql(z$counts, z$batch, z$group, reference = "1",
    verbosity = "quiet")
  d <- a@diagnostics$batch_contrasts

  expect_identical(a@counts, b@counts)
  expect_true(all(d$adaptive_weight >= 0 & d$adaptive_weight <= 1,
    na.rm = TRUE))
  expect_true(all(is.finite(a@counts) & a@counts >= 0 &
    a@counts == round(a@counts)))
  expect_identical(dim(a@counts), dim(z$counts))
  expect_identical(dimnames(a@counts), dimnames(z$counts))
  expect_identical(a@counts[, z$batch == "1", drop = FALSE],
    z$counts[, z$batch == "1", drop = FALSE])
  expect_identical(a@specification$transport, "negative-binomial mid-P")
  expect_true(any(a@diagnostics$dispersion$batch != "1"))

  expected <- matrix(c(
    10, 0, 1, 6, 3, 60, 12, 0, 2, 7, 4, 65,
    25, 0, 0, 15, 6, 151, 27, 0, 0, 17, 7, 161,
    8, 0, 3, 7, 2, 55, 9, 0, 2, 8, 3, 58,
    22, 0, 0, 19, 8, 141, 25, 0, 0, 20, 9, 156
  ), nrow = 6L)
  expect_identical(unname(a@counts), expected)
})
