test_that("startup banner is versioned and suppressible", {
  attach <- getFromNamespace(".onAttach", "combatrefql")
  expect_message(attach(.libPaths()[1L], "combatrefql"),
    "ComBat-refQL 0.0.1", fixed = TRUE)
  expect_silent(suppressPackageStartupMessages(
    attach(.libPaths()[1L], "combatrefql")))
})

test_that("fitting output reports only the major stages", {
  z <- combatrefql_test_fixtures()[[1L]]
  fractional <- z$counts
  fractional[1L, 1L] <- fractional[1L, 1L] + 0.25
  output <- paste(capture.output(
    fit <- combat_ref_ql(fractional, z$batch, z$group, reference = "1"),
    type = "message"), collapse = "\n")

  expect_match(output, "Rounded 1 fractional values across 1 genes")
  expect_match(output, "Prepared 6 genes × 8 samples")
  expect_match(output, "Selected reference: 1")
  expect_match(output, "Fitted preliminary model")
  expect_match(output, "Estimated hierarchical batch-specific NB dispersion")
  expect_match(output, "Fitted final model")
  expect_match(output, "Adjusted 4 | Unchanged 2 | Failed 0", fixed = TRUE)
  expect_match(output, "Completed in [0-9.]+ seconds")
  expect_false(any(grepl("Normalized with|edgeR v4 QL|Moderated batch|Constructed and validated", output)))
  expect_silent(combat_ref_ql(z$counts, z$batch, z$group, reference = "1",
    verbosity = "quiet"))
})

test_that("fit printing is compact and summary retains diagnostics", {
  z <- combatrefql_test_fixtures()[[1L]]
  fit <- combat_ref_ql(z$counts, z$batch, z$group, reference = "1",
    verbosity = "quiet")
  printed <- paste(capture.output(print(fit), type = "message"), collapse = "\n")
  expect_match(printed, "6 genes × 8 samples | 2 batches | 2 biological groups")
  expect_match(printed, "Batch effects +Evidence-adaptive EB")
  expect_match(printed, "Dispersion +Hierarchical NB")
  expect_match(printed, "Normalization +TMMwsp")
  expect_match(printed, "Transport +NB mid-P")
  expect_match(printed, "Reference +1")
  expect_match(printed, "Adjusted +4 \\(66.7%\\)")
  expect_match(printed, "Unchanged +2 \\(33.3%\\)")
  expect_match(printed, "Failed +0")
  expect_match(printed, "Runtime +[0-9.]+ s")

  diagnostic <- summary(fit)
  expect_true(all(c("samples", "coefficients", "rank", "residual_df",
    "condition_number") %in% names(diagnostic@design)))
  expect_true(all(c("reference_score", "second_best", "score_margin") %in%
    names(diagnostic@batches)))
  expect_true(all(c("multiplier", "status") %in% names(diagnostic@dispersion)))
  expect_true(all(c("median_weight", "status") %in% names(diagnostic@shrinkage)))
  expect_identical(diagnostic@gene_outcomes$outcome,
    c("adjusted", "unsupported", "failed"))
  expect_true(all(c("preparation", "normalization", "reference",
    "preliminary_ql_fit", "dispersion", "final_ql_fit",
    "contrast_moderation", "mapping", "total") %in% diagnostic@timings$stage))
})

test_that("adaptive EB preserves the fitted QL and dispersion model", {
  z <- combatrefql_test_fixtures()[[1L]]
  fit <- combat_ref_ql(z$counts, z$batch, z$group, reference = "1",
    verbosity = "quiet")
  expected_coefficients <- matrix(c(
    -2.07359417001567, -2.56167047599389, -3.37721570778924,
    -0.304828792768004, -1.74765553439272, -2.79208758333197,
    -3.54462386632599, -0.310573988562776, -0.0934468024984378,
    0.223664611963277, 0.163329486988363, -0.00760937495690332
  ), nrow = 4L)

  expect_equal(unname(fit@model$coefficients), expected_coefficients,
    tolerance = 1e-12)
  expect_equal(unname(fit@model$nb_dispersion),
    rep(9.76562563324275e-06, 4L), tolerance = 1e-12)
  expect_equal(unname(fit@model$ql_dispersion),
    c(0.0311771154128187, 0.055907752553489, 0.10374820379279,
      0.00708625931717186), tolerance = 1e-12)
  expect_equal(unname(fit@model$posterior_ql_dispersion),
    c(0.031172337929441, 0.0558843287582505, 0.103795489369206,
      0.00708642319654715), tolerance = 1e-12)
  expect_equal(unname(fit@model$batch_dispersion_shifts),
    c(0, -3.87215148833775e-10), tolerance = 1e-12)
  expect_equal(unname(fit@model$batch_dispersion_multipliers),
    c(1, 0.999999999612785), tolerance = 1e-12)
  expect_equal(unname(fit@model$reference_dispersion),
    rep(9.76562563324275e-06, 4L), tolerance = 1e-12)
})
