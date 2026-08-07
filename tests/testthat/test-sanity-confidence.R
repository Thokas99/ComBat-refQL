sanity_fixture <- function(batch = rep(c("A", "B"), each = 6L), seed = 901L) {
  set.seed(seed)
  counts <- matrix(rnbinom(120L * length(batch), mu = 40, size = 8), 120L,
    dimnames = list(paste0("g", 1:120), paste0("s", seq_along(batch))))
  list(counts = counts, batch = factor(batch))
}

test_that("zero-covariate workflow is first class and deterministic", {
  z <- sanity_fixture()
  first <- combat_ref_ql(z$counts, z$batch, verbose = FALSE)
  second <- combat_ref_ql(z$counts, z$batch, verbose = FALSE)
  expect_identical(first@counts, second@counts)
  expect_null(attr(first@counts, "confidence"))
  expect_true(all(first@diagnostics$mapping$reference_invariant))
  expect_true(all(first@counts[, z$batch == first@specification$reference_batch] ==
    z$counts[, z$batch == first@specification$reference_batch]))
})

test_that("fundamental count and metadata failures are typed", {
  z <- sanity_fixture()
  cases <- list(
    no_gene_ids = `rownames<-`(z$counts, NULL),
    duplicate_gene_ids = `rownames<-`(z$counts, rep("g", nrow(z$counts))),
    no_sample_ids = `colnames<-`(z$counts, NULL),
    empty = matrix(numeric(), 0, 0),
    non_numeric = matrix("1", 2, 2, dimnames = list(c("g1", "g2"), c("s1", "s2"))))
  for (counts in cases) expect_error(
    combat_ref_ql(counts, z$batch, verbose = FALSE),
    class = "combatrefql_input_error")
  for (value in c(NA_real_, Inf)) {
    bad <- z$counts; bad[1, 1] <- value
    expect_error(combat_ref_ql(bad, z$batch, verbose = FALSE),
      class = "combatrefql_input_error")
  }
  zero_sample <- z$counts; zero_sample[, 1] <- 0
  expect_error(combat_ref_ql(zero_sample, z$batch, verbose = FALSE),
    class = "combatrefql_input_error")
  expect_error(combat_ref_ql(z$counts, rep("A", ncol(z$counts)), verbose = FALSE),
    class = "combatrefql_input_error")
})

test_that("singleton and two-sample source batches use existing fallback conservatively", {
  singleton <- sanity_fixture(c(rep("A", 6), "B"))
  fit1 <- combat_ref_ql(singleton$counts, singleton$batch, reference = "A",
    verbose = FALSE)
  b1 <- fit1@diagnostics$batch_quality[fit1@diagnostics$batch_quality$batch == "B", ]
  c1 <- fit1@diagnostics$correction_confidence
  expect_false(b1$reference_eligible)
  expect_equal(b1$status, "low_confidence")
  expect_true(all(c1$confidence_label == "low"))
  expect_lte(max(c1$confidence_score), 0.25)

  selected <- combat_ref_ql(singleton$counts, singleton$batch, verbose = FALSE)
  expect_equal(selected@specification$reference_batch, "A")

  two <- sanity_fixture(c(rep("A", 6), "B", "B"))
  fit2 <- combat_ref_ql(two$counts, two$batch, reference = "A", verbose = FALSE)
  b2 <- fit2@diagnostics$batch_quality[fit2@diagnostics$batch_quality$batch == "B", ]
  expect_equal(b2$status, "limited_information")
  expect_lte(max(fit2@diagnostics$correction_confidence$confidence_score), 0.5)
  expect_false(any(fit2@diagnostics$correction_confidence$confidence_label == "high"))
})

test_that("design association warns heuristically while exact aliasing errors", {
  z <- sanity_fixture(rep(c("A", "B"), each = 24L))
  group <- factor(c(rep("ctrl", 23), "treated", "ctrl", rep("treated", 23)))
  expect_warning(fit <- combat_ref_ql(z$counts, z$batch, group,
    reference = "A", verbose = FALSE), "strongly associated")
  expect_gte(fit@diagnostics$input$maximum_entanglement, 0.8)
  expect_error(combat_ref_ql(z$counts, z$batch, z$batch, verbose = FALSE),
    class = "combatrefql_design_error")
})

test_that("confidence invariants and conservative gene aggregation hold", {
  z <- sanity_fixture(rep(c("A", "B", "C"), each = 6L))
  fit <- combat_ref_ql(z$counts, z$batch, reference = "A", verbose = FALSE)
  confidence <- fit@diagnostics$correction_confidence
  expect_true(all(is.finite(confidence$confidence_score)))
  expect_true(all(confidence$confidence_score >= 0 & confidence$confidence_score <= 1))
  minimum <- tapply(confidence$confidence_score, confidence$gene, min)
  expect_equal(fit@gene_status$confidence_score[match(names(minimum), fit@gene_status$gene)],
    as.numeric(minimum))
  expect_equal(combatrefql:::confidence_label(c(0.749, 0.75, 0.499, 0.5)),
    c("moderate", "high", "low", "moderate"))
  expect_equal(combatrefql:::confidence_label(1, FALSE), "failed")
  score <- combatrefql:::confidence_score
  expect_gte(score(0.2, 10, 6, TRUE, TRUE), score(1, 10, 6, TRUE, TRUE))
  expect_gte(score(0.2, 10, 6, TRUE, TRUE), score(0.2, 1, 6, TRUE, TRUE))
  expect_equal(score(0.2, 10, 6, TRUE, FALSE), 0)
})
