make_gene_status <- function(counts, batch, keep, mapping_valid) {
  zero_support <- vapply(levels(batch), function(x) rowSums(counts[, batch == x, drop = FALSE]) == 0, logical(nrow(counts)))
  if (is.null(dim(zero_support))) zero_support <- matrix(zero_support, ncol = 1L)
  zero_labels <- apply(zero_support, 1L, function(x) paste(levels(batch)[x], collapse = ","))
  status <- rep("unsupported", nrow(counts))
  status[keep] <- ifelse(mapping_valid, "adjusted", "failed")
  reason <- rep("all zero in at least one batch", nrow(counts))
  reason[status == "adjusted"] <- NA_character_
  reason[status == "failed"] <- "negative-binomial transport failed; original row retained"
  data.frame(gene = rownames(counts), status = status, reason = reason,
    adjusted = status == "adjusted", zero_batches = zero_labels,
    fit_valid = keep, mapping_valid = keep & status == "adjusted",
    stringsAsFactors = FALSE)
}

confidence_label <- function(score, mapping_valid = TRUE) {
  mapping_valid <- rep_len(mapping_valid, length(score))
  ifelse(!mapping_valid, "failed",
    ifelse(score >= 0.75, "high", ifelse(score >= 0.5, "moderate", "low")))
}

confidence_score <- function(ql_se, effective_information, samples,
                             dispersion_estimable, mapping_valid) {
  batch_cap <- ifelse(samples == 1L, 0.25, ifelse(samples == 2L, 0.5, 1))
  score <- pmin(1 / (1 + pmax(ql_se, 0)),
    pmax(effective_information, 0) / (1 + pmax(effective_information, 0)),
    batch_cap, ifelse(dispersion_estimable, 1, 0.25))
  score[!mapping_valid] <- 0
  score
}

make_correction_confidence <- function(shrinkage, dispersion, transported,
                                       batch, reference) {
  contrast <- shrinkage$diagnostics
  key <- paste(dispersion$diagnostics$gene, dispersion$diagnostics$batch)
  matched <- match(paste(contrast$gene, contrast$source_batch), key)
  info <- dispersion$diagnostics$weighted_residual_information[matched]
  estimable <- dispersion$diagnostics$estimable[matched]
  failures <- paste(transported$failures$gene, transported$failures$source_batch)
  mapping_valid <- !paste(contrast$gene, contrast$source_batch) %in% failures & contrast$valid
  samples <- as.integer(table(batch)[contrast$source_batch])
  score <- confidence_score(contrast$ql_se, info, samples, estimable,
    mapping_valid)
  data.frame(gene = contrast$gene, source_batch = contrast$source_batch,
    reference_batch = reference, confidence_score = score,
    confidence_label = confidence_label(score, mapping_valid),
    ql_se = contrast$ql_se, shrinkage_weight = contrast$shrinkage_weight,
    adaptive_weight = contrast$adaptive_weight,
    effective_information = info, dispersion_estimable = estimable,
    mapping_valid = mapping_valid,
    status = ifelse(!mapping_valid, "failed", ifelse(samples <= 2L,
      "low_information", ifelse(estimable, "ok", "dispersion_fallback"))),
    reason = ifelse(!mapping_valid, "mapping failed",
      ifelse(samples == 1L, "singleton source batch",
        ifelse(samples == 2L, "two-sample source batch",
          ifelse(!estimable, "source dispersion not independently estimable", NA_character_)))),
    stringsAsFactors = FALSE)
}

make_batch_quality <- function(selected, dispersion, batch) {
  sizes <- table(batch)
  rows <- lapply(names(sizes), function(value) {
    d <- dispersion$diagnostics[dispersion$diagnostics$batch == value, ]
    n <- as.integer(sizes[value])
    status <- if (n == 1L) "low_confidence" else if (n == 2L) "limited_information" else
      if (all(d$estimable)) "ok" else "dispersion_fallback"
    data.frame(batch = value, samples = n,
      effective_information = stats::median(d$weighted_residual_information),
      reference_eligible = n > 1L,
      selected_reference = value == selected$reference,
      confidence = if (n <= 2L) "low" else if (all(d$estimable)) "standard" else "low",
      status = status,
      message = switch(status, low_confidence = "Singleton: dispersion is shared/fallback and corrections are low-confidence.",
        limited_information = "Two samples provide limited independent information.",
        dispersion_fallback = "Batch-specific dispersion was not independently estimable.",
        "No batch-level limitation detected."), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

make_outcomes <- function(gene_status) {
  levels <- c("adjusted", "unsupported", "failed")
  genes <- tabulate(match(gene_status$status, levels), nbins = length(levels))
  data.frame(outcome = levels, genes = genes,
    proportion = genes / nrow(gene_status), stringsAsFactors = FALSE)
}

build_diagnostics <- function(selected, transported, timings, prepared, design,
                              warnings, gene_status, ql, shrinkage, dispersion,
                              confidence, batch_quality) {
  list(
    batches = selected$batches,
    mapping = transported$mapping,
    mapping_failures = transported$failures,
    timing = timings,
    input = cbind(prepared$input, design_rank = design$rank,
      coefficients = ncol(design$matrix), residual_df = design$residual_df,
      condition_number = design$condition_number,
      maximum_entanglement = design$maximum_entanglement,
      entanglement_warning = design$maximum_entanglement >= 0.8),
    design_entanglement = design$entanglement,
    batch_quality = batch_quality,
    correction_confidence = confidence,
    input_actions = prepared$input_actions,
    warnings = warnings,
    outcomes = make_outcomes(gene_status),
    ql = ql,
    batch_contrasts = shrinkage$diagnostics,
    dispersion = dispersion$diagnostics
  )
}

empty_warnings <- function() data.frame(stage = character(), class = character(), message = character())

timing_row <- function(stage, start) data.frame(stage = stage, seconds = proc.time()[["elapsed"]] - start)

`%||%` <- function(x, y) if (is.null(x)) y else x
