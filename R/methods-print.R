#' @export
S7::method(print, CombatRefQLFit) <- function(x, ...) {
  outcomes <- stats::setNames(x@diagnostics$outcomes$genes, x@diagnostics$outcomes$outcome)
  total <- x@diagnostics$timing$seconds[x@diagnostics$timing$stage == "total"]
  groups <- length(unique(stats::na.omit(x@samples$group)))
  cli::cli_h1("ComBat-refQL fit")
  cli::cli_h2("Data")
  cli::cli_text("{format(nrow(x@counts), big.mark = '')} genes \u00d7 {ncol(x@counts)} samples")
  cli::cli_text("{length(unique(x@samples$batch))} batches")
  cli::cli_text("Biological groups  {if (groups) groups else 'none'}")
  cli::cli_h2("Reference")
  cli::cli_text("{x@specification$reference_batch} ({x@specification$reference_selection})")
  cli::cli_h2("Outcome")
  labels <- c(adjusted = "Adjusted", unsupported = "Unchanged", failed = "Failed")
  for (key in names(labels)) {
    genes <- outcomes[[key]] %||% 0
    proportion <- if (nrow(x@counts)) 100 * genes / nrow(x@counts) else 0
    suffix <- if (key == "failed" && genes == 0) "" else
      sprintf(" (%.1f%%)", proportion)
    cli::cli_text("{labels[[key]]}  {genes}{suffix}")
  }
  cli::cli_h2("Confidence")
  confidence <- table(factor(x@diagnostics$correction_confidence$confidence_label,
    levels = c("high", "moderate", "low", "failed")))
  for (label in names(confidence)) cli::cli_text("{label}  {confidence[[label]]}")
  if (nrow(x@diagnostics$warnings)) {
    cli::cli_h2("Warnings")
    cli::cli_text("{nrow(x@diagnostics$warnings)} stored diagnostic warnings")
  }
  cli::cli_h2("Runtime")
  cli::cli_text("{format(round(total, 1), nsmall = 1)} s")
  invisible(x)
}

#' @export
S7::method(print, CombatRefQLSummary) <- function(x, ...) {
  cli::cli_h1("ComBat-refQL summary")
  cli::cli_h2("Design")
  design <- x@design[1L, ]
  cli::cli_text("Samples  {design$samples}")
  cli::cli_text("Coefficients  {design$coefficients}; rank {design$rank}")
  cli::cli_text("Residual df  {design$residual_df}")
  cli::cli_text("Condition number  {format(round(design$condition_number, 2), nsmall = 2)}")
  cli::cli_text("Biological groups  {if (x@overview$groups) x@overview$groups else 'none'}")
  cli::cli_text("Maximum batch association  {format(round(x@overview$maximum_entanglement, 3), nsmall = 3)}")
  if (x@overview$low_information_batches)
    cli::cli_alert_warning("{x@overview$low_information_batches} low-information batch(es)")
  cli::cli_h2("Reference")
  selected <- x@batches$batch[x@batches$selected][1L]
  cli::cli_text("Selected  {selected} ({x@specification$reference_selection})")
  for (i in seq_len(nrow(x@batches))) {
    score <- if (is.finite(x@batches$reference_score[i]))
      format(signif(x@batches$reference_score[i], 4)) else "not scored"
    cli::cli_text("{x@batches$batch[i]}  {x@batches$samples[i]} samples; score {score}")
  }
  if (x@specification$reference_selection == "automatic") {
    runner_up <- x@batches$batch[x@batches$second_best]
    cli::cli_text("Runner-up  {runner_up}; score margin {format(signif(x@batches$score_margin[1L], 4))}")
  }
  cli::cli_h2("Batch adjustment")
  for (i in seq_len(nrow(x@dispersion))) {
    if (x@dispersion$batch[i] == selected) next
    cli::cli_text("{x@dispersion$batch[i]}  dispersion x{format(signif(x@dispersion$multiplier[i], 3))} ({x@dispersion$status[i]})")
  }
  for (i in seq_len(nrow(x@shrinkage))) cli::cli_text(
    "{x@shrinkage$source_batch[i]}  median adaptive weight {format(signif(x@shrinkage$median_weight[i], 3))}")
  cli::cli_h2("Gene outcomes")
  labels <- c(adjusted = "Adjusted", unsupported = "Unchanged", failed = "Failed")
  for (i in seq_len(nrow(x@gene_outcomes))) cli::cli_text(
    "{labels[[x@gene_outcomes$outcome[i]]]}  {x@gene_outcomes$genes[i]} ({format(round(100*x@gene_outcomes$proportion[i], 1), nsmall=1)}%)")
  cli::cli_h2("Correction confidence")
  for (label in names(x@confidence$genes)) cli::cli_text("{label}  {x@confidence$genes[[label]]}")
  actions <- x@input_actions[x@input_actions$lossy, , drop = FALSE]
  if (nrow(actions)) {
    cli::cli_h2("Input actions")
    for (message in actions$message) cli::cli_alert_info(message)
  }
  if (nrow(x@warnings)) {
    cli::cli_h2("Warnings")
    cli::cli_text("{nrow(x@warnings)} stored diagnostic warnings")
  }
  cli::cli_h2("Runtime")
  total <- x@timings$seconds[x@timings$stage == "total"]
  cli::cli_text("{format(round(total, 1), nsmall = 1)} s")
  invisible(x)
}
