#' Fit ComBat-refQL reference-batch adjustment
#'
#' Fits model-based reference-counterfactual means with edgeR v4, moderates QL
#' batch contrasts with continuous evidence-adaptive EB, and adjusts
#' non-reference counts by deterministic negative-binomial mid-P transport.
#'
#' @param counts Numeric raw-count matrix with genes in rows and samples in columns.
#' @param batch Sample batch vector, optionally named for automatic reordering.
#' @param group Optional biological group vector.
#' @param covariates Optional sample-by-covariate data frame or matrix.
#' @param reference Optional reference batch. When `NULL`, the minimum-dispersion
#'   batch is selected from biologically adjusted within-batch QL fits.
#' @param normalization edgeR library normalization method. Only `"TMMwsp"` is supported.
#' @param fractional_counts Round with an audit record, or reject fractional values.
#' @param chunk_size Positive number of genes transported per chunk.
#' @param verbosity User output level.
#' @return A validated [CombatRefQLFit] object.
#' @examples
#' set.seed(1)
#' batch <- factor(rep(c("reference", "study"), each = 6))
#' group <- factor(rep(c("control", "treated"), 6))
#' counts <- matrix(rnbinom(1200, mu = 40, size = 8), nrow = 100)
#' rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
#' colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
#' counts[, batch == "study"] <- counts[, batch == "study"] * 2
#'
#' fit <- combat_ref_ql(counts, batch, group, reference = "reference",
#'   verbosity = "quiet")
#' fit
#' summary(fit)
#' corrected <- fit@counts
#' @export
combat_ref_ql <- function(counts, batch, group = NULL, covariates = NULL,
                          reference = NULL,
                          normalization = "TMMwsp",
                          fractional_counts = c("round", "error"), chunk_size = 5000L,
                          verbosity = c("normal", "quiet", "verbose")) {
  call <- match.call()
  fractional_counts <- match.arg(fractional_counts)
  verbosity <- match.arg(verbosity)
  if (!identical(normalization, "TMMwsp"))
    input_error("{.arg normalization} must be {.val TMMwsp}.")
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1 ||
      chunk_size != as.integer(chunk_size))
    input_error("{.arg chunk_size} must be one positive integer.")
  chunk_size <- as.integer(chunk_size)
  reporter <- new_reporter(verbosity)
  total_start <- proc.time()[["elapsed"]]
  start <- proc.time()[["elapsed"]]
  prepared <- prepare_inputs(counts, batch, group, covariates, reference, fractional_counts)
  timings <- timing_row("preparation", start)
  if (prepared$input_actions$affected_values[[1L]])
    reporter$info(prepared$input_actions$message[[1L]])
  reporter$success(sprintf("Prepared %s genes \u00d7 %s samples",
    nrow(prepared$counts), ncol(prepared$counts)))

  start <- proc.time()[["elapsed"]]
  dge <- normalize_counts(prepared, normalization)
  timings <- rbind(timings, timing_row("normalization", start))

  design <- build_design(prepared$batch, prepared$group, prepared$covariates)
  keep <- assess_support(prepared$counts, prepared$batch)
  start <- proc.time()[["elapsed"]]
  selected <- select_reference(dge, prepared$batch, prepared$group,
    prepared$covariates, prepared$reference, reporter$reference)
  timings <- rbind(timings, timing_row("reference", start))
  reporter$success(sprintf("Selected reference: %s", selected$reference))

  start <- proc.time()[["elapsed"]]
  supported <- dge[keep, , keep.lib.sizes = TRUE]
  preliminary <- fit_ql_model(supported, design$matrix)
  timings <- rbind(timings, timing_row("preliminary_ql_fit", start))
  reporter$success("Fitted preliminary model")

  start <- proc.time()[["elapsed"]]
  preliminary_phi <- matrix(as.numeric(preliminary$dispersion), sum(keep),
    ncol(prepared$counts))
  information <- effective_residual_information(design$matrix, prepared$batch,
    preliminary$fitted.values, preliminary_phi, preliminary$weights)
  dispersion <- estimate_batch_dispersion(supported, design$matrix,
    prepared$batch, selected$reference, preliminary, information)
  timings <- rbind(timings, timing_row("dispersion", start))
  reporter$success("Estimated hierarchical batch-specific NB dispersion")

  start <- proc.time()[["elapsed"]]
  fit <- fit_ql_model(supported, design$matrix, dispersion$dispersion)
  timings <- rbind(timings, timing_row("final_ql_fit", start))
  reporter$success("Fitted final model")

  start <- proc.time()[["elapsed"]]
  shrinkage <- shrink_batch_effects(fit, dispersion$dispersion,
    prepared$batch, selected$reference)
  timings <- rbind(timings, timing_row("contrast_moderation", start))

  start <- proc.time()[["elapsed"]]
  mapping_steps <- length(setdiff(levels(prepared$batch), selected$reference)) * ceiling(sum(keep) / chunk_size)
  reporter$start_mapping(mapping_steps)
  transported <- transport_counts(prepared$counts, keep, prepared$batch,
    selected$reference, fit, shrinkage, dispersion,
    chunk_size, reporter$mapping)
  timings <- rbind(timings, timing_row("mapping", start))
  timings <- rbind(timings, data.frame(stage = "total", seconds = proc.time()[["elapsed"]] - total_start))
  built <- build_combatrefql_result(prepared, dge, design, selected, fit,
    shrinkage, dispersion, transported, keep, timings, call, normalization,
    chunk_size)
  if (nrow(built$warnings)) reporter$warning(built$warnings$message[[1L]])
  outcomes <- stats::setNames(built$result@diagnostics$outcomes$genes,
    built$result@diagnostics$outcomes$outcome)
  reporter$outcome(outcomes[["adjusted"]], outcomes[["unsupported"]],
    outcomes[["failed"]])
  reporter$finish(timings$seconds[timings$stage == "total"])
  built$result
}
