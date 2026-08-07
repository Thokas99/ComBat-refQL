#' Fit ComBat-refQL reference-batch adjustment
#'
#' Fits model-based reference-counterfactual means with edgeR v4, moderates QL
#' batch contrasts with continuous evidence-adaptive EB, and adjusts
#' non-reference counts by deterministic negative-binomial mid-P transport.
#'
#' @param counts Numeric raw-count matrix with genes in rows and samples in columns.
#' @param batch Sample batch vector, optionally named for automatic reordering.
#' @param group Optional biological group vector to preserve. Exact batch/group
#'   aliasing is an error; strong association is reported as a warning.
#' @param covariates Optional sample-by-covariate data frame or matrix. Constant
#'   columns are recorded and omitted; redundant or aliased terms are errors.
#' @param reference Optional reference batch. When `NULL`, the minimum-dispersion
#'   batch is selected from biologically adjusted within-batch QL fits.
#' @param fractional_counts Round with an audit record, or reject fractional values.
#' @param chunk_size Positive number of genes transported per chunk.
#' @param verbose Logical. If `TRUE`, print concise progress, outcome, and
#'   runtime messages. Warnings and errors are not suppressed when `FALSE`.
#' @return A validated [CombatRefQLFit] object with adjusted counts and structured
#'   input, batch-quality, and correction-confidence diagnostics.
#' @examples
#' set.seed(1)
#' batch <- factor(rep(c("reference", "study"), each = 6))
#' group <- factor(rep(c("control", "treated"), 6))
#' counts <- matrix(rnbinom(1200, mu = 40, size = 8), nrow = 100)
#' rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
#' colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
#' counts[, batch == "study"] <- counts[, batch == "study"] * 2
#'
#' fit <- combat_ref_ql(counts, batch, group, reference = "reference")
#' fit
#' summary(fit)
#' corrected <- fit@counts
#' @export
combat_ref_ql <- function(counts, batch, group = NULL, covariates = NULL,
                          reference = NULL,
                          fractional_counts = c("round", "error"), chunk_size = 5000L,
                          verbose = TRUE) {
  call <- match.call()
  fractional_counts <- match.arg(fractional_counts)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1 ||
      chunk_size != as.integer(chunk_size))
    input_error("{.arg chunk_size} must be one positive integer.")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    input_error("{.arg verbose} must be {.val TRUE} or {.val FALSE}.")
  chunk_size <- as.integer(chunk_size)
  reporter <- new_reporter(verbose)
  reporter$start()
  total_start <- proc.time()[["elapsed"]]
  start <- proc.time()[["elapsed"]]
  prepared <- prepare_inputs(counts, batch, group, covariates, reference, fractional_counts)
  timings <- timing_row("preparation", start)
  if (prepared$input_actions$affected_values[[1L]])
    reporter$info(prepared$input_actions$message[[1L]])
  reporter$success(sprintf("Prepared %s genes \u00d7 %s samples",
    nrow(prepared$counts), ncol(prepared$counts)))

  start <- proc.time()[["elapsed"]]
  dge <- normalize_counts(prepared, "TMMwsp")
  timings <- rbind(timings, timing_row("normalization", start))

  design <- build_design(prepared$batch, prepared$group, prepared$covariates)
  if (nrow(design$warnings)) warn_combatrefql("combatrefql_design_warning",
    design$warnings$message[[1L]], stage = "design",
    maximum_entanglement = design$maximum_entanglement)
  keep <- assess_support(prepared$counts, prepared$batch)
  start <- proc.time()[["elapsed"]]
  selected <- select_reference(dge, prepared$batch, prepared$group,
    prepared$covariates, prepared$reference)
  timings <- rbind(timings, timing_row("reference", start))
  reference_label <- if (is.null(prepared$reference)) selected$reference else
    sprintf("%s (explicit)", selected$reference)
  reporter$success(sprintf("Reference: %s", reference_label))

  start <- proc.time()[["elapsed"]]
  supported <- dge[keep, , keep.lib.sizes = TRUE]
  preliminary <- fit_ql_model(supported, design$matrix)
  timings <- rbind(timings, timing_row("preliminary_ql_fit", start))

  start <- proc.time()[["elapsed"]]
  preliminary_phi <- matrix(as.numeric(preliminary$dispersion), sum(keep),
    ncol(prepared$counts))
  information <- effective_residual_information(design$matrix, prepared$batch,
    preliminary$fitted.values, preliminary_phi, preliminary$weights)
  dispersion <- estimate_batch_dispersion(supported, design$matrix,
    prepared$batch, selected$reference, preliminary, information)
  timings <- rbind(timings, timing_row("dispersion", start))

  start <- proc.time()[["elapsed"]]
  fit <- fit_ql_model(supported, design$matrix, dispersion$dispersion)
  timings <- rbind(timings, timing_row("final_ql_fit", start))

  start <- proc.time()[["elapsed"]]
  shrinkage <- shrink_batch_effects(fit, dispersion$dispersion,
    prepared$batch, selected$reference)
  timings <- rbind(timings, timing_row("contrast_moderation", start))
  reporter$success("Model fitted")

  start <- proc.time()[["elapsed"]]
  transported <- transport_counts(prepared$counts, keep, prepared$batch,
    selected$reference, fit, shrinkage, dispersion,
    chunk_size)
  timings <- rbind(timings, timing_row("mapping", start))
  timings <- rbind(timings, data.frame(stage = "total", seconds = proc.time()[["elapsed"]] - total_start))
  built <- build_combatrefql_result(prepared, dge, design, selected, fit,
    shrinkage, dispersion, transported, keep, timings, call, "TMMwsp",
    chunk_size)
  if (nrow(built$warnings)) reporter$warning(built$warnings$message[[1L]])
  outcomes <- stats::setNames(built$result@diagnostics$outcomes$genes,
    built$result@diagnostics$outcomes$outcome)
  reporter$outcome(outcomes[["adjusted"]], outcomes[["unsupported"]],
    outcomes[["failed"]])
  reporter$finish(timings$seconds[timings$stage == "total"])
  built$result
}
