download_cached <- function(url, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) utils::download.file(url, path, mode = "wb", quiet = TRUE)
  if (!file.exists(path) || file.info(path)$size == 0)
    stop("Download failed: ", url)
  path
}

log_cpm <- function(counts) edgeR::cpm(counts, log = TRUE, prior.count = 1)

alignment_metrics <- function(counts, batch, group) {
  y <- log_cpm(counts)
  design <- stats::model.matrix(~ group + batch)
  coefficients <- qr.solve(design, t(y))
  batch_columns <- grep("^batch", colnames(design))
  batch_effect <- coefficients[batch_columns, , drop = FALSE]
  residuals <- t(y) - design %*% coefficients
  residual_variance <- vapply(levels(batch), function(level)
    mean(apply(residuals[batch == level, , drop = FALSE], 2, stats::var)),
    numeric(1))

  top <- order(apply(y, 1, stats::var), decreasing = TRUE)[
    seq_len(min(500L, nrow(y)))]
  pca <- stats::prcomp(t(y[top, , drop = FALSE]))
  scores <- pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE]
  variance <- pca$sdev[seq_len(ncol(scores))]^2
  incremental <- vapply(seq_len(ncol(scores)), function(i) {
    reduced <- stats::lm(scores[, i] ~ group)
    full <- stats::lm(scores[, i] ~ group + batch)
    max(0, summary(full)$adj.r.squared - summary(reduced)$adj.r.squared)
  }, numeric(1))

  c(
    conditional_batch_r2 = stats::weighted.mean(incremental, variance),
    mean_alignment_rmse = sqrt(mean(batch_effect^2)),
    descriptive_residual_logcpm_variance_alignment =
      stats::sd(log(residual_variance))
  )
}

run_validation <- function(counts, batch, group, dataset) {
  batch <- droplevels(factor(batch))
  group <- droplevels(factor(group))
  keep <- rowSums(counts) > 0 & vapply(seq_len(nrow(counts)), function(i)
    all(vapply(levels(batch), function(level)
      any(counts[i, batch == level] > 3), logical(1))), logical(1))
  counts <- unname(as.matrix(counts[keep, , drop = FALSE]))
  storage.mode(counts) <- "double"
  dimnames(counts) <- list(paste0("g", seq_len(nrow(counts))),
    paste0("s", seq_len(ncol(counts))))

  fit <- combat_ref_ql(counts, batch, group, verbosity = "quiet")
  reference <- fit@specification$reference_batch
  outcomes <- stats::setNames(fit@diagnostics$outcomes$genes,
    fit@diagnostics$outcomes$outcome)
  raw_metrics <- alignment_metrics(counts, batch, group)
  data.frame(
    dataset = dataset,
    genes = nrow(counts), samples = ncol(counts), reference_batch = reference,
    adjusted_genes = outcomes[["adjusted"]],
    unchanged_genes = outcomes[["unsupported"]],
    failed_genes = outcomes[["failed"]],
    as.list(stats::setNames(raw_metrics, paste0("raw_", names(raw_metrics)))),
    as.list(alignment_metrics(fit@counts, batch, group)),
    transformation_magnitude = mean(abs(fit@counts - counts)),
    reference_invariant = identical(
      fit@counts[, batch == reference, drop = FALSE],
      counts[, batch == reference, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}
