build_design <- function(batch, group = NULL, covariates = NULL) {
  data <- data.frame(batch = batch)
  rhs <- "0 + batch"
  if (!is.null(group)) {
    data$group <- group
    rhs <- paste(rhs, "+ group")
  }
  if (!is.null(covariates)) {
    names(covariates) <- make.names(names(covariates), unique = TRUE)
    data <- cbind(data, covariates)
    rhs <- paste(rhs, "+", paste(names(covariates), collapse = " + "))
  }
  design <- stats::model.matrix(stats::as.formula(paste("~", rhs)), data)
  colnames(design)[seq_len(nlevels(batch))] <- paste0("batch", levels(batch))
  rank <- qr(design)$rank
  residual_df <- nrow(design) - rank
  if (rank < ncol(design) || residual_df <= 0L) {
    abort_combatrefql("combatrefql_design_error",
      "The batch, group, and covariate design is rank deficient or has no residual degrees of freedom.",
      "*" = "Remove confounded or redundant terms, or add replicates.",
      stage = "design", rank = rank, residual_df = residual_df)
  }
  batch_matrix <- stats::model.matrix(~ 0 + batch)
  assignment <- attr(design, "assign")
  labels <- attr(stats::terms(stats::as.formula(paste("~", rhs))), "term.labels")
  biology <- which(assignment > 1L)
  entanglement <- if (!length(biology)) data.frame(
    term = character(), columns = character(), batch_association = numeric(),
    status = character(), message = character()) else do.call(rbind,
    lapply(split(biology, assignment[biology]), function(index) {
      scores <- vapply(index, function(j) {
        y <- design[, j]
        fitted <- batch_matrix %*% stats::lm.fit(batch_matrix, y)$coefficients
        denominator <- sum((y - mean(y))^2)
        if (denominator == 0) 1 else max(0, min(1, 1 - sum((y - fitted)^2) / denominator))
      }, numeric(1))
      score <- max(scores)
      term <- labels[assignment[index[1L]]]
      strong <- score >= 0.8
      data.frame(term = term, columns = paste(colnames(design)[index], collapse = ","),
        batch_association = score, status = if (strong) "strong" else "ordinary",
        message = if (strong) sprintf("%s is strongly associated with batch (diagnostic R2 = %.3f).", term, score) else NA_character_,
        stringsAsFactors = FALSE)
    }))
  maximum <- if (nrow(entanglement)) max(entanglement$batch_association) else 0
  warnings <- if (maximum >= 0.8) data.frame(stage = "design",
    class = "combatrefql_design_warning",
    message = "Preserved biology/covariates are strongly associated with batch; corrections may be weakly identified.",
    stringsAsFactors = FALSE) else empty_warnings()
  list(matrix = design, rank = rank, residual_df = residual_df,
       condition_number = kappa(design), entanglement = entanglement,
       maximum_entanglement = maximum, warnings = warnings)
}

batch_contrast <- function(design, source, reference) {
  contrast <- stats::setNames(numeric(ncol(design)), colnames(design))
  source_column <- match(paste0("batch", source), names(contrast))
  reference_column <- match(paste0("batch", reference), names(contrast))
  if (!is.na(source_column)) contrast[source_column] <- 1
  if (!is.na(reference_column)) contrast[reference_column] <- -1
  if (is.na(source_column) && is.na(reference_column))
    abort_combatrefql("combatrefql_internal_error",
      "The fitted design does not identify the requested batch contrast.",
      stage = "contrast")
  contrast
}

effective_residual_information <- function(design, batch, fitted = NULL,
                                           dispersion = NULL, weights = NULL) {
  batches <- levels(batch)
  design_hat <- design %*% solve(crossprod(design), t(design))
  design_only <- vapply(batches, function(b)
    sum(1 - diag(design_hat)[batch == b]), numeric(1))
  if (is.null(fitted)) return(data.frame(batch = batches,
    design_residual_information = design_only))

  if (is.null(dim(dispersion))) dispersion <- matrix(dispersion,
    nrow(fitted), ncol(fitted))
  if (is.null(weights)) weights <- matrix(1, nrow(fitted), ncol(fitted))
  info <- matrix(0, nrow(fitted), length(batches),
    dimnames = list(rownames(fitted), batches))
  for (g in seq_len(nrow(fitted))) {
    working <- weights[g, ] * fitted[g, ] /
      (1 + dispersion[g, ] * fitted[g, ])
    root <- sqrt(pmax(working, 0))
    xw <- design * root
    inverse <- tryCatch(solve(crossprod(xw)), error = function(error) NULL)
    if (is.null(inverse)) next
    leverage <- rowSums((xw %*% inverse) * xw)
    info[g, ] <- vapply(batches, function(b)
      sum(1 - leverage[batch == b]), numeric(1))
  }
  info[info < sqrt(.Machine$double.eps)] <- 0
  precision <- ifelse(info > 0, 1 / trigamma(info / 2), 0)
  list(
    design_only = data.frame(batch = batches,
      design_residual_information = design_only),
    weighted = info,
    precision = precision,
    formula = "inverse_trigamma(weighted_residual_information / 2)"
  )
}
