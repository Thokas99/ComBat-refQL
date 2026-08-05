devtools::load_all(".", quiet = TRUE)
source("validation/real_data/common.R")
source("validation/real_data/gfrn/run.R")
source("validation/real_data/nasa/run.R")

metrics <- rbind(run_gfrn(), run_nasa())
anchors <- data.frame(
  dataset = c("GFRN", "NASA GeneLab"),
  mean_alignment_rmse = c(0.0812774, 0.1917808),
  conditional_batch_r2 = c(0.0311140, 0.1007872),
  descriptive_residual_logcpm_variance_alignment = c(0.1777319, 0.0929263)
)
matched <- anchors[match(metrics$dataset, anchors$dataset), ]
for (column in setdiff(names(anchors), "dataset"))
  stopifnot(isTRUE(all.equal(metrics[[column]], matched[[column]],
    tolerance = 5e-4)))
stopifnot(all(metrics$reference_invariant))

utils::write.table(metrics, "validation/real_data_metrics.tsv", sep = "\t",
  quote = FALSE, row.names = FALSE)
info <- capture.output({
  cat("ComBat-refQL version:", as.character(packageVersion("combatrefql")), "\n")
  cat("ComBat-ref source commit: f0c6d3f313e9f1f29de2c8c8f12b98113443be55\n")
  sessionInfo()
})
writeLines(info, "validation/session-info.txt")
print(metrics)
