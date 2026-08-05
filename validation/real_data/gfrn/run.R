run_gfrn <- function() {
  url <- paste0("https://raw.githubusercontent.com/xiaoyu12/Combat-ref/",
    "f0c6d3f313e9f1f29de2c8c8f12b98113443be55/",
    "real_data_application/signature_data.rds")
  object <- readRDS(download_cached(url,
    "validation/.cache/gfrn/signature_data.rds"))
  metadata <- as.data.frame(SummarizedExperiment::colData(object))
  group <- as.character(metadata$group)
  selected <- group %in% c("gfp_for_egfr", "gfp18", "gfp30", "her2",
    "egfr", "kraswt")
  group[group %in% c("gfp_for_egfr", "gfp18", "gfp30")] <- "gfp"
  run_validation(SummarizedExperiment::assay(object, "counts")[, selected,
    drop = FALSE], metadata$batch[selected], group[selected], "GFRN")
}
