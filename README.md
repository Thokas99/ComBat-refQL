# ComBat-refQL

ComBat-refQL adjusts RNA-seq count matrices toward one observed reference
batch. It combines edgeR quasi-likelihood modelling, evidence-adaptive
empirical-Bayes moderation, hierarchical batch-specific negative-binomial
dispersion, and deterministic count transport.

## Installation

```r
# install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")
```

ComBat-refQL requires R 4.3 or later and edgeR 4.2 or later.

## Quick start

```r
library(combatrefql)

fit <- combat_ref_ql(
  counts,
  batch = metadata$batch,
  group = metadata$condition,
  reference = "reference"
)

corrected <- fit@counts
summary(fit)
```

Genes are rows and samples are columns. Corrected counts are finite,
non-negative integer values with the original dimensions and dimnames.

## Method

The method fits a QL mean model, estimates hierarchical reference-centred NB
dispersion, and moderates each source-versus-reference batch contrast. Weak or
uncertain contrasts remain EB-stabilized; contrasts strongly separated from
the fitted prior move continuously toward raw QL. Counts are mapped to the
reference-counterfactual distribution by deterministic NB mid-P transport.

## Reference batch

Set `reference` when the target batch is known. Otherwise ComBat-refQL selects
a reference from biologically adjusted within-batch fits. Observed counts in
the selected reference batch are preserved exactly.

## Validation

The final implementation reproduces prespecified GFRN and NASA GeneLab
alignment metrics with exact reference invariance. These are empirical
diagnostics, not real-data ground truth. See [validation](vignettes/validation.Rmd).

## Documentation

- [Getting started](vignettes/getting-started.Rmd)
- [Method](vignettes/method.Rmd)
- [Function reference](man/combat_ref_ql.Rd)

## Citation and attribution

ComBat-refQL extends the reference-batch framework introduced by Xiaoyu Zhang
in *Highly effective batch effect correction method for RNA-seq count data*
([doi:10.1016/j.csbj.2024.12.010](https://doi.org/10.1016/j.csbj.2024.12.010)).
Run `citation("combatrefql")` and `citation("edgeR")` for citation metadata.

## License

ComBat-refQL is MIT licensed. The original ComBat-ref MIT notice is retained in
[LICENSES/ComBat-ref-MIT.txt](LICENSES/ComBat-ref-MIT.txt).
