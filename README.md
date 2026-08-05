# ComBat-refQL

[![R-CMD-check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml/badge.svg)](https://thokas99.github.io/ComBat-refQL/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

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

The real-data suite measures GFRN and NASA GeneLab alignment before and after
adjustment, checks prespecified regression anchors, and requires exact reference
invariance. These are empirical diagnostics, not ground truth. See the
[validation article](https://thokas99.github.io/ComBat-refQL/articles/validation.html)
or its [reproducible summary](validation/README.md).

## Documentation

- [Website](https://thokas99.github.io/ComBat-refQL/)
- [Getting started](https://thokas99.github.io/ComBat-refQL/articles/getting-started.html)
- [Method](https://thokas99.github.io/ComBat-refQL/articles/method.html)
- [Function reference](https://thokas99.github.io/ComBat-refQL/reference/)

## Citation and attribution

ComBat-refQL extends the reference-batch framework introduced by Xiaoyu Zhang
in *Highly effective batch effect correction method for RNA-seq count data*
([doi:10.1016/j.csbj.2024.12.010](https://doi.org/10.1016/j.csbj.2024.12.010)).
Its count-scale context includes ComBat-seq
([doi:10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078)); the NASA
validation follows the GeneLab benchmark
([doi:10.3389/fspas.2023.1200132](https://doi.org/10.3389/fspas.2023.1200132)).
Run `citation("combatrefql")` and `citation("edgeR")` for complete metadata.

## License

ComBat-refQL is MIT licensed. The original ComBat-ref MIT notice is retained in
[LICENSES/ComBat-ref-MIT.txt](LICENSES/ComBat-ref-MIT.txt).
