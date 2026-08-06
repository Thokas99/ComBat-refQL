# ComBat-refQL

[![R-CMD-check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml/badge.svg)](https://thokas99.github.io/ComBat-refQL/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

ComBat-refQL adjusts RNA-seq count matrices toward one observed reference
batch. It combines edgeR quasi-likelihood modelling, evidence-adaptive
empirical-Bayes moderation, hierarchical batch-specific negative-binomial
dispersion, and deterministic count transport.

**[Read the documentation](https://thokas99.github.io/ComBat-refQL/)** for the
getting-started guide, method details, validation results, and function
reference.

## Overview

Use ComBat-refQL when RNA-seq count data contain a known batch effect and one
observed batch is an appropriate target. The package returns integer counts,
preserves the selected reference batch exactly, and keeps the original matrix
dimensions and dimnames.

## Installation

```r
# install.packages("pak")
pak::pak("Thokas99/ComBat-refQL")
```

ComBat-refQL requires R 4.3 or later and edgeR 4.2 or later.

## Usage

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

Genes are rows and samples are columns. `batch` and `group` must follow the
sample order in `counts`; `group` represents the biological condition that
should be retained during adjustment.

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

## Learn more

- [Website](https://thokas99.github.io/ComBat-refQL/)
- [Getting started](https://thokas99.github.io/ComBat-refQL/articles/getting-started.html)
- [Input sanity and correction confidence](https://thokas99.github.io/ComBat-refQL/articles/input-sanity-and-confidence.html)
- [Method](https://thokas99.github.io/ComBat-refQL/articles/method.html)
- [Function reference](https://thokas99.github.io/ComBat-refQL/reference/)

## Getting help

If you find a bug, please open a [GitHub issue](https://github.com/Thokas99/ComBat-refQL/issues)
with a minimal reproducible example.

## Citation and attribution

ComBat-refQL builds on these methods and validation resources:

- **ComBat-ref** — Xiaoyu Zhang, *Highly effective batch effect correction
  method for RNA-seq count data*: [article](https://doi.org/10.1016/j.csbj.2024.12.010)
  · [source repository](https://github.com/xiaoyu12/Combat-ref)
- **ComBat-seq** — Yuqing Zhang, Giovanni Parmigiani, and W. Evan Johnson,
  *ComBat-seq: batch effect adjustment for RNA-seq count data*:
  [article](https://doi.org/10.1093/nargab/lqaa078)
  · [source repository](https://github.com/zhangyuqing/ComBat-seq)
- **NASA GeneLab validation benchmark** — Sanders et al., *Batch effect
  correction methods for NASA GeneLab transcriptomic datasets*:
  [article](https://doi.org/10.3389/fspas.2023.1200132)

Run `citation("combatrefql")` and `citation("edgeR")` for complete metadata.

## License

ComBat-refQL is MIT licensed. The original ComBat-ref MIT notice is retained in
[LICENSES/ComBat-ref-MIT.txt](LICENSES/ComBat-ref-MIT.txt).
