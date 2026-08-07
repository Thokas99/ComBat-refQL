# ComBat-refQL

[![R-CMD-check](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/ComBat-refQL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Thokas99/ComBat-refQL/actions/workflows/pkgdown.yaml/badge.svg)](https://thokas99.github.io/ComBat-refQL/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Reference-batch adjustment for bulk RNA-seq count data.**

ComBat-refQL models batch effects with edgeR quasi-likelihood methods,
moderates batch contrasts with empirical Bayes, and maps source batches toward
an observed reference batch while preserving modelled biological structure.

[Documentation](https://thokas99.github.io/ComBat-refQL/) ·
[Getting started](https://thokas99.github.io/ComBat-refQL/articles/getting-started.html) ·
[Method](https://thokas99.github.io/ComBat-refQL/articles/method.html) ·
[Validation](https://thokas99.github.io/ComBat-refQL/articles/validation.html)

## Highlights

- **Reference-based adjustment** — maps source batches toward one observed
  batch.
- **Biology-aware modelling** — preserves supplied biological groups and
  covariates in the design.
- **Count-scale output** — returns integer adjusted counts with original
  dimensions and identifiers.
- **Adaptive moderation** — combines quasi-likelihood evidence with
  empirical-Bayes stabilization.
- **Diagnostics** — reports batch quality, correction confidence, reference
  selection, and gene outcomes.
- **Reproducible fit objects** — retains detailed model and diagnostic
  information for inspection.

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
  group = metadata$condition
)

fit
summary(fit)
adjusted_counts <- fit@counts
fit@gene_status
```

Set `verbose = FALSE` to suppress fitting messages. Supply
`reference = "batch_name"` when the target batch is known; otherwise the
package selects a reference from biologically adjusted within-batch fits.

> **Input**
>
> Use raw or count-scale RNA-seq data with genes in rows and samples in
> columns. Batch labels and optional biological metadata must follow the
> sample order in `counts`.

## Choosing the right downstream workflow

The most useful representation depends on the scientific question.

> **Differential-expression inference:** use the original counts and include
> batch directly in the statistical design. This retains the observed count
> observations and propagates batch structure through the inferential model.
>
> **PCA, clustering, visualization, machine learning, and integrative
> analyses:** use ComBat-refQL adjusted counts when a harmonized expression
> space referenced to an observed batch is useful.

| Downstream task | Recommended data |
|---|---|
| Differential expression | Original counts + batch in the statistical model |
| PCA / UMAP / visualization | ComBat-refQL adjusted counts |
| Clustering and heatmaps | ComBat-refQL adjusted counts |
| Machine learning and classification | ComBat-refQL adjusted counts, with leakage-aware preprocessing |
| Cross-batch exploratory analysis | ComBat-refQL adjusted counts |

For differential expression, include batch directly in an edgeR, limma-voom,
DESeq2, or equivalent design. For predictive modelling, perform data-dependent
preprocessing within the training workflow and apply the learned procedure to
held-out data so evaluation samples remain independent of preprocessing
decisions.

## Working with adjusted counts

Adjusted counts are particularly useful when downstream methods operate
directly on an expression matrix and benefit from samples represented in a
common batch-adjusted space.

```r
log_cpm <- edgeR::cpm(fit@counts, log = TRUE, prior.count = 2)
pca <- prcomp(t(log_cpm), scale. = FALSE)
```

The selected reference-batch counts are preserved exactly. Source batches are
mapped toward their reference-counterfactual distribution, and the returned
matrix keeps the input dimensions and identifiers.

## Inputs and diagnostics

`combat_ref_ql()` accepts a raw count matrix, a discrete `batch` label, and
optional biological `group`, sample-level `covariates`, and explicit
`reference`. The package validates counts, identifiers, metadata alignment,
and design identifiability before fitting.

The fitted object records input actions, reference selection, batch-quality
diagnostics, correction confidence, gene outcomes, batch contrasts, dispersion
information, and full timing data. See [Input sanity and correction
confidence](https://thokas99.github.io/ComBat-refQL/articles/input-sanity-and-confidence.html)
for interpretation of replication, design association, and confidence labels.

## Reference selection

Choose `reference` explicitly when a target batch is scientifically preferred.
When it is omitted, ComBat-refQL selects a reference using biologically
adjusted within-batch modelling. The selected batch is reported in the fit and
its observed counts remain unchanged.

## Method and validation

ComBat-refQL combines quasi-likelihood modelling, hierarchical dispersion
modelling, empirical-Bayes moderation, and deterministic negative-binomial
count transport. The [method article](https://thokas99.github.io/ComBat-refQL/articles/method.html)
describes the estimand and fitting steps.

The validation suite evaluates cross-batch alignment on GFRN and NASA GeneLab
data, checks prespecified regression anchors, and requires exact reference
invariance. These are empirical diagnostics; the [validation
article](https://thokas99.github.io/ComBat-refQL/articles/validation.html)
provides the results and reproducible workflow.

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

## Getting help

Please open a [GitHub issue](https://github.com/Thokas99/ComBat-refQL/issues)
with a minimal reproducible example when you find a problem.

## License

ComBat-refQL is MIT licensed. The original ComBat-ref MIT notice is retained
in [LICENSES/ComBat-ref-MIT.txt](LICENSES/ComBat-ref-MIT.txt).
