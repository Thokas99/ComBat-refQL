# combatrefql 0.0.3

- Simplified `combat_ref_ql()` by internalizing the fixed TMMwsp normalization
  and replacing multi-level `verbosity` with logical `verbose`.
- Streamlined fitting output with concise `cli` progress, retained total
  runtime and important input actions, and simplified fit and summary printing.
- Refreshed documentation, vignettes, tests, and pkgdown for the public API.

# combatrefql 0.0.2

- Added structured input, design-entanglement, batch-quality, per-correction,
  and conservative per-gene confidence diagnostics without changing the fitted
  statistical formulation or count transport.
- Singleton source batches now use the existing shared-information fallback,
  are excluded from automatic reference selection, and are labelled
  low-confidence; two-sample batches are labelled limited-information.
- Added explicit warnings for strong batch/biology association and rank-based
  errors for exact aliasing or non-identifiable designs.
- Added first-class integration coverage for calls without groups or covariates.
- Added the input-sanity and correction-confidence vignette, concise confidence
  summaries, and pkgdown navigation for all guides.

# combatrefql 0.0.1

- Initial public implementation of ComBat-refQL.
- QL-based reference-batch modelling with evidence-adaptive EB moderation.
- Hierarchical batch-specific negative-binomial dispersion.
- Deterministic mid-P count transport with exact reference preservation.
- Validated S7 fit and summary objects with typed conditions and tidy diagnostic
  tables.
- Audited fractional-count rounding, deterministic chunked transport, and
  biologically adjusted automatic reference scoring.
- Reproducible real-data validation on GFRN and NASA GeneLab.
- Added the getting-started, method, and validation vignettes, package website,
  README installation guidance, citation metadata, and GitHub Pages deployment.
