# Real-data validation

The suite uses the archived ComBat-ref GFRN and NASA GeneLab count matrices.
GFRN retains 47 samples across six biological signatures. NASA retains 35 mouse
liver samples from OSD/GLDS 47, 48, 137, 173, 242, and 245; GLDS-168 is excluded
to match the ComBat-ref analysis. Condition is protected while study-mission is
treated as batch.

From the repository root, run:

```sh
Rscript validation/run_all.R
```

The runner downloads commit-pinned inputs into ignored `validation/.cache/`,
measures raw and adjusted alignment, checks adjusted regression anchors and
exact reference invariance, then writes:

- `real_data_metrics.tsv`: retained dimensions, selected reference, outcomes,
  alignment diagnostics, transformation magnitude, and reference invariance;
- `session-info.txt`: package, source commit, R, platform, and dependency state.

For conditional batch R2, mean-alignment RMSE, and residual logCPM variance
alignment, smaller is better. Raw columns make the direction and magnitude of
change auditable; adjusted columns are the prespecified regression anchors.
The residual-variance score is not negative-binomial dispersion, and neither
dataset provides uniform biological ground truth.

Sources: [ComBat-ref](https://doi.org/10.1016/j.csbj.2024.12.010),
[ComBat-seq](https://doi.org/10.1093/nargab/lqaa078), and the
[NASA GeneLab benchmark](https://doi.org/10.3389/fspas.2023.1200132).
