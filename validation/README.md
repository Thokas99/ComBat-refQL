# Real-data validation

The final ComBat-refQL method is validated on the established GFRN signature
dataset and NASA GeneLab liver collection. The sample definitions match the
original ComBat-ref workflow; NASA GLDS-168 is excluded as in that analysis.

From the repository root, run:

```sh
Rscript validation/run_all.R
```

The runner downloads inputs into ignored `validation/.cache/`, applies only the
final production method, checks prespecified regression anchors and exact
reference invariance, then writes:

- `real_data_metrics.tsv`: retained dimensions, selected reference, outcomes,
  alignment diagnostics, transformation magnitude, and reference invariance;
- `session-info.txt`: package, source commit, R, platform, and dependency state.

Conditional batch R2, mean-alignment RMSE, and residual logCPM variance
alignment are descriptive empirical diagnostics. The residual-variance score
is not negative-binomial dispersion, and the datasets do not provide uniform
biological ground truth.
