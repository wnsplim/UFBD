---
title: SBC
nav_order: 3
permalink: /sbc/
---

# SBC manual

`sbc` runs simulation-based calibration checking (Cook et al., 2006; Talts et al., 2020) of UFBD. For each replicate, it draws the parameters from the specified priors, simulates a dataset under the FBD process, analyzes the dataset with UFBD, and records the rank of each true value among the posterior samples. Building UFBD as in the [manual](manual.md#building) also builds `sbc`.

## Running

```
sbc -sbc_mode infer -sbc_reps <N> -sbc_out <prefix> [flags]
sbc -sbc_mode summarize -sbc_ranks <file,...> [-sbc_noderanks <file,...>] [flags]
sbc -help
```

`-sbc_mode simulate` only simulates the datasets, and `-sbc_mode emit` writes them as UFBD input files. `-sbc_ue <p>` makes each extant taxon an unsampled extant with probability `p`. The priors, the conditioning, the time bins and the MCMC settings take the same flags as `ufbd`.

`-sbc_mode summarize` prints the summary below for rank files, for example the per-run files of a batch merged into one; a node ranks file may carry a leading `chunk` column. Several comma-separated ranks files, with matching node ranks files, are summarized together, and each quantity is labeled by the last two directories of its file.

## Summary

For each quantity, `sbc` prints the ECDF statistic γ (Säilynoja et al., 2022; Modrák et al., 2025), its p-value against uniform ranks, and the Benjamini–Hochberg (`p_BH`; Benjamini & Hochberg, 1995) and Benjamini–Yekutieli (`p_BY`; Benjamini & Yekutieli, 2001) adjusted p-values over all quantities in the summary. A quantity fails at level q when its adjusted p-value is below q. γ_BH and γ_BY are the quantiles of γ under uniform ranks at the first rejection thresholds of the two procedures, printed with log(γ/γ_BH) and log(γ/γ_BY). The quantities are the columns of the ranks files, the origin (`x0`), the crown (`x1`), and one internal node per replicate chosen at random. `-sbc_q <q>` sets q (default 0.05), `-sbc_level <p>` adds γ_p, the `p` quantile, with log(γ/γ_p), and `-sbc_null <N>` sets the number of uniform rank sets.

## Output

* `<prefix>_ranks.tsv` — one row per replicate, with the rank of the true value of each parameter, `logLik` and `logPrior`
* `<prefix>_noderanks.tsv` — one row per node per replicate. `x0` is the origin and `x1` ... `xK` are the backbone nodes
* `<prefix>_rep<k>.log` and `<prefix>_rep<k>.trees` — the chain of replicate `k`

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful
approach to multiple testing. *Journal of the Royal Statistical Society Series B: Statistical Methodology,
57*(1), 289–300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>

Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate in multiple testing
under dependency. *The Annals of Statistics, 29*(4), 1165–1188. <https://doi.org/10.1214/aos/1013699998>

Cook, S. R., Gelman, A., & Rubin, D. B. (2006). Validation of software for Bayesian models using
posterior quantiles. *Journal of Computational and Graphical Statistics, 15*(3), 675–692.
<https://doi.org/10.1198/106186006X136976>

Modrák, M., Moon, A. H., Kim, S., Bürkner, P., Huurre, N., Faltejsková, K., Gelman, A., & Vehtari, A.
(2025). Simulation-based calibration checking for Bayesian computation: The choice of test
quantities shapes sensitivity. *Bayesian Analysis, 20*(2), 461–488. <https://doi.org/10.1214/23-BA1404>

Säilynoja, T., Bürkner, P.-C., & Vehtari, A. (2022). Graphical test for discrete uniformity and its
applications in goodness-of-fit evaluation and multiple sample comparison. *Statistics and
Computing, 32*(2), 32. <https://doi.org/10.1007/s11222-022-10090-6>

Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2020). Validating Bayesian
inference algorithms with simulation-based calibration. *arXiv*. <https://doi.org/10.48550/arXiv.1804.06788>
