---
title: SBC
nav_order: 3
permalink: /sbc/
---

# SBC manual

`sbc` runs simulation-based calibration of UFBD. For each replicate, it draws the parameters from the specified priors, simulates a dataset under the FBD process, analyzes the dataset with UFBD, and records the rank of each true value among the posterior samples. Building UFBD as in the [manual](manual.md#building) also builds `sbc`.

## Running

```
sbc -sbc_mode infer -sbc_reps <N> -sbc_out <prefix> [flags]
sbc -sbc_mode summarize -sbc_ranks <file,...> [-sbc_noderanks <file,...>] [flags]
sbc -help
```

`-sbc_mode simulate` only simulates the datasets, and `-sbc_mode emit` writes them as UFBD input files. `-sbc_ue <p>` makes each extant taxon an unsampled extant with probability `p`. The priors, the conditioning, the time bins and the MCMC settings take the same flags as `ufbd`.

`-sbc_mode summarize` prints the summary below for rank files, for example the per-run files of a batch merged into one; a node ranks file may carry a leading `chunk` column. Several comma-separated ranks files, with matching node ranks files, are summarized together, and each quantity is labeled by the last two directories of its file.

## Summary

For each quantity, `sbc` prints the ECDF statistic γ of Säilynoja, Bürkner and Vehtari (2022), its p-value against uniform ranks, and the Benjamini–Hochberg (`p_BH`) and Benjamini–Yekutieli (`p_BY`) adjusted p-values over all quantities in the summary. A quantity fails at level q when its adjusted p-value is below q. γ_BH and γ_BY are the quantiles of γ under uniform ranks at the first rejection thresholds of the two procedures, printed with log(γ/γ_BH) and log(γ/γ_BY). The quantities are the columns of the ranks files, the origin (`x0`), the crown (`x1`), and one internal node per replicate chosen at random. `-sbc_q <q>` sets q (default 0.05), `-sbc_level <p>` adds γ_p, the `p` quantile, with log(γ/γ_p), and `-sbc_null <N>` sets the number of uniform rank sets.

## Output

* `<prefix>_ranks.tsv` — one row per replicate, with the rank of the true value of each parameter, `logLik` and `logPrior`
* `<prefix>_noderanks.tsv` — one row per node per replicate. `x0` is the origin and `x1` ... `xK` are the backbone nodes
* `<prefix>_rep<k>.log` and `<prefix>_rep<k>.trees` — the chain of replicate `k`
