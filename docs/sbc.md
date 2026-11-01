---
title: SBC
nav_order: 3
permalink: /sbc/
---

# sbc manual

`sbc` runs simulation-based calibration of UFBD. For each replicate, it draws the parameters from the specified priors, simulates a dataset under the FBD process, analyzes the dataset with UFBD, and records the rank of each true value among the posterior samples. Building UFBD as in the [manual](manual.md#building) also builds `sbc`.

## Running

```
sbc -sbc_mode infer -sbc_reps <N> -sbc_out <prefix> [flags]
sbc -help
```

`-sbc_mode simulate` only simulates the datasets, and `-sbc_mode emit` writes them as UFBD input files. `-sbc_ue <p>` makes each extant taxon an unsampled extant with probability `p`. The priors, the conditioning, the time bins and the MCMC settings take the same flags as `ufbd`.

## Output

* `<prefix>_ranks.tsv` — one row per replicate, with the rank of the true value of each parameter, `logLik` and `logPrior`
* `<prefix>_noderanks.tsv` — one row per node per replicate. `x0` is the origin and `x1` ... `xK` are the backbone nodes
* `<prefix>_rep<k>.log` and `<prefix>_rep<k>.trees` — the chain of replicate `k`
