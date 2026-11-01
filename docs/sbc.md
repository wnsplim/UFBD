---
title: sbc
nav_order: 3
permalink: /sbc/
---

# sbc manual

`sbc` runs simulation-based calibration of UFBD. It draws the parameters from specified priors, simulates a dataset under the FBD process, analyzes the dataset with one UFBD chain, and records the rank of true value falls among the posterior samples. The build in the [manual](manual.md#building) also produces `sbc`.

## Running

```
sbc -sbc_mode infer -sbc_reps <N> -sbc_out <prefix> [flags]
sbc -help
```

`-sbc_mode infer` simulates and analyzes `-sbc_reps` datasets. `-sbc_mode simulate` prints the sizes of the simulated datasets, and `-sbc_mode emit` writes the simulated datasets as UFBD input files.

`-sbc_ue` gives the probability that each extant taxon is an unsampled extant, and the default is `0`. The priors, the conditioning, the time bins, `mode = ou` and the MCMC settings use the `ufbd` flags.

## Output

`-sbc_mode infer` writes:

* `<prefix>_ranks.tsv` — one row per replicate: the randomized normalized rank of the true value of each rate, each OU hyperparameter, `logLik` and `logPrior`, then `bulkESS`, `tailESS`, `maxRhat` and the true values `t_lambda`, `t_mu`, `t_psi`, `t_logLik` and `t_logPrior`. `t_lambda`, `t_mu` and `t_psi` are the rates of the youngest bin
* `<prefix>_noderanks.tsv` — one row per node per replicate, with the columns `rep`, `node` and `rank`. `x0` is the origin, `originAge` in the replicate log, and `x1` ... `xK` are the backbone nodes
* `<prefix>_rep<k>.log` and `<prefix>_rep<k>.trees` — the chain of replicate `k`

## Settings of Lim et al. (2026)

Lim et al. (2026) run each setting below as 100 runs of 50 replicates. Run `k` uses `-seed <base>+k`, and `-fbd_model UFBD` and `-fbd_model HEA14` use the same seeds. Every run uses

```
-sbc_mode infer -sbc_reps 50 -chain_length 5000000 -thinning 1000 -burn_in 0.3 -parallel_chains 1 -cores 1 -rho 1
```

and adds:

```
base 110000   -lambda_prior unif:0.03,0.06 -mu_prior unif:0.005,0.025 -psi_prior unif:0.08,0.12
              -conditioning origin lognormal:4.25,0.12 -sbc_ue 0
base 111000   base 110000 with -sbc_ue 0.5
base 112000   base 110000 with -sbc_ue 1
base 113000   -lambda_prior unif:0.06,0.08 -mu_prior unif:0.04,0.06 -psi_prior unif:0.3,0.5
              -conditioning extinct lognormal:4.25,0.12 -sbc_ue 1
base 114000   base 111000 with -lambda_skyline_times 20,40 -mu_skyline_times 20,40 -psi_skyline_times 20,40
base 115000   base 114000 with -lambda_prior_mode ou -mu_prior_mode ou -lambda_ou_theta 0.045,0.15
              -mu_ou_theta 0.015,0.15 -lambda_ou_sd 20,100 -mu_ou_sd 20,100 in place of -lambda_prior and -mu_prior
```
