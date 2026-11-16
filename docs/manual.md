---
title: Manual
nav_order: 2
permalink: /manual/
---

# UFBD manual

## Building

The repository holds a static Linux x86-64 binary, `ufbd`, in `releases/`. On other platforms, you can build UFBD from the source, using e.g. CMake 3.16 or later and a C++20 compiler. UFBD links jemalloc only if jemalloc is present.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

`CMakeLists.txt` uses the flag `-march=x86-64-v3`. An ARM compiler rejects this flag, so Apple Silicon you'd likely need it to be removed. The build uses GCC and Clang flags only. On Windows, use WSL, MinGW or Clang instead of MSVC.

## Running

```
ufbd -config config_file.cfg
ufbd -help
```

You can set the options in a config file, or as flags on the command line. If you give both, the config file is prioritized. UFBD reads every path relative to the directory you run the command from.

Two options are available only as a command-line flag:

`-resume` continues an interrupted run from its checkpoint. Run the same command again and add `-resume`.

`-no_latent_log` suppresses writing `<prefix>_latent.log` file. This saves (tiny bit of) disk space. The latent parameters in that file are usually of less interest, and they are rarely a mixing bottleneck.

## Probability distribution grammar

```
exp:10            exponential:10          # rate
gamma:2,3                                 # shape, rate
lognormal:4.6,0.3                         # mu, sigma of the log
unif:0.03,0.06    uniform:0.03,0.06       # a, b
normal:55,2                               # mean, sd
truncnormal:55,2  truncnorm:55,2          # mean, sd, truncated at 0
improper                                  # improper Unif[0,inf) distribution
empirical:samples.txt                     # file of samples, one per line
0.05                                      # a bare number for fixed value
normal:log(0.2),1.6                       # any number may be written log(x) (treated as natural logarithm)
```

Every distribution family except `unif` and `empirical` accepts an offset. The offset moves the support of the distribution. For example, `exp:0.1,50.0` starts at 50.0, `truncnorm:55,2,50` is truncated at 50, and `improper:100` is flat above 100. `empirical` uses a file of samples as the prior. This lets you pass in the posterior of an earlier analysis with no parametric fit. Under `-age_offset`, UFBD shifts the samples with the rest of the model.

## Input files

`clade_def` — Tab-separated, with an optional header. A clade is defined as every taxon tipward of the MRCA of the listed taxa. Thus, two taxa are enough to define a clade. Use this file only with a backbone tree. With no backbone tree, UFBD makes one clade *whole* that holds every taxon in the analysis.

```
clade	taxa
whole	taxa1,taxa2
genus1	taxa1,taxa3
```

`fossils` — Tab-separated, with an optional header. `assignment` sets where a fossil can attach, relative to its `clade`. A `CROWN` fossil attaches inside the crown group. A `TOTAL` fossil attaches anywhere in the total group. A `STEM` fossil attaches only on the stem stalk. A row with `max_age` 0 is treated as an unsampled extant, and takes only `CROWN` or `TOTAL`. The optional sixth column gives the fossil preservation type.

```
taxon	min_age	max_age	clade	assignment	type
taxa1	45.9000	48.8000	whole	STEM	lithic
taxa3	4.5000	4.9000	genus1	TOTAL	amber
```

`unsampled_extants` — Tab-separated, with an optional header. Each row gives a count of unsampled extants instead of a name for each one. UFBD names the taxa of a row `<clade>_<crown|total>_1` through `<clade>_<crown|total>_<number>`. A clade can appear once as `CROWN` and once as `TOTAL`. This is a good replacement input format when you have many (say, thousands of) unsampled extants.

```
clade	assignment	number
genus1	TOTAL	400
genus1	CROWN	1000
genus2	CROWN	500
```

`backbone_tree` — NEWICK or NEXUS format. The tree must be rooted and bifurcating. UFBD holds the backbone topology fixed for the whole analysis. A backbone tip whose name matches a `fossils` file row becomes a non-contemporaneous tip. UFBD then draws its age uniformly between the `min_age` and the `max_age` of that fossil. A fossil in the backbone tree cannot become a sampled ancestor.

If you want to fix where a fossil sits in the backbone tree (say, if fossil-only clade's position is very well known with regards to extant clades), put it in `backbone_tree` file as a tip and add it in a row in `fossils` file under the same name. Its `clade` must be defined in `clade_def` file, and its `assignment` follows the same rules as for an unresolved fossil. A clade in `clade_def` file can be defined only with fossil tips. Under `extinct` conditioning, the backbone tree (if given) can hold fossil tips only, with no sequence data.

`sequence` — FASTA, PHYLIP or NEXUS format.

`partition` — NEXUS `charset` blocks. Add a `charpartition clock` block to set the clock partitions by hand. The `partition` entry of `[substitution]` names this file.

```
#NEXUS
begin sets;
    charset gene1 = 1-679;
    charset gene2 = 680-1377;
    charset gene3 = 1378-2044;
    charpartition clock = fast: gene1 gene2, slow: gene3;
end;
```

In `charpartition clock`, a `name:` starts a clock partition. You can also set the clock partitions in the config file or with a flag (see below).

`hessian` — A Hessian file, the `in.BV` file of PAML. PAML (Yang, 2007), IQ-TREE (Demotte et al., 2025) and phyloHessian (Wang & Meade, 2026) all compute this file format. Its topology and taxon labels must match `backbone_tree`. Supply this file for an approximate-likelihood analysis. When using aprpoximate likelihood path, the sequence partition definition come from the Hessian file. Note in this case partitions do not carry labels.

## The config file

An absent entry keeps its default.

### Input and output

```ini
fossils       = fossils.tsv
clade_def     = clades.tsv
# unsampled_extants = ue.tsv      # counts of unsampled extants per clade/assignment
backbone_tree = backbone.tree
sequence      = align.fasta
# hessian     = in.BV             # for approximate-likelihood path, replacing sequence entry.
log_output    = output/prefix
tree_output   = output/prefix
```

`log_output` and `tree_output` give a prefix. UFBD adds the extensions below to that prefix.

With `parallel_chains > 1`, each chain writes four files, indexed from 0:

* `<prefix>_chainN.log`
* `<prefix>_chainN.trees`
* `<prefix>_chainN_latent.log` — attachment ages, attachment zones, fossil ages and per-branch clock rates
`<prefix>_chainN_zones.tsv` — clade and the assignment definitions
* `<prefix>_chainN.log.ckp` — the checkpoint file

At the end of the run, UFBD merges these logs into `<prefix>.log` and `<prefix>.trees` with burn-in dropped. UFBD also writes the posterior mean tree `<prefix>.tree`. With `coupled_chains > 1`, UFBD logs the cold chain only. When no backbone tree is given, UFBD does not write `<prefix>.trees`, `<prefix>.tree`, and `<prefix>_chainN_zones.tsv`.

### MCMC

```ini
chain_length    = auto       # or a generation count
parallel_chains = 4          # number of parallel chains
coupled_chains  = 4          # number of Metropolis-coupled chains
thinning        = 1000
burn_in         = 0.25
cores           = 16
seed            = 42
```

With `coupled_chains > 1`, you can also set the MC3 parameters:

```ini
delta_temperature = 0.1      # initial chain spacing. Self-adapts
swap_interval     = 1000     # generations between swap attempts
```

### Convergence

This block applies only when `chain_length = auto`. The run stops when every parameter reaches `min_ess` bulk-ESS and `min_ess` tail-ESS, and its R-hat falls to `rhat` or lower. The run also stops at `max_gen`.

```ini
max_gen = 1000000000
min_ess =                    # default 100 x parallel_chains
rhat    = 1.01
```

### FBD model

```ini
conditioning = origin exp:0.1   # origin | crown | anysample | extinct, then an age prior for origin (x0) 
                                # or crown (x1; only for crown conditioning)
rho          = 1                # fraction of extant species sampled
```

`origin` and `crown` condition on survival to the present (either at the origin or crown node of the clade), so they need one or more extant taxa. `crown` always needs a backbone tree to be defined. `extinct` requires no extant taxa, and forces `rho` to 1. With `extinct`, the backbone tree (if given) can hold fossil tips only, each named in a `fossils` row. `anysample` is the loosest option, and needs one or more extant or fossil taxa.

If you give no age prior after the keyword, the conditioning point gets a flat improper prior. This usually mixes poorly.

### `[lambda]`, `[mu]`, `[psi]`

Without `time_bins`, constant rate is assumed for that parameter:

```ini
[lambda]
prior = exp:5
```

`time_bins` sets piecewise-constant rates. Each bin has a name, and `|` separates the bins. The bins collectively must cover `[t_min, inf)` with no gaps and no overlaps, where `t_min` is the youngest bin edge. `+` joins intervals into one bin, so `recent_old` below uses one rate over [20,40) and [80,inf). `t_min` is 0 in most analyses. Use `t_min > 0` only under `extinct` conditioning, when you are sure the clade went extinct before `t_min`. All three FBD rates must share the same `t_min`.

```ini
[lambda]
time_bins = early:(0,20) | recent_old:(20,40)+(80,inf) | mid:(40,80)
prior     = early:exp:5 | recent_old:exp:5 | mid:unif:0,1
mode      = indep
```

Below sets same `prior` to every bin:

```ini
[psi]
time_bins = b0:(0,20) | b1:(20,40) | b2:(40,inf)
prior     = exp:2
mode      = indep
```

`mode = ou` smooths the log-rates across bins with a log-OU process:

```ini
[mu]
time_bins = early:(0,30) | mid:(30,60) | old:(60,inf)
mode      = ou
ou_theta  = normal:log(0.2),1.6  # prior on the OU level theta (a log rate)
ou_sd     = exp:3.3333   # prior on the OU log-SD
ou_nu     = exp:2        # prior on the reversion rate per Myr
```

`ou` replaces the per-bin priors, so `prior` and `mode = ou` together give an error. `ou` also errors when joined bins (e.g., `recent_old:(20,40)+(80,inf)` above) are set, as this makes piecewise smoothing incompatible.

`ou_theta` sets the magnitude of the long-term rate. `normal:log(m),s` puts e^θ within m·e^(±2s) with 95% probability. The default `normal:log(0.2),1.6` give you about 0.009 to 4.6 events per lineage per Myr.

`ou_sd` sets how far a bin's rate can sit with regards to the long-term rate e^θ. About 95% of bins lie within a factor e^(2σ) of it: σ = 0.35 allows a factor of 2, σ = 0.55 a factor of 3. The default `exp:3.3333` gives σ a mean of 0.3.

`ou_nu` sets how strongly bins are tied. The log rates of two bins whose midpoints are Δ Myr apart have correlation e^(−νΔ). For a median correlation c between neighboring bins of width Δ, use `ou_nu = exp:r` with r = 0.693·Δ/(−ln c). For c = 0.7 this is r ≈ 1.94·Δ, close to the default `exp:2`. The data rarely inform ν, so this prior largely sets the smoothing.

### `[psi <type>]`

The sixth optional column of the `fossils` file gives an optional preservation type. When multiple preservation types are set, write one section for each type, with a matching label. Each type has its own `time_bins`, so the types can use different bins:

```ini
[psi amber]
time_bins = early:(0,30) | old:(30,inf)
prior     = early:exp:5 | old:exp:20
mode      = indep

[psi lithic]
prior = exp:5
```

Only the preservation rate splits into types. With one type, leave `<type>` out and write `[psi]` in the same way as `[lambda]` and `[mu]`.

### `[substitution]`

```ini
[substitution]
partition      = parts.nex   # NEXUS charsets; blank = single partition
datatype       = nt          # nt | aa
n_states       = 4           # for hessian (approximate-likelihood) path only
ctmc_model     = gtr         # gtr, or a PAML matrix file for amino acids
ctmc_gamma_cat = 4
ctmc_inv       = off         # on | off
ctmc_freq      = model       # model | empirical | estimated
```

`ctmc_freq` sets the source of the equilibrium frequencies. `estimated` samples them as free parameters. `empirical` holds them at the frequencies in the alignment. Under `gtr`, `model` samples them. Under an empirical amino-acid matrix, `model` takes them from the matrix file and holds them fixed.

The approximate-likelihood path does not use the sequence model. When this path is taken, UFBD ignores `partition`, `ctmc_model`, `ctmc_gamma_cat`, `ctmc_inv` and `ctmc_freq`, and takes the partitions from the Hessian file. Under approximate-likelihood path the size of the state space of the sequence data must be given. UFBD reads that size from `n_states`, or from `datatype` if `n_states` is absent. For example, Hessian file generated from the amino-acid sequneces need `n_states = 20` or `datatype = aa`.

### `[clock]`

```ini
[clock]
clock_partitions = fast:(gene1,gene2) | slow:(gene3)
clock_model      = fast:gbm | slow:ucln
rgene_gamma      = 2,2000,1
sigma2_gamma     = 1,10,1
sigma2_param     = fast:nc | slow:c
```

A clock partition is a set of sequence partitions that share one clock model. With two or more clock partitions, the columns for the clock parameters in the log file carry a `_<index>` suffix. The index starts at 0, in the order you declare the clock partitions. `rgene_gamma` and `sigma2_gamma` follow `shape,rate,concentration` gamma–Dirichlet priors (dos Reis et al., 2014). With single clock partition, each prior reduces to a gamma.

`clock_partitions` takes `all`, `single`, or named groups that you write yourself. `all` is the default, and gives one clock to each sequence partition, named after that partition. `single` gives one clock to every sequence partition. In named groups, `|` separates the clock partitions and `,` separates the members of one partition. A member is a charset name or a sequence-partition index that starts at 0. Every sequence partition must appear in exactly one group. A `charpartition clock` block in the NEXUS file does the same thing. If you give both, `clock_partitions` is prioritized and UFBD ignores the `charpartition clock` block in the NEXUS file.

If no explicit partition labels are given, the index counts the sequence partitions in the order you declare them, and starts at 0. For example:

```
charset gene1 = 1-679;       -> 0
charset gene2 = 680-1377;    -> 1
charset gene3 = 1378-2044;   -> 2
```

Under the approximate-likelihood path, the index follows the order of the blocks in the Hessian file:

```
 72                          -> 0
((Slon24923_Stenella_longirostris: 0.000185, ... );
  0.000287  0.000546 ...
  0.185047  3.047834 ...
Hessian
 -3.354e+09   -2.7e+07 ...

 72                          -> 1
...
```

Summary of the four ways to group the sequence partitions into clock partitions:

```ini
clock_partitions = all                    # three clocks, named gene1, gene2, gene3
clock_partitions = single                 # one clock for all
clock_partitions = fast:(gene1,gene2) | slow:(gene3)
clock_partitions = fast:(0,1) | slow:(2)  # the same grouping by index, the only specifiable form under approximate-likelihood path
```

You can use the clock-partition names above to set `clock_model` and `sigma2_param`.

```ini
clock_model = gbm                    # every clock partition gets GBM clock model
clock_model = fast:gbm | slow:ucln   # per clock partition clock models
```

`clock_model` is `ucln` or `gbm`, and the default is `ucln`.

`sigma2_param` is `c` or `nc`, the centered and the non-centered parameterizations (Papaspiliopoulos et al., 2007). The default is `c`. For a short partition, up to about 10 kb, `nc` is much more efficient. The exact limit depends on how informative the sequence is. For long sequences and genome-scale Hessian files, `c` mixes better.

## Output log columns

```
n posterior likelihood prior
x0                                             origin age, with sequence data or a Hessian file; skipped when -conditioning crown
x1 ... xK                                      backbone node ages
originAge                                      origin age, without sequence data or a Hessian file; skipped when -conditioning crown
nSA                                            number of sampled ancestors
lambda mu psi                                  a rate with when no time_bins are set
lambda_early lambda_mid ...                    a rate with time_bins set, one per bin label
ou_lambda_theta ou_lambda_sdEq ou_lambda_nu    added by mode = ou
psi_amber_early psi_amber_old psi_lithic       psi split by preservation type
ou_psi_amber_theta ...                         mode = ou on one preservation type
clockMean clockSigma2                          when only one clock partition is set
clockMean_0 clockSigma2_0 ...                  when multiple clock partitions are set
exch<k>_* freq<k>_* alpha<k> pinv<k>           sequence model of sequence partition k
```

Bins that share a rate give one column under the shared label. For example, `x:(0,10)+(25,45) | y:(10,25) | z:(45,inf)` gives `lambda_x`, `lambda_y` and `lambda_z`.

`<prefix>_latent.log` holds the latent variables of each sample. `y_<taxon>` is a fossil attachment age, `z_` an attachment position, `az_` an attachment zone, `rate_` a per-branch clock rate and `sa_` a sampled-ancestor indicator.

## Defaults

```ini
chain_length      = 1000000
parallel_chains   = 4
coupled_chains    = 1
thinning          = 1000
burn_in           = 0.25
cores             = 1
seed              =              # drawn from the system RNG
max_gen           = 1000000000
min_ess           =              # 100 x parallel_chains
rhat              = 1.01
rho               = 1
delta_temperature = 0.1
swap_interval     = 1000
```

The same defaults hold in `[lambda]`, `[mu]` and every `[psi <type>]`:

```ini
prior    = exp:5
mode     = indep
ou_theta = normal:log(0.2),1.6
ou_sd    = exp:3.3333
ou_nu    = exp:2
```

```ini
[substitution]
partition      =                 # one sequence partition
datatype       = nt
n_states       = 4
ctmc_model     = gtr
ctmc_gamma_cat = 4
ctmc_inv       = off
ctmc_freq      = model
```

```ini
[clock]
clock_partitions =               # the charpartition clock if present, otherwise all
clock_model      = ucln
sigma2_param     = c
rgene_gamma      = 2,2000,1
sigma2_gamma     = 1,10,1
```

## References

Demotte, P., Panchaksaram, M., Kumarasinghe, H., Ly-Trong, N., dos Reis, M., & Minh, B. Q. (2025).
IQ2MC: A new framework to infer phylogenetic time trees using IQ-TREE 3 and MCMCTree with mixture
models. *EcoEvoRxiv*. <https://doi.org/10.32942/X2CD2X>

dos Reis, M., Zhu, T., & Yang, Z. (2014). The impact of the rate prior on Bayesian estimation of
divergence times with multiple loci. *Systematic Biology, 63*(4), 555–565.
<https://doi.org/10.1093/sysbio/syu020>

Heath, T. A., Huelsenbeck, J. P., & Stadler, T. (2014). The fossilized birth–death process for
coherent calibration of divergence-time estimates. *Proceedings of the National Academy of
Sciences, 111*(29), E2957–E2966. <https://doi.org/10.1073/pnas.1319091111>

Papaspiliopoulos, O., Roberts, G. O., & Sköld, M. (2007). A general framework for the
parametrization of hierarchical models. *Statistical Science, 22*(1), 59–73.
<https://doi.org/10.1214/088342307000000014>

Wang, S., & Meade, A. (2026). Molecular clock dating using complex mixture models: Applied to
ancient symbionts. *Molecular Biology and Evolution, 43*(3), msag039.
<https://doi.org/10.1093/molbev/msag039>

Yang, Z. (2007). PAML 4: Phylogenetic analysis by maximum likelihood. *Molecular Biology and
Evolution, 24*(8), 1586–1591. <https://doi.org/10.1093/molbev/msm088>
