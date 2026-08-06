# UFBD

Fossilized birth-death (FBD) analysis under the unresolved formulation of the FBD process (Heath et al., 2014, Lim et al. 2026). UFBD can infer divergence times and FBD rates (speciation, extinction, and sampling rates) when given (fixed) backbone topology, or only FBD rates when only given fossil records without any backbone topology.

Please cite following paper if you found this program to be useful:

* Lim, W., Raskin, L. Y., Li, J. K., Huelsenbeck, J., & Nielsen, R. (2026). Estimating divergence times and diversification rates with unresolved fossilized birth-death process. *TBA*

For any comments, help or bug report, please file a GitHub issue, or contact: <david9456@berkeley.edu>

## Building

The repository ships a static Linux x86-64 binary, `ufbd` (see releases/). On other platforms, build from source with CMake 3.16 or newer and a C++20 compiler. Eigen is bundled and jemalloc is linked only if present.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

`CMakeLists.txt` compiles with `-march=x86-64-v3`, which an ARM compiler will reject, so the flag should be dropped on Apple Silicon. The build uses GCC/Clang flags throughout, so on Windows use WSL, MinGW or Clang rather than MSVC.

## Running

```
ufbd -config config_file.cfg
ufbd -help
```

Settings can be specified using a config file or on the command line as flags. If both are supplied, the config file is prioritized. Paths are resolved against the directory the command is run from.

Two user-facing options are available only as a command line flag:

`-resume` continues an interrupted run from the checkpoint. To resume, re-run the same command with `-resume` added.

`-no_latent_log` skips writing `<prefix>_latent.log`. Might save some storage as the latent parameters logged in these files are of less interest, and often are not a mixing bottleneck.

## Probability distribution grammar

```
exp:10            exponential:10          # rate
gamma:2,3                                 # shape, rate
lognormal:4.6,0.3                         # mu, sigma of the log
unif:0.03,0.06    uniform:0.03,0.06       # a, b
normal:55,2                               # mean, sd
truncnormal:55,2  truncnorm:55,2          # mean, sd, truncated at 0
improper                                  # improper Unif[0,inf) distribution
0.05                                      # a bare number for fixed value
```

Every family except `unif` can take an offset that shifts its support. For example, `exp:0.1,50.0` starts at 50.0, `truncnorm:55,2,50` is truncated below at 50, and `improper:100` is flat above 100.

## Input files

`clade_def` — Tab-separated, with an optional header. A clade is defined as every taxon under the MRCA of the listed taxa, so two taxa are enough to define the clade. This input is only valid when a backbone tree file is supplied. When no backbone tree is provided, a single clade *whole*, including all taxa in the analysis, is automatically generated even when not specified.

```
clade	taxa
whole	taxa1,taxa2
genus1	taxa1,taxa3
```

`fossils` — Tab-separated, with an optional header. `assignment` defines where the fossils can attach relative to their assigned clade. `CROWN` fossils can attach inside the crown group, `TOTAL` can attach anywhere on the total group, and `STEM` can only attach on the stem stalk. A row with `max_age` 0 is read as an unsampled extant, and can only be assigned to `CROWN` or `TOTAL`. The optional 6th column gives the preservation type (see below).

```
taxon	min_age	max_age	clade	assignment	type
taxa1	45.9000	48.8000	whole	STEM	lithic
taxa3	4.5000	4.9000	genus1	TOTAL	amber
```

`backbone_tree` — NEWICK or NEXUS formats are accepted, which must be rooted and bifurcating. Backbone tree topology is assumed to be fixed throughout the analysis. A backbone tip whose name matches a `fossils` row becomes a non-contemporaneous tip, with its age sampled uniformly within that fossil's `min_age`–`max_age` range. Note that when one include fossils in the backbone tree, the program will suppress that fossil from becoming a sampled ancestor.

`sequence` — FASTA, PHYLIP or NEXUS formats are accepted. Must be supplied for full-likelihood analysis.

`partition` — NEXUS `charset` blocks, optionally with a `charpartition clock` for manual clock partitioning. This file is named by the `partition` entry of `[substitution]` (see below).

```
#NEXUS
begin sets;
    charset gene1 = 1-679;
    charset gene2 = 680-1377;
    charset gene3 = 1378-2044;
    charpartition clock = fast: gene1 gene2, slow: gene3;
end;
```

In `charpartition clock`, a `name:` starts a clock partition definition and every charset after it belongs to that partition until the next clock partition name. Commas and spaces both separate members, so `fast` partition above includes gene1 and gene2. Clock partitions could be manually specified using config or flag as well (see below).

`hessian` — Hessian file (`in.BV` of PAML), which can be computed by programs like PAML (Yang, 2007), IQ-TREE (Demotte et al., 2025) or phyloHessian (Wang & Meade, 2026). Its topology and taxon labels must match `backbone_tree`. Must be supplied for approximate-likelihood analysis. The sequence partitions come from the Hessian file itself and carry no partition labels.

## The config file

Entries before any `[section]` are global, and an absent value leaves the setting at its default. Every default is listed at the end.

### Input and output

```ini
fossils       = fossils.tsv
clade_def     = clades.tsv
backbone_tree = backbone.tree
sequence      = align.fasta
# hessian     = in.BV             # for approximate-likelihood path, replacing sequence entry.
log_output    = output/prefix
tree_output   = output/prefix
```

`log_output` and `tree_output` give a prefix, to which the extensions below are appended.

Output files are written as follows. When `parallel_chains > 1` is invoked, each parallel chain writes `<prefix>_chainN.log`, `<prefix>_chainN.trees`, `<prefix>_chainN_latent.log` (attachment ages, attachment zones, fossil ages, and per-branch clock rates) and `<prefix>_chainN.log.ckp`, indexed from 0. After the run finishes, `<prefix>.log` and `<prefix>.trees` are written by merging the logs for each `parallel_chains` with burn-in dropped and the generation column renumbered as a running counter, also generating `<prefix>.tree`, the posterior mean tree. `<prefix>.console.txt` logs the terminal console output. When an unresolved fossil has more than one candidate attachment zone, `<prefix>_chainN_zones.tsv` gives the clade and assignment behind each assignment zone (`az_`) column of the latent log. With `parallel_chains = 1` the per-replicate suffix is dropped, and since no merging step runs, `<prefix>.log` keeps its burn-in. For `coupled_chains > 1`, only the cold chain is logged. `<prefix>.trees`, `<prefix>.tree` and the zone files are not written when no backbone tree is supplied.

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

With `coupled_chains > 1`, one can also optionally set MC3 parameters:

```ini
delta_temperature = 0.1      # initial chain spacing. Self-adapts
swap_interval     = 1000     # generations between swap attempts
```

A chain shorter than `swap_interval` never attempts a swap.

### Convergence

This block is valid only when `chain_length = auto`. The run stops once every parameter reaches `min_ess` bulk- and tail-ESS and falls to `rhat` or below, or reaches `max_gen`. Fixed parameters, discrete variables (`nSA`), and `posterior`, `likelihood` and `prior` are excluded from the calculation of R-hat.

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

`crown` places the conditioning point at the root of the backbone and needs at least two backbone tips that define the crown node of the clade. The other three place the conditioning at the origin of the clade. `origin` and `crown` condition on survival to the present, so they need at least one extant taxon. `crown` cannot be used when there is no backbone. `extinct` requires no extant samples, and forces `rho` to 1. `anysample` is the most loose conditioning, and requires at least one of either extant or fossil taxon. Under `crown` no fossil may be `STEM` on the *whole* clade.

Omitting he age prior after the keyword leaves the conditioning point with a flat improper prior (which generally leads to poor mixing).

### `[lambda]`, `[mu]`, `[psi]`

Without `time_bins` a rate is shared across all time points:

```ini
[lambda]
prior = exp:5
```

`time_bins` enables setting of piecewise-constant rates. Bins are named, separated by `|`, and must tile `[t_min, inf)` without gaps or overlaps, where `t_min` is the youngest bin edge. `+` unions intervals into one bin, so `recent_old` below assigns a single rate over [20,40) and [80,inf). For most of the analysis, `t_min = 0`, and `t_min != 0` only when one conditions on extinction, and is confident that the clade went extinct before `t_min`. `t_min` should be shared between all three FBD rates. 

```ini
[lambda]
time_bins = early:(0,20) | recent_old:(20,40)+(80,inf) | mid:(40,80)
prior     = early:exp:5 | recent_old:exp:5 | mid:unif:0,1
mode      = indep
```

An unnamed `prior` applies to every bin giving simpler representation:

```ini
[psi]
time_bins = b0:(0,20) | b1:(20,40) | b2:(40,inf)
prior     = exp:2
mode      = indep
```

`mode = ou` smooths the log-rates across bins through log-OU process:

```ini
[mu]
time_bins = early:(0,30) | mid:(30,60) | old:(60,inf)
mode      = ou
ou_theta  = 0.2,0.5      # theta ~ Normal(log(median), sd)
ou_sd     = 6,5          # gamma(shape,rate) on the OU log-SD
ou_nu     = 4,20         # gamma(shape,rate) on the reversion rate
```

`ou` replaces the per-bin priors, so setting `prior` and `mode = ou` together gives an error. `ou` also rejects union bins, as piecewise smoothing becomes logically incompatible under this formulation. With only one bin it falls back to `indep` with a warning as there is nothing to smooth. The three hyperparameters are optional and each takes both of its numbers. The default `ou_nu` rate is `4 * dt / -log(0.7)`, where `dt` is the median spacing between consecutive bin midpoints, which puts the prior mean reversion rate at `-log(0.7) / dt` so that the correlation between adjacent bins is 0.7.

### `[psi <type>]`

For each fossil preservation type (optionally) presented in the 6th column of the fossil table, one should write one section per type with matching labels. Each type carries its own `time_bins`, so the types need not share a bin structure:

```ini
[psi amber]
time_bins = early:(0,30) | old:(30,inf)
prior     = early:exp:5 | old:exp:20
mode      = indep

[psi lithic]
prior = exp:5
```

For now, only preservation rate could be split into multiple types. `<type>` argument can be omitted when there are no multiple preservation types, and `[psi]` could be set as same as `[lambda]` and `[mu]` following above. Once more than one type is declared, every fossils must carry a type matching a section name.

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

`ctmc_freq` chooses where the equilibrium frequencies come from. `estimated` samples them as free parameters. `empirical` fixes them at the frequencies observed in the alignment. `model` samples them under `gtr`, and under an empirical amino-acid matrix takes them from the matrix file and holds them fixed.

The sequence model is unused under `hessian`, so `partition`, `ctmc_model`, `ctmc_gamma_cat`, `ctmc_inv` and `ctmc_freq` are all ignored when approximate-likelihood path is used, and the partitions come from the Hessian file. Approximate-likelihood path requires size of the state space, and it is read from `n_states`, or from `datatype` when `n_states` is absent. Amino-acid Hessian data therefore need `n_states = 20` or `datatype = aa`.

### `[clock]`

```ini
[clock]
clock_partitions = fast:(gene1,gene2) | slow:(gene3)
clock_model      = fast:gbm | slow:ucln
rgene_gamma      = 2,2000,1
sigma2_gamma     = 1,10,1
sigma2_param     = fast:nc | slow:c
```

A clock partition is a set of sequence partitions sharing the clock model. With more than one clock partition the MCMC log columns of clock model parameters carry a `_<index>` suffix, indexed from 0 in the order the clock partitions were declared. `rgene_gamma` and `sigma2_gamma` are `shape,rate,concentration` gamma–Dirichlet priors (dos Reis et al., 2014) on the mean rates and the clock-rate variances across clock partitions, reducing to a gamma when there is only one.

`clock_partitions` takes `all` (default: one clock per sequence partition, named after it), `single` (one clock for every sequence partitions), or manually set named groups. When the clock partitions are manually set, `|` separates the clock partitions and `,` separates its members, where a member is a charset name or a 0-based sequence-partition index. Every sequence partition must appear in exactly one group. A `charpartition clock` in the NEXUS file does the same and supplies the clock partition labels. If both are given, `clock_partitions` is prioritized and the `charpartition clock` is ignored. Under `hessian` there are no charset labels, so members must be indices and the `all` groups are named by index.

The index counts the sequence partitions in the order they are declared, starting at 0. For the `partition` file shown earlier that is

```
charset gene1 = 1-679;       -> 0
charset gene2 = 680-1377;    -> 1
charset gene3 = 1378-2044;   -> 2
```

and under `hessian` it is the order the blocks appear in the Hessian file:

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

The four ways of grouping the sequence partition into clock partitions:

```ini
clock_partitions = all                    # three clocks, named gene1, gene2, gene3
clock_partitions = single                 # one clock for all
clock_partitions = fast:(gene1,gene2) | slow:(gene3)
clock_partitions = fast:(0,1) | slow:(2)  # the same grouping by index, the only specifiable form under approximate-likelihood path
```

`clock_model` and `sigma2_param` are set using the above defined clock-partition name. If supplied with a single value same setting applies to all clock partitions, or one per clock partition.

```ini
clock_model = gbm                    # every clock partition gets GBM clock model
clock_model = fast:gbm | slow:ucln   # per clock partition clock models
```

`clock_model` is `ucln` (the default) or `gbm`, and may differ between clock partitions.

`sigma2_param` is `c` (the default) or `nc`, the centered and non-centered parameterizations (see Papaspiliopoulos et al., 2007). For a partition with short total sequence length (up to ~10 kb; depends on the informativeness of the sequence) non-centered (`nc`) parameterization is much more efficient, and for long sequences and genome-scale Hessian files, centered (`c`) parameterizations tends to work better.

## Output log columns

Column labels in `<prefix>.log` are built from the model, using the bin and preservation-type labels given in the config.

```
n posterior likelihood prior
x0                                             origin age, with a backbone tree; skipped when -conditioning crown
x1 ... xK                                      backbone node ages
originAge                                      origin age, when ran without a backbone tree
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

Bins that share a rate are one column under the shared label, so `x:(0,10)+(25,45) | y:(10,25) | z:(45,inf)` gives `lambda_x`, `lambda_y` and `lambda_z`. Under `mode = ou` the per-bin rates are written alongside the three hyperparameters.

Every column label is mapped to its time interval at the start of the run, under `Rate intervals:`, and kept in `<prefix>.console.txt`.

`<prefix>_latent.log` holds the per-sample latent variables: `y_<taxon>` fossil attachment ages, `z_` attachment positions, `az_` attachment zones, `rate_` per-branch clock rates and `sa_` sampled-ancestor indicators.

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
ou_theta = 0.2,0.5
ou_sd    = 6,5
ou_nu    =                       # shape 4, rate 4 * dt / -log(0.7) where dt is median bin length
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
models. *EcoEvoRxiv*. https://doi.org/10.32942/X2CD2X

dos Reis, M., Zhu, T., & Yang, Z. (2014). The impact of the rate prior on Bayesian estimation of
divergence times with multiple loci. *Systematic Biology, 63*(4), 555–565.
https://doi.org/10.1093/sysbio/syu020

Heath, T. A., Huelsenbeck, J. P., & Stadler, T. (2014). The fossilized birth–death process for
coherent calibration of divergence-time estimates. *Proceedings of the National Academy of
Sciences, 111*(29), E2957–E2966. https://doi.org/10.1073/pnas.1319091111

Papaspiliopoulos, O., Roberts, G. O., & Sköld, M. (2007). A general framework for the
parametrization of hierarchical models. *Statistical Science, 22*(1), 59–73.
https://doi.org/10.1214/088342307000000014

Wang, S., & Meade, A. (2026). Molecular clock dating using complex mixture models: Applied to
ancient symbionts. *Molecular Biology and Evolution, 43*(3), msag039.
https://doi.org/10.1093/molbev/msag039

Yang, Z. (2007). PAML 4: Phylogenetic analysis by maximum likelihood. *Molecular Biology and
Evolution, 24*(8), 1586–1591. https://doi.org/10.1093/molbev/msm088
