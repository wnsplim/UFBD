# UFBD

UFBD performs fossilized birth-death (FBD) analysis under the unresolved formulation of the FBD process (Heath et al., 2014; Lim et al., 2026). With a fixed backbone tree, UFBD infers divergence times and the FBD rates. The FBD rates are the speciation rate, the extinction rate and the sampling rate. With only fossil records and no backbone tree, UFBD infers the FBD rates alone.

Please cite the following paper if you found this program to be useful:

* Lim, W., Raskin, L. Y., Li, J. K., Huelsenbeck, J., & Nielsen, R. (2026). Estimating divergence times and diversification rates with unresolved fossilized birth-death process. *TBA*

For comments, help or a bug report, file a GitHub issue or write to <david9456@berkeley.edu>.

## Building

The repository holds a static Linux x86-64 binary, `ufbd`, in `releases/`. On other platforms, build UFBD from the source. You need CMake 3.16 or later and a C++20 compiler. Eigen is included. UFBD links jemalloc only if jemalloc is present.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

`CMakeLists.txt` uses the flag `-march=x86-64-v3`. An ARM compiler rejects this flag, so remove it on Apple Silicon. The build uses GCC and Clang flags only. On Windows, use WSL, MinGW or Clang instead of MSVC.

## Running

```
ufbd -config config_file.cfg
ufbd -help
```

Set the options in a config file, or as flags on the command line. If you give both, the config file wins. UFBD reads every path relative to the directory you run the command from.

Two options are available only as a command-line flag:

`-resume` continues an interrupted run from its checkpoint. Run the same command again and add `-resume`.

`-no_latent_log` does not write `<prefix>_latent.log`. This saves disk space. The latent parameters in that file are usually of less interest, and they are rarely a mixing bottleneck.

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
```

Every family except `unif` and `empirical` accepts an offset. The offset moves the support of the distribution. For example, `exp:0.1,50.0` starts at 50.0, `truncnorm:55,2,50` is truncated at 50, and `improper:100` is flat above 100. `empirical` uses a file of samples as the prior. This lets you pass in the posterior of an earlier analysis with no parametric fit. Under `-age_offset`, UFBD shifts the samples with the rest of the model.

## Input files

`clade_def` — Tab-separated, with an optional header. A clade is every taxon tipward of the MRCA of the listed taxa. Two taxa are therefore enough to define a clade. Use this file only with a backbone tree. With no backbone tree, UFBD makes one clade *whole* that holds every taxon in the analysis.

```
clade	taxa
whole	taxa1,taxa2
genus1	taxa1,taxa3
```

`fossils` — Tab-separated, with an optional header. `assignment` sets where a fossil can attach, relative to its clade. A `CROWN` fossil attaches inside the crown group. A `TOTAL` fossil attaches anywhere in the total group. A `STEM` fossil attaches only on the stem stalk. A row with `max_age` 0 is an unsampled extant, and takes only `CROWN` or `TOTAL`. The optional sixth column gives the preservation type.

```
taxon	min_age	max_age	clade	assignment	type
taxa1	45.9000	48.8000	whole	STEM	lithic
taxa3	4.5000	4.9000	genus1	TOTAL	amber
```

`unsampled_extants` — Tab-separated, with an optional header. Each row gives a count of unsampled extants instead of a name for each one. One row is equal to that many `fossils` rows with `min_age` and `max_age` both 0. `assignment` takes only `CROWN` or `TOTAL`. UFBD names the taxa of a row `<clade>_<crown|total>_1` through `<clade>_<crown|total>_<number>`. These names must not already belong to a backbone tip or a `fossils` row. A clade can appear once as `CROWN` and once as `TOTAL`.

```
clade	assignment	number
genus1	TOTAL	400
genus1	CROWN	1000
genus2	CROWN	500
```

`backbone_tree` — NEWICK or NEXUS format. The tree must be rooted and bifurcating. UFBD holds the backbone topology fixed for the whole analysis. A backbone tip whose name matches a `fossils` row becomes a non-contemporaneous tip. UFBD then draws its age uniformly between the `min_age` and the `max_age` of that fossil. A fossil in the backbone tree cannot become a sampled ancestor.

`sequence` — FASTA, PHYLIP or NEXUS format. Supply this file for a full-likelihood analysis.

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

In `charpartition clock`, a `name:` starts a clock partition. Every charset after it belongs to that partition, until the next clock partition name. Commas and spaces both separate members, so `fast` above holds gene1 and gene2. You can also set the clock partitions in the config file or with a flag.

`hessian` — A Hessian file, the `in.BV` file of PAML. PAML (Yang, 2007), IQ-TREE (Demotte et al., 2025) and phyloHessian (Wang & Meade, 2026) all compute this file. Its topology and taxon labels must match `backbone_tree`. Supply this file for an approximate-likelihood analysis. The sequence partitions come from the Hessian file, and they carry no labels.

## The config file

Entries before the first `[section]` are global. An absent entry keeps its default. The last section lists every default.

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
* `<prefix>_chainN.log.ckp` — the checkpoint

At the end of the run, UFBD merges these logs into `<prefix>.log` and `<prefix>.trees`. The merge drops the burn-in, and renumbers the generation column as a running counter. UFBD also writes `<prefix>.tree`, the posterior mean tree.

With `parallel_chains = 1`, UFBD drops the `_chainN` suffix. No merge runs, so `<prefix>.log` keeps its burn-in. With `coupled_chains > 1`, UFBD logs the cold chain only.

`<prefix>.console.txt` holds the console output. If an unresolved fossil has more than one candidate attachment zone, `<prefix>_chainN_zones.tsv` gives the clade and the assignment behind each `az_` column of the latent log. With no backbone tree, UFBD does not write `<prefix>.trees`, `<prefix>.tree` or the zone files.

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

A chain shorter than `swap_interval` never attempts a swap.

### Convergence

This block applies only when `chain_length = auto`. The run stops when every parameter reaches `min_ess` bulk-ESS and `min_ess` tail-ESS, and its R-hat falls to `rhat` or lower. The run also stops at `max_gen`. UFBD leaves fixed parameters, the discrete variable `nSA`, and `posterior`, `likelihood` and `prior` out of the R-hat calculation.

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

`crown` puts the conditioning point at the root of the backbone tree. It needs two or more backbone tips to define the crown node of the clade. The other three options put the conditioning point at the origin of the clade. `origin` and `crown` condition on survival to the present, so they need one or more extant taxa. `crown` needs a backbone tree. `extinct` allows no extant taxa, and forces `rho` to 1. `anysample` is the loosest option, and needs one or more extant or fossil taxa. Under `crown`, no fossil can be `STEM` on the *whole* clade.

If you give no age prior after the keyword, the conditioning point gets a flat improper prior. This usually mixes poorly.

### `[lambda]`, `[mu]`, `[psi]`

Without `time_bins`, one rate applies at all times:

```ini
[lambda]
prior = exp:5
```

`time_bins` sets piecewise-constant rates. Each bin has a name, and `|` separates the bins. The bins must cover `[t_min, inf)` with no gaps and no overlaps, where `t_min` is the youngest bin edge. `+` joins intervals into one bin, so `recent_old` below uses one rate over [20,40) and [80,inf). `t_min` is 0 in most analyses. Use `t_min > 0` only under `extinct` conditioning, when you are sure the clade went extinct before `t_min`. All three FBD rates must share the same `t_min`.

```ini
[lambda]
time_bins = early:(0,20) | recent_old:(20,40)+(80,inf) | mid:(40,80)
prior     = early:exp:5 | recent_old:exp:5 | mid:unif:0,1
mode      = indep
```

An unnamed `prior` applies to every bin:

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
ou_theta  = 0.2,0.5      # theta ~ Normal(log(median), sd)
ou_sd     = 6,5          # gamma(shape,rate) on the OU log-SD
ou_nu     = 4,20         # gamma(shape,rate) on the reversion rate
```

`ou` replaces the per-bin priors, so `prior` and `mode = ou` together give an error. `ou` also rejects union bins, because piecewise smoothing is not compatible with them. With one bin only, UFBD uses `indep` and gives a warning, because there is nothing to smooth. The three hyperparameters are optional, and each one takes both of its numbers. The default `ou_nu` rate is `4 * dt / -log(0.7)`, where `dt` is the median distance between the midpoints of consecutive bins. This puts the prior mean reversion rate at `-log(0.7) / dt`, so the correlation between adjacent bins is 0.7.

### `[psi <type>]`

The sixth column of the `fossils` file gives an optional preservation type. Write one section for each type, with a matching label. Each type has its own `time_bins`, so the types can use different bins:

```ini
[psi amber]
time_bins = early:(0,30) | old:(30,inf)
prior     = early:exp:5 | old:exp:20
mode      = indep

[psi lithic]
prior = exp:5
```

Only the preservation rate splits into types. With one type, leave `<type>` out and write `[psi]` in the same way as `[lambda]` and `[mu]`. With two or more types, every fossil must carry a type that matches a section name.

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

The approximate-likelihood path does not use the sequence model. It ignores `partition`, `ctmc_model`, `ctmc_gamma_cat`, `ctmc_inv` and `ctmc_freq`, and takes the partitions from the Hessian file. This path needs the size of the state space. UFBD reads that size from `n_states`, or from `datatype` if `n_states` is absent. Amino-acid Hessian data therefore need `n_states = 20` or `datatype = aa`.

### `[clock]`

```ini
[clock]
clock_partitions = fast:(gene1,gene2) | slow:(gene3)
clock_model      = fast:gbm | slow:ucln
rgene_gamma      = 2,2000,1
sigma2_gamma     = 1,10,1
sigma2_param     = fast:nc | slow:c
```

A clock partition is a set of sequence partitions that share one clock model. With two or more clock partitions, the log columns for the clock parameters carry a `_<index>` suffix. The index starts at 0, in the order you declare the clock partitions. `rgene_gamma` and `sigma2_gamma` are `shape,rate,concentration` gamma–Dirichlet priors (dos Reis et al., 2014). They apply to the mean rates and to the clock-rate variances across clock partitions. With one clock partition, each prior reduces to a gamma.

`clock_partitions` takes `all`, `single`, or named groups that you write yourself. `all` is the default, and gives one clock to each sequence partition, named after that partition. `single` gives one clock to every sequence partition. In named groups, `|` separates the clock partitions and `,` separates the members of one partition. A member is a charset name or a sequence-partition index that starts at 0. Every sequence partition must appear in exactly one group. A `charpartition clock` block in the NEXUS file does the same, and supplies the labels. If you give both, `clock_partitions` wins and UFBD ignores the `charpartition clock` block. Under `hessian` there are no charset labels, so members must be indices, and UFBD names the `all` groups by index.

The index counts the sequence partitions in the order you declare them, and starts at 0. For the `partition` file above:

```
charset gene1 = 1-679;       -> 0
charset gene2 = 680-1377;    -> 1
charset gene3 = 1378-2044;   -> 2
```

Under `hessian`, the index follows the order of the blocks in the Hessian file:

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

There are four ways to group the sequence partitions into clock partitions:

```ini
clock_partitions = all                    # three clocks, named gene1, gene2, gene3
clock_partitions = single                 # one clock for all
clock_partitions = fast:(gene1,gene2) | slow:(gene3)
clock_partitions = fast:(0,1) | slow:(2)  # the same grouping by index, the only specifiable form under approximate-likelihood path
```

Use the clock-partition names above to set `clock_model` and `sigma2_param`. One value applies to every clock partition. You can also give one value for each clock partition.

```ini
clock_model = gbm                    # every clock partition gets GBM clock model
clock_model = fast:gbm | slow:ucln   # per clock partition clock models
```

`clock_model` is `ucln` or `gbm`, and the default is `ucln`. Clock partitions can use different models.

`sigma2_param` is `c` or `nc`, the centered and the non-centered parameterizations (Papaspiliopoulos et al., 2007). The default is `c`. For a short partition, up to about 10 kb, `nc` is much more efficient. The exact limit depends on how informative the sequence is. For long sequences and genome-scale Hessian files, `c` works better.

## Output log columns

UFBD builds the column labels in `<prefix>.log` from the model. It uses the bin labels and the preservation-type labels from the config file.

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

Bins that share a rate give one column under the shared label. For example, `x:(0,10)+(25,45) | y:(10,25) | z:(45,inf)` gives `lambda_x`, `lambda_y` and `lambda_z`. Under `mode = ou`, UFBD writes the per-bin rates and the three hyperparameters.

At the start of the run, UFBD maps every column label to its time interval under `Rate intervals:`. This map stays in `<prefix>.console.txt`.

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
