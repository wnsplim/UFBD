# UFBD
## Input

Three files:

- **Tree** — rooted, fully-resolved newick (-t)
- **Clades** (TSV) — one clade per row: `name <tab> taxon1,taxon2,…` (≥2 backbone taxa; the clade is their MRCA). (-c)
- **Fossils** (TSV) — one fossil per row: `taxon <tab> min_age <tab> max_age <tab> clade <tab> assignment`, where `assignment` is `CROWN` or `TOTAL`. (-f)