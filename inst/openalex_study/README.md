# OpenAlex Hawkes study

Inhomogeneous HawkesNet fit on an OpenAlex citation/collaboration network, plus temporal Hawkes fit, KS tests, and ERGM-style goodness-of-fit (degree, ESP, geodesic distance, waiting times between triangle / 2-star / 3-star formations).

## Run from package root

From the **package root** (directory containing `R/`, `inst/`, `DESCRIPTION`):

```bash
Rscript inst/openalex_study/openalex_hawkes_study.R
```

## SLURM

Submit from the package root:

```bash
sbatch inst/openalex_study/run_openalex.slurm
```

Defaults: 8 CPUs, 32G mem, 24h. Edit `run_openalex.slurm` to change.

### One-shot: pull, install, submit

From package root or from `inst/`:

```bash
./inst/pull_install_submit.sh <BA|CS|openalex>
```

- **openalex** – OpenAlex Hawkes study (`run_openalex.slurm`)
- **BA** – Simulation study, Barabási–Albert (`inst/simulation_study/run_BA.slurm`)
- **CS** – Simulation study, Change Statistics (`inst/simulation_study/run_CS.slurm`)

The script runs: `git pull` → `module load R` → `R --no-save -e "devtools::install()"` → `sbatch <chosen slurm>`. Example: `./inst/pull_install_submit.sh openalex`

## Runtime (cloud / SLURM)

Yes. With the default settings the job should finish in **reasonable time** (typically well under the 24h limit). Rough breakdown:

- **OpenAlex fetch**: A few minutes (depends on `OPENALEX_PAGES` and API latency).
- **Inhomogeneous fit**: Dominant cost; scales with number of events (edges). The script uses multiple cores (from `SLURM_CPUS_PER_TASK`) for the intensity pre-computation, so 8 CPUs help. For hundreds of events expect tens of minutes; for a few thousand, possibly 1–3 hours.
- **Temporal fit + KS**: Usually a few minutes.
- **GOF**: 50 simulations from the fitted model; cost scales with event count. Can be on the order of 30 min to a few hours for large networks.

**If you want a quicker run**: set `OPENALEX_PAGES=10` or `20` (smaller network), or in the script set `N_GOF <- 20L` and/or `MAX_ITER <- 1000L`. For a full run with 50 pages and 50 GOF sims, 24h and 8 CPUs are generally sufficient.

## Requirements

The **hawkesGrowthNet** package must be installed (e.g. `devtools::install()` from the package root). It includes the inhomogeneous (KDE) fit and all helpers; no external scripts are needed. After cloning the repo, install the package and run the study from the package root.

## Environment

- **`OPENALEX_EMAIL`**  
  Email for OpenAlex API (required by OpenAlex).  
  Default: `duncan-clark@outlook.com`.

- **`OPENALEX_STRING`**  
  Search string for OpenAlex (e.g. author/topic).  
  Default: `Hawkes`.

- **`OPENALEX_PAGES`**  
  Number of pages to fetch (100 works per page).  
  Default: `50`.

Example:

```bash
export OPENALEX_EMAIL=your@email.com
export OPENALEX_STRING=Hawkes
export OPENALEX_PAGES=30
sbatch inst/openalex_study/run_openalex.slurm
```

## Output and rehydration

- Full state is saved to  
  `inst/openalex_study/results_openalex_full.RDS`  
  (network, inhomogeneous fit, temporal fit, KS p-value, GOF results, config).

- Rehydration: run the same script with `PAPER_OUTPUT <- TRUE` (default). It loads the RDS and prints fit summaries and GOF plots (degree, ESP, geodesic ECDF, waiting-time boxplots) without re-fitting or re-fetching data.

To only rehydrate (e.g. after a SLURM run), ensure the RDS path is correct and run:

```bash
Rscript inst/openalex_study/openalex_hawkes_study.R
```

from the package root; the script will load the RDS and produce the paper output.

## GOF metrics

- **Degree distribution**  
  Counts by degree (0 to `max_deg`); observed vs simulated (boxplots).

- **ESP (edge-wise shared partners)**  
  Distribution of shared partners per edge (0 to `k_max`); observed vs simulated.

- **Geodesic distance**  
  Distribution of pairwise distances (excluding unreachable); ECDF observed vs simulated.

- **Waiting times between formations**  
  For triangles, 2-stars (node degree ≥ 2), and 3-stars (node degree ≥ 3): the time *gaps* between consecutive formation events (each time the count of triangles / 2-stars / 3-stars increases). Observed vs simulated distributions (boxplots).

Simulated networks are generated from the fitted inhomogeneous HawkesNet using **`cond_intensity_inhom`** with time-varying background rate (same formula and mark PMF as the fit). The background rate `mu(t)` is computed from the KDE fit (`inhom_bg$mu_fit$mu_fun`) at each event time during simulation, ensuring the simulations match the fitted inhomogeneous model exactly.
