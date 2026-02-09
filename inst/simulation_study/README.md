# Simulation studies (BA and CS)

This directory contains the main simulation studies for the **Barabási–Albert (BA)** and **change-statistic (CS)** mark models, plus SLURM job scripts to run them on a cluster.

## After `git pull`

1. **Install or load the package** from the package root (the directory that contains `R/`, `inst/`, `DESCRIPTION`):

   - From the terminal (one-time install):
     ```bash
     cd /path/to/hawkesNet
     R -e "devtools::install()"
     ```
   - Or in RStudio: open the project at the package root, then run `devtools::install()` or `devtools::load_all()`.

2. **Run a study** using one of the options below.

---

## Option A: Interactive RStudio (cloud or local)

1. Set the working directory to the **package root** (e.g. in RStudio: *Session → Set Working Directory → To Project Directory*, or `setwd("/path/to/hawkesNet")`).
2. Load the package:
   ```r
   devtools::load_all()   # development load, or
   library(hawkesNet)  # if already installed
   ```
3. Source the study script:
   ```r
   source("inst/simulation_study/simulation_study_BA.R")   # BA model
   # or
   source("inst/simulation_study/simulation_study_CS.R")  # CS model
   ```

Results and saved state will be written in the **current working directory** (e.g. `results_BA_full.RDS` / `results_CS_full.RDS`). For long runs (consistency study), consider using the SLURM option instead.

---

## Option B: SLURM (cluster batch job)

Submit from the **package root** so that `Rscript` can find the package and the script path is correct.

1. `cd` to the package root:
   ```bash
   cd /path/to/hawkesNet
   ```
2. Submit the job:
   ```bash
   sbatch inst/simulation_study/run_BA.slurm   # BA study
   sbatch inst/simulation_study/run_CS.slurm  # CS study
   ```

**SLURM scripts:**

| File | Model | Output / error logs |
|------|--------|----------------------|
| `run_BA.slurm` | BA (degree-weighted attachment) | `BA_sim_<jobid>.out`, `BA_sim_<jobid>.err` |
| `run_CS.slurm` | CS (ERGM-style change statistics) | `CS_sim_<jobid>.out`, `CS_sim_<jobid>.err` |

Defaults: 25 CPUs, 32 GB RAM, 50 h walltime. The R scripts read `SLURM_CPUS_PER_TASK` to set the number of cores; adjust `#SBATCH --cpus-per-task` and `--time` / `--mem` in the `.slurm` files if needed.

Outputs (e.g. `results_BA_full.RDS`, `results_CS_full.RDS`) are written in the **directory where the job was run** (package root).

---

## Scripts and outputs

- **simulation_study_BA.R** – BA mark model: main study, consistency (time windows), explosive regime; saves `results_BA_full.RDS`.
- **simulation_study_CS.R** – CS mark model: same structure with CS-specific options; saves `results_CS_full.RDS`.

Toggle at the top of each script: `SIMULATE`, `RUN_CONSISTENCY`, `RUN_EXPLOSIVE`, `PAPER_OUTPUT`. Full-state save and paper output (tables/plots) run at the end; see comments in the scripts for runtime estimates.
