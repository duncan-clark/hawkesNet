# Archived models and studies

The final package snapshot with the older growth, change-statistic, RHEM and
dyadic Hawkes interfaces is `hawkesNet` 0.2.0 at commit `1484437`. The snapshot
contains the original source, help, tests and checkout-based study scripts.
[Browse the archive](https://github.com/duncan-clark/hawkesNet/tree/1484437a3a491b875de3ad5d48ad4d9af2fbb393).

The pilots and development scripts retired in the earlier cleanup are in the
[preceding snapshot](https://github.com/duncan-clark/hawkesNet/tree/476f48df57e0389dad8a854c355a7c06241ec80b).

| Workflow | Location |
|---|---|
| Older preferential-attachment and CS simulation studies | `inst/simulation_study/` in the snapshot |
| Hypertext, OpenAlex and college messages | `inst/hypertext_conference/`, `inst/openalex_study/` and `inst/message_network/` in the snapshot |
| RHEM, valued repeated hits and joint dyadic Hawkes models | `inst/RHEM_TIMENET.md`, `R/rhem.R`, `R/mark_PMF.R` and `R/joint_marked_hawkes.R` in the snapshot |
| Preliminary CRAN fit, forecast and report scripts | `inst/cran/` in the snapshot |
| AMLSim and Ethereum applications | `studies/amlsim/` and `studies/ethereum/` in the research workspace |

## Isolated installation

From the current package checkout, create a detached worktree and a separate
R library. The archived dependencies below pin the valued `ernm` backend and
the timestamped `timeNet` backend. The setup requires R, a C++ compiler and Git.

```sh
git worktree add --detach ../hawkesNet-legacy 1484437a3a491b875de3ad5d48ad4d9af2fbb393
mkdir -p ../hawkesNet-legacy-library
export HAWKESNET_LEGACY_LIBRARY="$(cd ../hawkesNet-legacy-library && pwd)"
export R_LIBS="$HAWKESNET_LEGACY_LIBRARY${R_LIBS:+:$R_LIBS}"

Rscript --vanilla - <<'RSCRIPT'
lib <- Sys.getenv("HAWKESNET_LEGACY_LIBRARY")
cran <- "https://cloud.r-project.org"
packages <- c("remotes", "data.table", "hash", "sna", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = cran, lib = lib)
options(repos = c(CRAN = cran))
remotes::install_github("duncan-clark/ernm@cb7256d70b38c5c345d174b39e1a22c85a85be00",
                        lib = lib, dependencies = NA, upgrade = "never")
remotes::install_github("duncan-clark/timeNet@40fc945d61027e745b4eb61da582e9e1f4980032",
                        lib = lib, dependencies = NA, upgrade = "never")
RSCRIPT

R CMD INSTALL -l "$HAWKESNET_LEGACY_LIBRARY" ../hawkesNet-legacy
Rscript --vanilla -e 'stopifnot(packageVersion("hawkesNet") == "0.2.0"); print(find.package("hawkesNet"))'
```

Keep the legacy `R_LIBS` setting in the study shell and use a separate shell
for the current package. For archived fits, use the study's original data and
dependency versions where recorded.

## Running the studies

For studies inside the snapshot, follow their README files from the archived
worktree. Set `HAWKESNET_OUTPUT_DIR` to a separate results directory. The
archived NeSI helper `pull_install_submit.sh` pulls a moving branch and is
unsuitable for the detached worktree. Submit the relevant `.slurm` script
directly after configuring the R library, output and scheduler log paths.

Some historical CRAN diagnostics expect the package source at `hawkesNet/`
relative to the working directory. Create a separate project directory and
name the pinned worktree `hawkesNet` inside it. Copy the required diagnostics
and inputs into their expected relative paths before running the study. This
keeps its table and report outputs separate from the active manuscript.

For the sibling applications, open a shell at the research workspace root
with the archive library still first in `R_LIBS`. Follow the corresponding
`studies/` README for input paths and run limits. The source and data for these
applications live in that workspace. The study runners use the installed
archive in this shell.

```sh
AMLSIM_TX_PATH=/path/to/transactions.csv Rscript --vanilla studies/amlsim/run_amlsim_study.R
ETH_TRANSFERS_CSV=/path/to/transfers.csv Rscript --vanilla studies/ethereum/run_ethereum_anomaly_study.R
```
