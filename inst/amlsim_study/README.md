# AMLSim Valued HawkesNet Study

This directory contains an exact-valued HawkesNet study scaffold for IBM
AMLSim data.

Primary question:

- Can an exact valued repeated-hit HawkesNet model identify suspicious/SAR
  transactions competitively with standard ML baselines?

Secondary question:

- How far does exact full-dyad evaluation scale as the number of actors and
  events increases?

## Data

The importer supports common AMLSim transaction formats:

- converted output: `transactions.csv` with fields like `orig_acct`,
  `bene_acct`, `base_amt`, `tran_timestamp`, `is_sar`, `alert_id`
- sample/simple output: `tx.csv` with fields like `ACCOUNT_ID`,
  `COUNTER_PARTY_ACCOUNT_NUM`, `TXN_AMOUNT_ORIG`, `start`
- graph-generator output: files with fields like `id`, `src`, `dst`, `ttype`
  plus optional alert membership files such as `alert_members.csv`

The standardized repeated-hit format is:

- `t`: strict event time in `[0, 1]`
- `i`, `j`: integer sender/receiver actor IDs
- `weight`: transaction value, defaulting to rescaled `log1p(amount)`
- `is_sar`: SAR/suspicious transaction label

AMLSim often has multiple transactions in one simulation step. The study code
preserves gaps between simulator timestamps and only spreads same-step ties
within the smallest observed time gap because the current HawkesNet likelihood
expects one observed mark per event time.

## Running

From the package root:

```r
pkgload::load_all()
source("inst/amlsim_study/run_amlsim_study.R")
```

For the timestamped study, generate Java simulator output with AMLSim and use
`external/AMLSim/outputs/sample`. You can override paths with environment
variables:

```bash
AMLSIM_TX_PATH=/path/to/transactions.csv \
AMLSIM_MAX_ACTORS=100 \
AMLSIM_MAX_EVENTS=1000 \
AMLSIM_HAWKES_MAXIT=50 \
Rscript inst/amlsim_study/run_amlsim_study.R
```

Outputs are written to `amlsim_study_results/` unless
`AMLSIM_STUDY_OUTPUT_DIR` is set.

## Models

The HawkesNet SAR score is an anomaly score:

1. train exact valued RHEM-HawkesNet on the non-SAR training prefix
2. score test transactions using negative log fitted conditional intensity
3. evaluate whether high anomaly scores recover SAR labels

The current study optimizes the Hawkes ground-process parameters `mu`, `K`, and
`beta_overall` along with the RHEM coefficients. It keeps `beta_edges` fixed at
zero in the timestamped rerun so the first pass isolates the ground-process
timing contribution from the valued repeated-hit mark distribution.

The exact valued RHEM formula is:

```r
edgeValue + recipValue + senderValueActivity + receiverValueActivity +
  transitiveValue + cycleValue + commonSourceValue + commonTargetValue
```

The ML baselines are supervised classifiers over comparable transaction-history
features:

- logistic regression (`glm_history_features`)
- classification tree if `rpart` is installed (`rpart_history_features`)

## Metrics

The current runner reports:

- ROC AUC
- average precision
- number of test transactions
- number of SAR positives
- exact directed dyad risk-set size
- fit and scoring timings

## Exact Scalability

Use `amlsim_run_exact_scalability_grid()` to run increasing actor/event caps.
The exact risk set is `n_actors * (n_actors - 1)` per event, so large AMLSim
runs are expected to become expensive quickly. For the valued formula above,
HawkesNet now uses ERNM's fast Rcpp valued-hit backends, which compute exact
candidate changes and streaming log probabilities from the valued adjacency
matrix rather than recalculating global statistics candidate by candidate. On
the timestamped simulator data, the 50- and 99-node/1000-event exact runs
completed in about 19 and 25 seconds end-to-end with `AMLSIM_HAWKES_MAXIT=10`
and full-intensity scoring. The full timestamped ML baseline completed on
122,108 events in about 20 seconds; full exact HawkesNet remains too expensive
because one likelihood pass over the 72,686 non-SAR training events is estimated
to take hours.
