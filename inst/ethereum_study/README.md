# Ethereum Stablecoin HawkesNet Study

This directory contains a public-data pipeline for applying valued repeated-hit
HawkesNet to Ethereum stablecoin transfers.

The target event format is:

- `t`: strict event time in `[0, 1]`
- `i`, `j`: sender and receiver entity IDs
- `weight`: rescaled `log1p(amount)`
- `amount`: token amount after token decimals
- `token_symbol`: `USDC`, `USDT`, or `DAI`

## Why This Data

Stablecoin transfers are public, timestamped, directed, and valued financial
transactions. Collapsing addresses into labeled entities gives a network of
institution-like actors: exchanges, treasuries, DeFi protocols, bridges,
contracts, and high-activity wallets. This makes the data a practical public
analogue to proprietary interbank transaction networks such as e-MID.

## Data Access

The preferred public source is Google BigQuery:

`bigquery-public-data.crypto_ethereum.token_transfers`

Relevant columns are `block_timestamp`, `transaction_hash`, `log_index`,
`token_address`, `from_address`, `to_address`, and `value`.

The utilities query canonical Ethereum mainnet stablecoins:

- USDC: `0xA0b86991c6218b36c1d19D4a2e9Eb0cE3606eB48`
- USDT: `0xdAC17F958D2ee523a2206206994597C13D831ec7`
- DAI: `0x6B175474E89094C44Da98b954EedeAC495271d0F`

BigQuery requires a Google Cloud project for billing public-data queries. Set:

```bash
GCP_PROJECT=my-gcp-project
ETH_START_DATE=2024-03-01
ETH_END_DATE=2024-03-08
ETH_MAX_ENTITIES=1000
ETH_MAX_EVENTS=10000
Rscript inst/ethereum_study/run_ethereum_study.R
```

If you already exported a CSV with the BigQuery schema, bypass BigQuery:

```bash
ETH_TRANSFERS_CSV=/path/to/token_transfers.csv \
ETH_MAX_ENTITIES=1000 \
ETH_MAX_EVENTS=10000 \
Rscript inst/ethereum_study/run_ethereum_study.R
```

When BigQuery credentials are unavailable, a smaller public-RPC scrape is also
available:

```bash
ETH_RPC_URL=https://eth-mainnet.public.blastapi.io \
ETH_RPC_BLOCK_WINDOW=200 \
ETH_STUDY_OUTPUT_DIR=ethereum_study_results_rpc \
python3 inst/ethereum_study/scrape_stablecoin_rpc.py
```

Some public RPC providers restrict `eth_getLogs` ranges. The scraper defaults to
10-block chunks via `ETH_RPC_LOG_CHUNK=10`.

## Labels and Node Covariates

The pipeline can use public Etherscan label snapshots from:

`https://github.com/brianleect/etherscan-labels`

Set `ETH_DOWNLOAD_LABELS=true` to download the combined Ethereum JSON label file
into the output directory, or pass your own entity file:

```bash
ETH_LABELS_CSV=/path/to/address_labels.csv \
Rscript inst/ethereum_study/run_ethereum_study.R
```

Expected label columns are flexible but should include an address column and,
when possible, entity/name/category columns. Unlabeled addresses are kept as
address-level residual actors unless filtered out by the active-entity cap.

The node covariate output includes:

- category/entity label when available
- in/out transaction counts
- in/out stablecoin amounts
- total activity and value
- net received amount
- activity share
- first/last observed activity time and active span

These are saved to `ethereum_node_covariates.csv`.

## Outputs

By default outputs are written to `ethereum_study_results/`:

- `ethereum_transfers_prepared.csv`
- `ethereum_hits.csv`
- `ethereum_actor_map.csv`
- `ethereum_node_covariates.csv`
- `ethereum_summary.csv`
- `ethereum_hawkesnet_data.rds`
- `ethereum_model_scores.csv`
- `ethereum_model_metrics.csv`

By default the runner computes cheap held-out mark-prediction baselines:

- dyad frequency
- dyad recency
- sender-receiver activity
- fitted history softmax

Set `ETH_RUN_RHEM=true` to fit the no-excitation HawkesNet special case:

```bash
ETH_RUN_RHEM=true \
ETH_HAWKES_MAXIT=10 \
ETH_MAX_ENTITIES=200 \
ETH_MAX_EVENTS=2000 \
Rscript inst/ethereum_study/run_ethereum_study.R
```

This model fixes `K = 0`, so the ground process is a homogeneous Poisson
process, but the mark distribution is still `PMF_mark_RHEM`. It is the direct
RHEM comparison nested inside HawkesNet.

Set `ETH_RUN_HAWKES=true` to run the exact valued HawkesNet fit with temporal
self-excitation:

```bash
ETH_RUN_HAWKES=true \
ETH_HAWKES_MAXIT=10 \
ETH_MAX_ENTITIES=200 \
ETH_MAX_EVENTS=2000 \
Rscript inst/ethereum_study/run_ethereum_study.R
```

Exact runtime scales with `n_events * n_entities * (n_entities - 1)`. Start with
small caps and increase gradually.

Set `ETH_RHEM_FORMULA=simple` to use the built-in repeated-hit features
`repetition`, `reciprocity`, `sender_activity`, and `receiver_activity`. The
default `ETH_RHEM_FORMULA=valued` uses the full valued ERNM formula terms.

Set `ETH_PREFER_LABELED=true` only when labels are more important than activity
coverage. The default keeps the most active entities first and attaches labels
when available.

## Scientific Target

The primary evaluation is held-out mark likelihood:

Given that a stablecoin transfer happened at time `t`, how well does HawkesNet
predict the sender-receiver entity pair using valued repeated-transfer history?

This separates the Hawkes ground process, which captures bursts in stablecoin
activity, from the network mark model, which captures institutional routing,
reciprocity, fan-in/fan-out, repeated exchange flows, and liquidity cycles.
