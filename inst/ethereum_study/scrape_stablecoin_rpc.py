#!/usr/bin/env python3
"""Scrape Ethereum stablecoin Transfer logs through JSON-RPC.

This is a lightweight fallback when BigQuery credentials are unavailable. It
uses eth_getLogs in small block chunks and writes the CSV schema consumed by
run_ethereum_study.R.
"""

import csv
import json
import os
import time
import urllib.request


RPC_URL = os.environ.get("ETH_RPC_URL", "https://eth-mainnet.public.blastapi.io")
BLOCK_WINDOW = int(os.environ.get("ETH_RPC_BLOCK_WINDOW", "200"))
CHUNK_SIZE = int(os.environ.get("ETH_RPC_LOG_CHUNK", "10"))
OUTDIR = os.environ.get("ETH_STUDY_OUTPUT_DIR", "ethereum_study_results_rpc")
OUTFILE = os.environ.get(
    "ETH_RPC_TRANSFERS_CSV",
    os.path.join(OUTDIR, "ethereum_stablecoin_rpc_transfers.csv"),
)

TOKENS = {
    "USDC": "0xa0b86991c6218b36c1d19d4a2e9eb0ce3606eb48",
    "USDT": "0xdac17f958d2ee523a2206206994597c13d831ec7",
    "DAI": "0x6b175474e89094c44da98b954eedeac495271d0f",
}
TRANSFER_TOPIC = "0xddf252ad1be2c89b69c2b068fc378daa952ba7f163c4a11628f55a4df523b3ef"


def rpc(method, params, tries=4):
    payload = {"jsonrpc": "2.0", "id": 1, "method": method, "params": params}
    data = json.dumps(payload).encode()
    last_error = None
    for attempt in range(tries):
        request = urllib.request.Request(
            RPC_URL,
            data=data,
            headers={
                "Content-Type": "application/json",
                "User-Agent": "hawkesNet-study/0.1",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                answer = json.loads(response.read().decode())
            if "error" in answer:
                raise RuntimeError(answer["error"])
            return answer["result"]
        except Exception as exc:  # noqa: BLE001 - command-line retry wrapper
            last_error = exc
            time.sleep(0.5 * (attempt + 1))
    raise last_error


def decode_transfer_log(token_address, log):
    topics = log.get("topics", [])
    if len(topics) < 3:
        return None
    return {
        "block_number": int(log["blockNumber"], 16),
        "transaction_hash": log["transactionHash"],
        "log_index": int(log["logIndex"], 16),
        "token_address": token_address,
        "from_address": ("0x" + topics[1][-40:]).lower(),
        "to_address": ("0x" + topics[2][-40:]).lower(),
        "raw_value": str(int(log.get("data", "0x0"), 16)),
    }


def main():
    latest = int(rpc("eth_blockNumber", []), 16)
    from_block = max(0, latest - BLOCK_WINDOW + 1)
    print(f"rpc {RPC_URL} latest {latest} from {from_block} window {BLOCK_WINDOW}")

    rows = []
    for symbol, token_address in TOKENS.items():
        token_count = 0
        for start in range(from_block, latest + 1, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE - 1, latest)
            logs = rpc(
                "eth_getLogs",
                [{
                    "fromBlock": hex(start),
                    "toBlock": hex(end),
                    "address": token_address,
                    "topics": [TRANSFER_TOPIC],
                }],
            )
            token_count += len(logs)
            for log in logs:
                row = decode_transfer_log(token_address, log)
                if row is not None:
                    rows.append(row)
        print(f"{symbol} {token_count} logs")

    rows.sort(key=lambda row: (row["block_number"], row["log_index"]))
    blocks = sorted({row["block_number"] for row in rows})
    print(f"total_logs {len(rows)} blocks {len(blocks)}")

    block_timestamps = {}
    for block_number in blocks:
        block = rpc("eth_getBlockByNumber", [hex(block_number), False])
        block_timestamps[block_number] = time.strftime(
            "%Y-%m-%d %H:%M:%S",
            time.gmtime(int(block["timestamp"], 16)),
        )
    for row in rows:
        row["block_timestamp"] = block_timestamps[row["block_number"]]

    os.makedirs(os.path.dirname(OUTFILE) or ".", exist_ok=True)
    fields = [
        "block_timestamp",
        "transaction_hash",
        "log_index",
        "block_number",
        "token_address",
        "from_address",
        "to_address",
        "raw_value",
    ]
    with open(OUTFILE, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})
    print(f"wrote {OUTFILE}")


if __name__ == "__main__":
    main()
