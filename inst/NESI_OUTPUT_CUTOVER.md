# NeSI output cutover (zero-loss)

Package code no longer writes into `PKG_ROOT/cluster_output`. Use an external directory.

## Target paths

| Role | Path |
|------|------|
| Package clone | `/nesi/project/uoo04008/Duncan/hawkes_net/hawkesNet` |
| Durable outputs | `/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output` |
| Env var | `export HAWKESNET_OUTPUT_DIR=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output` |

## Cutover steps (copy first, delete never until verified)

```bash
# 1) Create durable store
mkdir -p /nesi/project/uoo04008/Duncan/hawkes_net/cluster_output/{runs,logs,paper_figures,diagnostics}

# 2) Union-copy from old in-clone location (do not delete source yet)
OLD=/nesi/project/uoo04008/Duncan/hawkes_net/hawkesNet/cluster_output
NEW=/nesi/project/uoo04008/Duncan/hawkes_net/cluster_output
rsync -avh --ignore-existing "$OLD"/ "$NEW"/

# 3) Inventory / checksum both trees
(cd "$OLD" && find . -type f ! -name '.*' | sort | xargs shasum -a 256) > /tmp/old_sha.txt
(cd "$NEW" && find . -type f ! -name '.*' | sort | xargs shasum -a 256) > /tmp/new_sha.txt
# Manually inspect any basename conflicts (same name, different hash).
# Keep both under distinct names if needed.

# 4) Point jobs at the new dir
export HAWKESNET_OUTPUT_DIR="$NEW"
# or: source /path/to/hawkesNet/inst/cluster_env.sh

# 5) Update package clone to main / v0.1.0-jasa
cd /nesi/project/uoo04008/Duncan/hawkes_net/hawkesNet
git fetch origin
git checkout main
git pull --ff-only origin main
git checkout v0.1.0-jasa   # optional pinned tag for reproduction

# 6) Only after local laptop sync also has the files:
#    mv "$OLD" "${OLD}__archived_YYYYMMDD"
```

## Laptop sync

rsync/scp the durable NeSI `cluster_output` into the project-root folder:

`.../0003_hawkes_net/cluster_output`

Do not put that folder inside the package git repo.
