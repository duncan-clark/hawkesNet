# Cluster Compatibility Checklist for OpenAlex Study

## ✅ Verified Items

### 1. Package Loading
- ✅ Package loading with fallback to `devtools::load_all()` for development
- ✅ Error handling if package not found
- ✅ All required functions exported in NAMESPACE:
  - `PMF_mark_CS` ✓
  - `cond_intensity_inhom` ✓
  - `gof()` ✓

### 2. Dependencies
- ✅ Core dependencies: network, sna, ernm, parallel (in DESCRIPTION)
- ✅ ggplot2: Optional, checked with `requireNamespace()` before use
- ✅ dplyr: Required, checked before loading
- ✅ get_network_openalex.R dependencies: httr, jsonlite, dplyr, purrr, tidyr

### 3. File Paths
- ✅ Uses `PKG_ROOT <- getwd()` (works from package root)
- ✅ `cluster_output/` directory created with `dir.create(..., recursive = TRUE)`
- ✅ Relative paths: `file.path(PKG_ROOT, "cluster_output", ...)`
- ✅ SLURM script creates `cluster_output/` directory

### 4. Environment Variables
- ✅ `SLURM_CPUS_PER_TASK` used for `N_CORES` (defaults to 7)
- ✅ `OPENALEX_EMAIL`, `OPENALEX_STRING`, `OPENALEX_PAGES` can be set via env

### 5. Error Handling
- ✅ All major sections wrapped in `tryCatch()`
- ✅ GOF function handles missing ggplot2 gracefully
- ✅ Paper output checks for ggplot2 before plotting
- ✅ Falls back to legacy plotting if gof() plots not available

### 6. Parallelization
- ✅ Uses `parallel::mclapply()` (works on cluster)
- ✅ Respects `N_CORES` from SLURM
- ✅ All GOF statistics computed in parallel

### 7. Output Files
- ✅ Saves to `cluster_output/results_openalex_full.RDS`
- ✅ SLURM output/error go to `cluster_output/openalex_%j.out/err`
- ✅ All outputs use relative paths

## ⚠️ Potential Issues to Verify

### 1. Package Installation on Cluster
- **Action**: Ensure `hawkesGrowthNet` is installed on cluster
- **Check**: Run `Rscript -e "library(hawkesGrowthNet)"` before submitting job

### 2. R Module Version
- **Action**: Verify R module version is compatible
- **Check**: `module load R` should load appropriate version

### 3. Network Access
- **Action**: Ensure cluster nodes can access OpenAlex API (api.openalex.org)
- **Check**: Test with `curl https://api.openalex.org/works?per_page=1`

### 4. Memory Requirements
- **Current**: 32G requested
- **Check**: Monitor actual usage; may need adjustment for full 100 pages

### 5. Time Limit
- **Current**: 24 hours
- **Estimate**: 
  - Data fetch: ~5-10 min
  - Fitting (2-stage): ~30-60 min
  - GOF (50 sims): ~1-2 hours
  - **Total**: Should complete well within 24 hours

### 6. ggplot2 Installation
- **Status**: Optional but recommended for plots
- **Action**: Install ggplot2 if plots needed: `install.packages("ggplot2")`

## 📋 Pre-Submission Checklist

- [ ] Package `hawkesGrowthNet` is installed on cluster
- [ ] R module is loaded (`module load R`)
- [ ] `cluster_output/` directory exists or will be created
- [ ] Network access to OpenAlex API is available
- [ ] Optional: ggplot2 installed for plots
- [ ] SLURM script paths are correct (run from package root)
- [ ] Environment variables set if needed (OPENALEX_EMAIL, etc.)

## 🚀 Submission Command

```bash
cd /path/to/hawkesGrowthNet
mkdir -p cluster_output
sbatch inst/openalex_study/run_openalex.slurm
```

## 📊 Expected Output Files

1. `cluster_output/openalex_<JOBID>.out` - Standard output
2. `cluster_output/openalex_<JOBID>.err` - Error output  
3. `cluster_output/results_openalex_full.RDS` - Full results (if successful)
