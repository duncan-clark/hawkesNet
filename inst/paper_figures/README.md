# Paper figures

JASA manuscript PDFs are **not** stored as large RDS in git.

- **Manuscript copies:** `paper_1/assets/` (local project, outside this package)
- **Synced evidence copies:** `$HAWKESNET_OUTPUT_DIR/paper_figures/`
- **Lightweight BA set in-repo:** `paper_output/` (boxplots / RMSE / explosive)

Regenerate Hypertext GOF (title-stripped) with:

```bash
export HAWKESNET_OUTPUT_DIR=/path/to/cluster_output
Rscript inst/hypertext_conference/save_gof_plots_paper.R
```

See the project-level `cluster_output/FIGURE_MANIFEST.md` for the full figure → RDS → script map.
