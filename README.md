# biomed

<!-- badges: start -->
[![R-CMD-check](https://github.com/jiegla/biomed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jiegla/biomed/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`biomed` is an R package for reusable, auditable clinical and biomedical
research workflows. Its planned scope includes cohort preparation, efficacy
and survival endpoints, biomarker exploration, omics analysis, and
publication-ready tables and figures.

> [!IMPORTANT]
> `biomed` is for research use only. It is not medical advice and is not a
> validated clinical decision-support system.

## Built-in colors and survival export

`biomed_colors()` returns the 14 package colors in their supplied order,
including the `CC` alpha channel (80% opacity). `biomed_colors(n)` selects
colors, or interpolates when more than 14 are needed. Default categorical
plots use this palette; explicit colors and named ggsci palettes still work.

```r
biomed_colors()
# p <- surv_fig_hr("group", d, time = "time", status = "status",
#                  output_dir = NULL)
# ggplot2::ggsave("survival.pdf", plot = p, width = 8, height = 7)
# ggplot2::ggsave("survival.png", plot = p, width = 8, height = 7, dpi = 300)
```

The package registers `grid.draw.ggsurvplot()` automatically. Saving the whole
`p` object includes its risk table; saving `p$plot` saves only the curve panel.
You do not need to define a helper in your workspace.

## Function naming (0.3.0)

Public functions each have a matching `R/function_name.R` file. Shared internal
helpers remain in `R/internal_utils.R`; package documentation remains in
`R/biomed-package.R`. Existing standardized interfaces such as `batch_anova`,
`batch_survival`, `plot_survival` and `plot_correlation` keep their names.

The following functions were renamed without changing their arguments or
calculations. Update old scripts using this table; the old exports are removed.

| Previous name | New name |
| --- | --- |
| `best_cutoff_jgl()` | `best_cutoff()` |
| `batch_surv_jgl()` | `batch_surv()` |
| `get_cor_jgl()` | `get_cor()` |
| `vennjgl()` | `venn_plot()` |
| `batch_ANOVA()` | `anova_batch()` |
| `donutPie()` | `donut_pie()` |
| `cophx_batch()` | `cox_batch()` |

`anova_batch.R` and `batch_anova.R` now also work on case-insensitive file systems.

## Design principles

- Reproducible analyses with explicit inputs, endpoints, and assumptions.
- Safe defaults for patient-level data and credentials.
- Publication-ready outputs in PDF/PNG and XLSX where appropriate.
- Small, tested functions that can be combined into study-specific pipelines.
- Compatibility with established R ecosystems rather than reimplementing them.

## Installation

Install the development version from GitHub:

```r
pak::pak("jiegla/biomed")
```

## Standardized API

New code should use the interfaces below. All original entry points remain
exported and retain their result column names.

| Original name | Standard name |
|---|---|
| `makegroup()` | `make_group()` |
| `anova_batch()` | `batch_anova()` |
| `chi_square_batch_df()` | `batch_chi_square()` |
| `cox_batch()` | `batch_cox()` |
| `roc_batch()` | `batch_roc()` |
| `plot_roc_batch()` | `plot_roc()` |
| `donut_pie()` | `plot_donut()` |
| `draw_stack_barplot()` | `plot_stacked_bar()` |
| `venn_plot()` | `plot_venn()` |

The new interfaces consistently use `data`, `variables`, `output_dir`,
`colors` and `remove_na` where applicable. Analysis tables use
`variable`, `p_value`, `p_adjust` and other snake_case columns.
Grouping returns augmented data; donut and ROC plotting return named plot lists;
stacked bars return `plots` and `statistics`. Venn returns its drawing,
partitions and file paths.

`output_dir = NULL` writes no files. Supply a directory to export batch
analysis tables as XLSX and figures as PDF/PNG (Venn also supports TIFF).
Existing output filenames are overwritten. Options not renamed by a new wrapper
can still be passed using the original parameter names through `...`.

```r
d <- data.frame(
  group = rep(c("A", "B"), each = 6),
  marker = c(1:6, 5:10),
  response = rep(c("No", "Yes"), 6)
)
batch_anova(d, variables = "marker", group = "group")
batch_roc(d, "marker", response = "response", positive_class = "Yes")
plot_stacked_bar(d, variables = "response", group = "group")
# batch_anova(d, "marker", output_dir = "results") # XLSX
```

### Compatibility and statistical choices

- The nine original names remain callable. New wrappers rename columns; the
  original functions retain their historical package column names.
- `feature_manipulation = TRUE` again filters incomplete, non-numeric,
  infinite and constant features. It no longer needs IOBR.
- Three-group cutoffs use the original 0.33/0.66 quantiles. Tied cutoffs keep all
  rows and the Low/Middle/High levels; some groups can legitimately be empty.
- Automatic categorical tests use Fisher whenever any expected cell count is
  below five. Large exact tests can fail; batch results retain the error.
- Cox requires 0 = censored and 1 = event. Multi-level factors return all contrasts.
- Standard ROC interfaces infer the positive class only when omitted and record
  it in the result. Set it explicitly in research scripts. The original
  `roc_batch()` retains its default positive class of `"1"`.
- ROC plots label false-positive rate correctly. Tied optimal cutoffs select the
  first finite observed optimum; fixed-sensitivity/specificity grouping selects
  a realizable threshold meeting at least the target, rather than interpolating.
  Low/High labels always describe predictor magnitude, regardless of ROC direction.
- ROC confidence intervals for perfect separation can warn; this is retained.
  Optimized thresholds and same-data AUCs require independent validation.
- Two to five sets use Venn geometry. Six sets use an exact membership matrix
  because the VennDiagram backend supports at most five sets.
- Legacy stacked bars save only with `save = TRUE`; their original CSV export
  remains available through `statistics_format = "csv"`. XLSX is the default.
  Legacy Venn retains its automatic file output; new plotting wrappers do not
  write unless `output_dir` is supplied.

## Included utilities

| Function | Purpose |
|---|---|
| `makegroup()` | Reproducible median, survival, or ROC-derived cutoffs |
| `anova_batch()` | Batch one-way ANOVA with BH-adjusted p-values |
| `chi_square_batch_df()` | Batch chi-square/Fisher association tests |
| `cox_batch()` | Batch univariable Cox regression |
| `roc_batch()` | Batch ROC metrics and confidence intervals |
| `plot_roc_batch()` | ROC figures returned as reusable ggplot objects |
| `donut_pie()` | Donut charts for categorical variables |
| `draw_stack_barplot()` | Proportional stacked bars with association tests |
| `venn_plot()` | Venn diagrams and intersection-member tables |

## Quick start

```r
library(biomed)

cohort <- data.frame(
  patient_id = c("P001", "P002"),
  response = c("PR", "SD")
)

biomed_check_columns(cohort, c("patient_id", "response"))

analysis_data <- data.frame(
  response = factor(c("No", "No", "Yes", "Yes")),
  biomarker = c(0.8, 1.2, 2.7, 3.1)
)

grouped <- makegroup(analysis_data, "biomarker", method = "median")
donut_pie(grouped, "biomarker_binary")
```

ROC, maximally selected survival cutoffs, Venn diagrams, and XLSX output use
suggested packages that are installed on demand by users who need those
features. Functions do not require patient identifiers, and examples and tests
use synthetic data only.

## Survival and correlation tools (0.2.0)

| Original function | Standard interface | Main purpose |
| --- | --- | --- |
| `best_cutoff()` | `find_survival_cutoff()` | Survival cutoff with valid-group checks and median fallback |
| `batch_surv()` | `batch_survival()` | Batch Cox models, eligibility thresholds, failures and FDR |
| `surv_fig_hr()` | `plot_survival()` | Kaplan-Meier curves, risk table, medians and directed HR |
| Internal pairwise helper | `pairwise_survival()` | Named pairwise HRs and log-rank tests |
| `get_cor()` | `plot_correlation()` | Correlation statistics, scatter plots, labels and exports |
| `sanitize_filename()` | `sanitize_filename()` | Safe output file names without a pipe dependency |

The new survival interfaces default to **0 = censored, 1 = event**. Specify
`status_encoding = "12"` for 1 = censored, 2 = event, or `"labels"` for words
such as `alive`/`dead`. Legacy interfaces retain automatic detection, which
treats an all-1 cohort as all events; use an explicit encoding for ambiguous data.
HRs always describe the named group relative to the named reference. No HR is
inverted merely because it exceeds 1. All fits use the full eligible follow-up;
`max_time` only limits the plot display.

```r
set.seed(42)
d <- data.frame(
  time = rexp(120), status = rep(c(1, 1, 0), 40),
  group = factor(rep(c("A", "B"), 60)), marker = rnorm(120)
)
analysis <- batch_survival(d, c("marker", "group"))
analysis$results
analysis$failed

# Optional packages: install.packages(c("survminer", "ggsci", "openxlsx"))
km <- plot_survival(d, "group", reference = "A")
print(km)
attr(km, "statistics")
pairwise_survival(d, "group", reference = "A")

p <- plot_correlation(d, "marker", "time", scale = FALSE)
attr(p, "correlation")
print(p)

# Optional exports:
# batch_survival(d, "marker", output_dir = "results")
# plot_survival(d, "group", output_dir = "results")
```

`plot_survival()` exports the curve **and risk table** together as PDF/PNG plus
a statistics text file. New interfaces write nothing unless `output_dir` is
provided. Legacy `surv_fig_hr()` retains its `Surv_Output` text log; disable it
with `output_dir = NULL`. Legacy `get_cor(save_plot = TRUE, path = ...)`
retains the analyzed-data RData export; disable it with `save_data = FALSE`.
The new `plot_correlation()` requires `save_data = TRUE` to export data.

Correlation removes original zero values and missing groups before scaling or
testing, so the reported sample matches the plot. Grouped regression lines
remain group-specific; the annotated correlation is for all retained rows.
Cutoff selection is exploratory: downstream Cox P values and confidence
intervals do not account for selecting the cutoff on the same data, and BH
correction does not remove this selection bias.

## Roadmap modules

- `clinical`: cohort cleaning, endpoint derivation, response and survival.
- `biomarker`: subgroup, interaction, multivariable, and validation workflows.
- `genomics`: mutation profiles, MAF utilities, and biomarker screening.
- `singlecell`: sample-aware summaries and pathway-score comparisons.
- `reporting`: journal-ready tables, figures, and reproducibility manifests.

The exact public API will be developed incrementally from validated study
workflows.

## Patient data and secrets

Never commit identifiable or patient-level source data. Keep protected data
outside the repository, use local paths supplied at runtime, and commit only
synthetic examples. Store API credentials in environment variables or GitHub
Actions secrets; never place tokens in R scripts, notebooks, configuration
files, examples, or logs. See [Security and data handling](SECURITY.md).

## Development

```r
pak::local_install_deps(dependencies = TRUE)
devtools::document()
devtools::test()
devtools::check()
```

## License

MIT © 2026 Guangling Jie

## Single-cell and Seurat workflows

All public single-cell functions start with `sce_`; each R filename matches its function.

| Original script/function | biomed function |
| --- | --- |
| add_seurat_metadata | sce_add_seurat_metadata |
| analyze_annotated_seurat | sce_analyze_annotated_seurat |
| compare_gene_expression_groups | sce_compare_gene_expression_groups |
| run_sce_scoring | sce_run_scoring |
| run_scenic script | sce_run_scenic |
| sce_cell_cycle | sce_cell_cycle |

```r
sce <- biomed::sce_add_seurat_metadata(sce, "clinical.csv", file_encoding = "UTF-8")
sce <- biomed::sce_cell_cycle(sce, species = "human", assay = "RNA")
sce <- biomed::sce_run_scoring(sce, gene_sets, method = "UCell")
comparison <- biomed::sce_compare_gene_expression_groups(
  sce, genes = c("CD3D", "MS4A1"), statistic_group = "group",
  split_by_celltype = "celltype", statistic_level = "samples")
summary <- biomed::sce_analyze_annotated_seurat(
  sce, celltype = "celltype", clinical_var = "group")
```

Install Seurat and the reporting dependencies for these workflows. Scoring uses optional
Bioconductor packages `UCell`, `GSVA` (>= 1.50.0), `AUCell`, and `BiocParallel`:

```r
install.packages(c("Seurat", "dplyr", "tidyr", "tidyselect", "readxl", "readr",
                   "patchwork", "pheatmap", "gplots", "openxlsx", "withr"))
BiocManager::install(c("UCell", "GSVA", "AUCell", "BiocParallel"))
```

Metadata files use the first column as sample IDs. All Seurat samples must match by default;
extra file samples are ignored. `strict_match = TRUE` requires equal sample sets.
`allow_missing = TRUE` explicitly restores the uploaded script's NA-fill behavior.
Existing columns receive `.new` suffixes unless `replace_sample = TRUE`.
CSV leading-zero identifiers are preserved. `cols` retains tidyselect syntax.

Expression comparisons default to **samples as independent units**. Cell-level tests remain
available, but do not account for correlations between cells from the same sample.
The reporting functions retain PDF/PNG figures, XLSX tables and the contributed CSV/RDS
analysis exports. The annotation workflow defaults to `biomed_colors()`.
See each function's R help for marker, correlation, export and aggregation controls.

### SCENIC requirements

Install [SCENIC](https://github.com/aertslab/SCENIC) and its dependencies separately, and
obtain species-compatible cisTarget ranking databases before calling:

```r
remotes::install_github(c("aertslab/RcisTarget", "aertslab/SCENIC"))
result <- biomed::sce_run_scenic(sce, db_dir = "/path/to/cisTarget",
  dbs = c("compatible-ranking-database.feather"), species = "human",
  output_dir = "scenic_output", n_cores = 8)
sce <- result$sce
```

SCENIC exports regulon AUC, the updated object, metadata-name mapping, feature plots and
cell-type mean AUC heatmaps. It never runs while loading biomed. The database directory
and expression layer are explicit; matching gene symbols and motif annotations are required.
The full workflow can require substantial memory because SCENIC uses a dense expression matrix.
CI checks input handling; full regulatory-network validation requires external ranking
and annotation resources and is not covered by the small synthetic tests.
GSVA uses the current [parameter-object API](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html).


## Additional single-cell and spatial workflows

These six contributed workflows use matching `sce_` function and R filenames.
The existing six single-cell functions are unchanged.

| Function | Preserved outputs |
| --- | --- |
| `sce_calculate_cell_fraction()` | Cell counts/fractions and default XLSX export |
| `sce_gene_expression_statistic()` | Sample/cell summaries, omnibus/pairwise tests, XLSX/CSV, RDS, PNG/PDF |
| `sce_run_two_group_deg_enrichment()` | Seurat DEG, GMT ORA/GSEA, per-cell-type Excel/RDA objects and plots |
| `sce_run_scTenifoldKnk_KO()` | Single/multiple virtual KOs, perturbation tables, figures, TXT summaries and RDS bundles |
| `sce_plot_spatial_continuous()` | Smoothed expression/metadata surfaces, H&E, contours, masks and plot data |
| `sce_plot_pathway_rds_continuous()` | Batch pathway maps across RDS files/images, summary and plotting logs |

The uploaded `vis_plot_*.R` scripts are exposed as `sce_plot_*` to follow the
package's single-cell prefix convention. Existing argument names and default
exports are retained. Use `?function_name` for full parameter documentation.

```r
fractions <- sce_calculate_cell_fraction(
  sce, sample_col = "orig.ident", celltype_col = "celltype_major"
)
expression <- sce_gene_expression_statistic(
  sce, genes = c("CD3D", "MS4A1"), sample_col = "orig.ident",
  celltype_col = "celltype_major", group_col = "condition",
  output_dir = "expression_report"
)
deg <- sce_run_two_group_deg_enrichment(
  sce, group_col = "condition", ident_1 = "treated", ident_2 = "control",
  celltype_col = "celltype_major", output_dir = "deg_report",
  run_ora = FALSE, run_gsea = FALSE
)
# To enable ORA/GSEA, supply gmt_files or gmt_dir and enable the corresponding flags.
# Example configuration:
system.file("examples", "sce_deg_config.yml", package = "biomed")

ko <- sce_run_scTenifoldKnk_KO(
  sce, gKO = "CASP4", subset_col = "celltype_major",
  subset_values = "Macrophage", layer = "counts", nCores = 2
)
p <- sce_plot_spatial_continuous(
  spatial, feature = "CD3D", colors = biomed_colors(),
  save_path = "spatial/CD3D.png"
)
maps <- sce_plot_pathway_rds_continuous(
  spatial, score_dir = "gmt_pathway_scores", keep_plots = TRUE
)
```

Dependencies are optional and checked at use time. Install `scTenifoldKnk` for
virtual KO, `FNN`/`viridisLite` for spatial smoothing, `writexl` for the
gene-expression Excel report (otherwise CSV), and `openxlsx` for the other Excel
exports. GMT enrichment requires `clusterProfiler` from Bioconductor;
gene-ID conversion additionally needs `AnnotationDbi` and the configured OrgDb.
YAML and QS input require `yaml` and `qs`, respectively.

Sample-level gene-expression tests use one summary per sample. Cell-level tests
and Seurat FindMarkers do not model within-sample dependence. With
`remove_zero = TRUE`, filtering precedes summaries and minimum-cell thresholds.
Virtual KO requires raw counts; its FC statistic is a perturbation measure, not
signed expression log2FC. The original QC arguments are adapted to the
[scTenifoldKnk API](https://cran.r-project.org/web/packages/scTenifoldKnk/refman/scTenifoldKnk.html);
version 1.1 uses a cell fraction for the minimum gene-detection QC threshold.
Enrichr requests occur only when enrichment or network annotation is enabled.
Inspect returned DEG status tables and pathway plotting logs for skipped/failed steps.

