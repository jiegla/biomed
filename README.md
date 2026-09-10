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
| `batch_ANOVA()` | `batch_anova()` |
| `chi_square_batch_df()` | `batch_chi_square()` |
| `cophx_batch()` | `batch_cox()` |
| `roc_batch()` | `batch_roc()` |
| `plot_roc_batch()` | `plot_roc()` |
| `donutPie()` | `plot_donut()` |
| `draw_stack_barplot()` | `plot_stacked_bar()` |
| `vennjgl()` | `plot_venn()` |

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
| `batch_ANOVA()` | Batch one-way ANOVA with BH-adjusted p-values |
| `chi_square_batch_df()` | Batch chi-square/Fisher association tests |
| `cophx_batch()` | Batch univariable Cox regression |
| `roc_batch()` | Batch ROC metrics and confidence intervals |
| `plot_roc_batch()` | ROC figures returned as reusable ggplot objects |
| `donutPie()` | Donut charts for categorical variables |
| `draw_stack_barplot()` | Proportional stacked bars with association tests |
| `vennjgl()` | Venn diagrams and intersection-member tables |

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
donutPie(grouped, "biomarker_binary")
```

ROC, maximally selected survival cutoffs, Venn diagrams, and XLSX output use
suggested packages that are installed on demand by users who need those
features. Functions do not require patient identifiers, and examples and tests
use synthetic data only.

## Survival and correlation tools (0.2.0)

| Original function | Standard interface | Main purpose |
| --- | --- | --- |
| `best_cutoff_jgl()` | `find_survival_cutoff()` | Survival cutoff with valid-group checks and median fallback |
| `batch_surv_jgl()` | `batch_survival()` | Batch Cox models, eligibility thresholds, failures and FDR |
| `surv_fig_hr()` | `plot_survival()` | Kaplan-Meier curves, risk table, medians and directed HR |
| Internal pairwise helper | `pairwise_survival()` | Named pairwise HRs and log-rank tests |
| `get_cor_jgl()` | `plot_correlation()` | Correlation statistics, scatter plots, labels and exports |
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
with `output_dir = NULL`. Legacy `get_cor_jgl(save_plot = TRUE, path = ...)`
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
