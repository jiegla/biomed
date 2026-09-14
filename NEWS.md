# biomed 0.4.2

- Include the survival endpoint column before the grouping variable in timestamped
  `surv_fig_hr()` text filenames, for example `PFS_treatment_Output.txt`.

# biomed 0.4.1

- Report empty and missing survival column names explicitly in `surv_fig_hr()`
  and shared column validation, including the available data columns.

# biomed 0.4.0

- Add six single-cell workflows with consistent `sce_` names and matching R files.
- Preserve contributed analysis reports, figures, sample-level tests and cell-level options.
- Add explicit CSV encoding and sample-matching controls, preserving leading-zero IDs.
- Align expression and score matrices by cell names and respect assay/layer selection.
- Update GSVA scoring to its parameter-object API; support UCell and AUCell.
- Convert SCENIC from a fixed-path script to a function with configurable databases,
  correct filtering and correlation steps, and unique regulon metadata names.

# biomed 0.3.1

* Add biomed_colors() with the supplied 14 RGBA colors, preserving CC alpha.
  Default categorical palettes, including correlation and Venn plots, use it.
  Explicit custom colors and named ggsci palettes remain available.
* Register grid.draw.ggsurvplot in surv_fig_hr.R so grid.draw and ggsave can
  draw complete survival plots with risk tables without a global helper.

# biomed 0.3.0

* Give each exported function a matching R source filename; split bundled APIs
  and keep shared private helpers in internal_utils.R.
* Remove personal suffixes: best_cutoff_jgl -> best_cutoff,
  batch_surv_jgl -> batch_surv, get_cor_jgl -> get_cor, vennjgl -> venn_plot.
* Normalize batch_ANOVA -> anova_batch, donutPie -> donut_pie, and
  cophx_batch -> cox_batch. Avoid case-only source filename collisions.
* Update calls, exports, help pages and tests. Function arguments and algorithms
  are unchanged. Old renamed exports are removed; see README migration table.

# biomed 0.2.0

* Added survival cutoff selection, screened batch Cox analysis, Kaplan-Meier
  curves with risk tables, directed pairwise survival comparisons, correlation
  plots, and portable filename sanitization.
* Preserved best_cutoff_jgl, batch_surv_jgl, surv_fig_hr, get_cor_jgl and
  sanitize_filename; added standardized snake_case interfaces.
* Fixed reversed HR annotations, filtering after scaling/correlation, hidden
  palette dependencies, and collisions with existing binary predictor columns.
* Added explicit outcome encodings, failure reports, optional exports and
  synthetic regression tests. Unknown outcome labels now fail explicitly.
* The re-uploaded plot_roc_batch source was identical to the original upload;
  retained the existing hardened ROC implementation.

# biomed 0.1.0

* Added standardized snake_case interfaces while retaining all nine original names.
* Fixed ANOVA result binding, Venn colour validation, ROC axes and threshold ties.
* Restored ANOVA preprocessing, sparse-table Fisher tests and 0.33/0.66 grouping.
* Added XLSX export to batch analyses and PDF/PNG export to donut plots.
* Added exact six-set intersection matrix output and 0/1 survival-event validation.
* See README for compatibility notes and the new API mapping.

# biomed 0.0.0.9000

* Created the initial R package structure.
* Added column validation, automated tests, CI, and secure data-handling
  guidance.
* Added nine contributed utilities for grouping, ANOVA, categorical tests,
  Cox regression, ROC analysis, donut and stacked-bar charts, and Venn
  diagrams.
* Hardened binary-outcome handling, missing-data checks, sparse-table tests,
  multi-level Cox contrasts, and opt-in file output.
