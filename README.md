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

The package is under active development. After the repository is available:

```r
pak::pak("jiegla/biomed")
```

## First utility

```r
library(biomed)

cohort <- data.frame(
  patient_id = c("P001", "P002"),
  response = c("PR", "SD")
)

biomed_check_columns(cohort, c("patient_id", "response"))
```

## Proposed modules

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
