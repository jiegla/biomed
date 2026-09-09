# Security and data handling

## Clinical and research data

- Do not commit protected health information, direct identifiers, linkage
  tables, raw exports from clinical systems, or unreviewed patient-level data.
- Keep sensitive inputs outside the repository. Pass their locations to R at
  runtime through function arguments or local environment variables.
- Use synthetic or formally de-identified data for examples and automated
  tests. A study identifier alone does not guarantee de-identification.
- Review plots, tables, logs, error messages, caches, and serialized R objects
  before sharing; they can contain patient-level values or filesystem paths.

## Credentials and tokens

- Local development: keep secrets in the user-level `.Renviron`, a credential
  manager, or another ignored local secret store. Never commit `.Renviron`.
- R code: read credentials only at the point of use with `Sys.getenv()` and
  fail clearly when they are absent. Never provide a real token as a default.
- GitHub Actions: use repository or environment secrets and reference them as
  `${{ secrets.NAME }}`. Grant the workflow token the minimum permissions
  required for that job.
- Requests: send tokens only to the intended HTTPS origin, normally in an
  `Authorization` header. Never add tokens to query strings or print request
  headers in logs.
- Rotation: revoke and replace a credential immediately if it appears in a
  commit, issue, workflow log, artifact, or release. Removing it from the latest
  commit does not remove it from Git history.

## GitHub Actions token model

The included R-CMD-check workflow uses GitHub's short-lived
`secrets.GITHUB_TOKEN`. GitHub creates it for the workflow run and invalidates
it after the job finishes. The workflow declares `contents: read`, so package
checks can clone the repository without write access.

## Reporting a problem

For a suspected credential leak or exposure of clinical data, do not open a
public issue containing the material. Revoke access first, preserve only the
minimum evidence needed for investigation, and contact the repository owner
privately.
