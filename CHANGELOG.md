# Changelog

## v1.0.0 (2026-11-04)

- Two entry points: `snv_genotyping.nf` (PASS SNV union, fragment-level pysam counting) and
  `sv_genotyping.nf` (BRASS/GRIDSS/SvABA PASS SV union, svtyper fork @27a0427).
- GRIDSS somatic filter on matched-normal VF; BRASS missing FILTER treated as PASS.
- Explicit, manifest-checked scatter/gather; csverve-style output YAML; QC JSON with a
  positive control (target's own calls).
- Profiles: iris, slurm, singularity, docker, conda, test, test_snv, test_sv, stub.
- `--python`: interpreter for `bin/` scripts (default `python3`); `procps` comes from apt in the
  Dockerfile so `environment.yml` also solves on macOS.
- Container image: record the pushed digest here (`<registry>/nf-wgs-genotyping:1.0.0@sha256:…`).
