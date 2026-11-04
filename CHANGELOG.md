# Changelog

## v1.0.1 (2026-11-04)

- `params.container = 'mattmyers/nf-wgs-genotyping:1.0.1'` (Docker Hub; Nextflow adds `docker://`
  for singularity). Record the pushed digest here: `mattmyers/nf-wgs-genotyping:1.0.1@sha256:…`.
- `--python`: interpreter for `bin/` scripts (default `python3`); `procps` comes from apt in the
  Dockerfile so `environment.yml` also solves on macOS.
- Test: `params.container` must be a pullable `[registry/]<namespace>/<image>:<tag>` reference.

## v1.0.0 (2026-11-04)

- Two entry points: `snv_genotyping.nf` (PASS SNV union, fragment-level pysam counting) and
  `sv_genotyping.nf` (BRASS/GRIDSS/SvABA PASS SV union, svtyper fork @27a0427).
- GRIDSS somatic filter on matched-normal VF; BRASS missing FILTER treated as PASS.
- Explicit, manifest-checked scatter/gather; csverve-style output YAML; QC JSON with a
  positive control (target's own calls).
- Profiles: iris, slurm, singularity, docker, conda, test, test_snv, test_sv, stub.
- `params.container` was a placeholder; this tag cannot pull an image on iris.
