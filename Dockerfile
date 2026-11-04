# nf-wgs-genotyping runtime image (bin/ scripts are supplied by Nextflow at run time).
#   docker build -t <registry>/nf-wgs-genotyping:1.0.0 .
#   docker push <registry>/nf-wgs-genotyping:1.0.0
FROM mambaorg/micromamba:1.5.10-jammy

LABEL org.opencontainers.image.title="nf-wgs-genotyping" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.description="Bulk WGS SNV/SV re-genotyping (pysam, svtyper fork 27a0427)"

# Nextflow runs `ps` inside the container to record task CPU/memory metrics (trace/report);
# the micromamba base image does not ship it. Installed from apt (Linux-only) rather than
# environment.yml so the conda env still solves on macOS (osx-arm64 has no procps-ng build).
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes \
    && /opt/conda/bin/python -c "import pysam, pandas, yaml, svtyper.classic"

ENV PATH=/opt/conda/bin:$PATH
