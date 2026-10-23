# Source of the published image docker.io/microds/autocycler (tag and digest
# pinned in conf/modules.config). Every package is pinned to an exact version and
# the base image to a tag and digest, so a rebuild resolves the same versions; a
# rebuilt image MUST still be pushed under a new tag and the digest pin updated.
# modules/local/autocycler/environment.yml lists the same versions (conda profile).
# dnaapler is a Plassembler dependency, pinned because the solver otherwise picks 0.8.1.
FROM --platform=linux/amd64 mambaorg/micromamba:2.8.1-debian12-slim@sha256:c8198d53228ad7cfd7adcc0704e8837f9d1c9327fb363c6a7d3c5b1a51a4b561

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

USER root

RUN micromamba install -y -n base \
    -c conda-forge -c bioconda \
    "autocycler=0.8.0" \
    "canu=2.3" \
    "dnaapler=1.4.0" \
    "flye=2.9.6" \
    "metamdbg=1.4" \
    "miniasm=0.3" \
    "minimap2=2.31" \
    "minipolish=0.2.1" \
    "myloasm=0.7.0" \
    "necat=0.0.1_update20200803" \
    "nextdenovo=2.5.2" \
    "nextpolish=1.4.1" \
    "plassembler=1.8.5" \
    "racon=1.5.0" \
    "raven-assembler=1.8.3" \
    && micromamba clean --all --yes

RUN apt-get update && apt-get install -y --no-install-recommends \
    gawk \
    grep \
    procps \
    sed \
    coreutils \
    && rm -rf /var/lib/apt/lists/*

SHELL ["micromamba", "run", "-n", "base", "/bin/bash", "-lc"]

CMD ["micromamba", "run", "-n", "base", "autocycler", "--help"]
