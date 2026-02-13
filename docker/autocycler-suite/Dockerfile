FROM --platform=linux/amd64 mambaorg/micromamba:latest

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH

USER root

RUN micromamba install -y -n base \
    -c conda-forge -c bioconda \
    python>=3.10 \
    autocycler>=0.6.0 \
    canu>=2.3 \
    flye>=2.9.6 \
    metamdbg>=1.0 \
    miniasm>=0.3 \
    minimap2>=2.28 \
    minipolish>=0.2.0 \
    myloasm>=0.1.0 \
    necat>=0.0.1_update20200803 \
    nextdenovo>=2.5.2 \
    nextpolish>=1.4.1 \
    plassembler>=1.8.0 \
    racon>=1.5.0 \
    raven-assembler>=1.8.3 \
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
