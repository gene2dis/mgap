# Source of the published image docker.io/microds/medaka (digest pinned in conf/modules.config).
# The medaka 2.2.2 biocontainer, by digest, plus the model files its automatic model selection
# can choose and the bacterial model (--bacteria). The biocontainer ships only four models and
# downloads any other on first use, which fails in a pipeline task: tasks run as a non-root user
# with HOME=/ (no writable model store), on nodes that may have no network.
# A rebuilt image MUST be pushed under a new tag and the digest pin updated.
# modules/local/medaka/environment.yml lists the same medaka version (conda profile; there the
# models are downloaded into ~/.medaka on first use).
FROM --platform=linux/amd64 quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0@sha256:18d5d73aa7f0aa70d0f623c4d34ce83ca791ba9c14c751bb3ed5d3452a8fb4c3

# Consensus models of every basecaller model medaka 2.2.2 recognises (medaka.options.basecaller_models),
# and r1041_e82_400bps_bacterial_methylation; about 1.5 MB each
RUN medaka tools download_models --models \
        r103_hac_g507 \
        r1041_e82_260bps_hac_g632 \
        r1041_e82_260bps_hac_v4.0.0 \
        r1041_e82_260bps_hac_v4.1.0 \
        r1041_e82_260bps_sup_v4.0.0 \
        r1041_e82_260bps_sup_v4.1.0 \
        r1041_e82_400bps_hac_g615 \
        r1041_e82_400bps_hac_g632 \
        r1041_e82_400bps_hac_v4.0.0 \
        r1041_e82_400bps_hac_v4.1.0 \
        r1041_e82_400bps_hac_v4.2.0 \
        r1041_e82_400bps_hac_v4.3.0 \
        r1041_e82_400bps_hac_v5.0.0 \
        r1041_e82_400bps_hac_v5.2.0 \
        r1041_e82_400bps_hac_v6.0.0 \
        r1041_e82_400bps_sup_g615 \
        r1041_e82_400bps_sup_v4.0.0 \
        r1041_e82_400bps_sup_v4.1.0 \
        r1041_e82_400bps_sup_v4.2.0 \
        r1041_e82_400bps_sup_v4.3.0 \
        r1041_e82_400bps_sup_v5.0.0 \
        r1041_e82_400bps_sup_v5.2.0 \
        r941_min_fast_g507 \
        r941_min_hac_g507 \
        r941_min_sup_g507 \
        r1041_e82_400bps_bacterial_methylation \
    && ls /usr/local/lib/python3.12/site-packages/medaka/data
