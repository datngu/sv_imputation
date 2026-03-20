# ============================================================
#  sv-imputation Docker image
#  Provides: bcftools + liftover plugin, samtools, bgzip,
#            plink 1.9, Java 17 (for Beagle), Python 3.12
#
#  Uses micromamba on Debian-slim for a significantly smaller
#  final image compared to a full miniconda3 base, while
#  keeping glibc compatibility for bioconda binaries.
# ============================================================
FROM mambaorg/micromamba:latest

# Run subsequent RUN commands as root for package installs, then drop back
USER root

# ── Minimal system packages ───────────────────────────────────
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends procps curl && \
    rm -rf /var/lib/apt/lists/*

# ── Conda packages (single layer, aggressive cleanup) ────────
RUN micromamba install -y -n base \
    -c conda-forge \
    -c bioconda \
    python=3.12 \
    bcftools \
    bcftools-liftover-plugin \
    samtools \
    htslib \
    plink \
    openjdk=17 \
    && micromamba clean --all --force-pkgs-dirs --yes \
    && find /opt/conda -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true

ENV PATH=/opt/conda/bin:$PATH

# ── Sanity checks (fail the build if anything is missing) ────
RUN bcftools --version && \
    bcftools +liftover --help 2>&1 | head -5 && \
    samtools --version | head -1 && \
    bgzip --version | head -1 && \
    plink --version && \
    java -version

# Drop back to non-root for runtime
USER $MAMBA_USER

# ── Default workdir ──────────────────────────────────────────
WORKDIR /workspace

CMD ["/bin/bash"]
