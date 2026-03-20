#!/usr/bin/env bash
# =============================================================================
# install_bcftools_liftover.sh
#
# Builds from C/C++ sources and installs to $HOME/software:
#   1. htslib   — shared library
#   2. samtools
#   3. bcftools  (with all bundled plugins)
#   4. bcftools +liftover plugin  — from github.com/freeseek/score
#
# Every stdout/stderr line is tee'd to a timestamped master log in
# $HOME/software/logs/  and each build phase has its own sub-log.
#
# Usage:
#   bash install_bcftools_liftover.sh
#
# Optional env overrides:
#   HTSLIB_VERSION    (default: 1.21)
#   SAMTOOLS_VERSION  (default: 1.21)
#   BCFTOOLS_VERSION  (default: 1.21)
#   INSTALL_PREFIX    (default: $HOME/software)
# =============================================================================

set -euo pipefail

# ── Versions ──────────────────────────────────────────────────────────────────
HTSLIB_VERSION="${HTSLIB_VERSION:-1.21}"
SAMTOOLS_VERSION="${SAMTOOLS_VERSION:-1.21}"
BCFTOOLS_VERSION="${BCFTOOLS_VERSION:-1.21}"
# score commit/tag that contains liftover.c (no versioned releases; use main)
SCORE_REF="${SCORE_REF:-main}"

# ── Directories ───────────────────────────────────────────────────────────────
INSTALL_PREFIX="${INSTALL_PREFIX:-$HOME/software}"
SRC_DIR="${INSTALL_PREFIX}/src"
LOG_DIR="${INSTALL_PREFIX}/logs"
PLUGIN_DIR="${INSTALL_PREFIX}/libexec/bcftools"

mkdir -p "${SRC_DIR}" "${LOG_DIR}" "${PLUGIN_DIR}"

# ── Master log (tee everything) ───────────────────────────────────────────────
INSTALL_LOG="${LOG_DIR}/install_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "${INSTALL_LOG}") 2>&1

# ── Homebrew lib/include hints (macOS) ────────────────────────────────────────
_brew_prefix() { brew --prefix "$1" 2>/dev/null || echo /usr; }
BZIP2_PREFIX="$(_brew_prefix bzip2)"
ZLIB_PREFIX="$(_brew_prefix zlib)"
CURL_PREFIX="$(_brew_prefix curl)"
XZ_PREFIX="$(_brew_prefix xz)"
SSL_PREFIX="$(brew --prefix openssl@3 2>/dev/null || brew --prefix openssl 2>/dev/null || echo /usr)"

# Only add paths that actually exist (avoids 'ld: warning: search path … not found')
_add_flag() {
    local flag="$1" dir="$2"
    [[ -d "$dir" ]] && echo "$flag$dir" || true
}

CPPFLAGS=""
LDFLAGS=""
for p in "$BZIP2_PREFIX" "$ZLIB_PREFIX" "$CURL_PREFIX" "$XZ_PREFIX" "$SSL_PREFIX"; do
    [[ -d "$p/include" ]] && CPPFLAGS+=" -I${p}/include"
    [[ -d "$p/lib"     ]] && LDFLAGS+=" -L${p}/lib"
done
export CPPFLAGS LDFLAGS

PKG_CONFIG_PATH=""
for p in "$BZIP2_PREFIX" "$ZLIB_PREFIX" "$CURL_PREFIX" "$XZ_PREFIX" "$SSL_PREFIX"; do
    [[ -d "$p/lib/pkgconfig" ]] && PKG_CONFIG_PATH+="${p}/lib/pkgconfig:"
done
export PKG_CONFIG_PATH="${PKG_CONFIG_PATH%:}${PKG_CONFIG_PATH:+:}${PKG_CONFIG_PATH:-}"

# ── Helpers ───────────────────────────────────────────────────────────────────
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

log_step() {
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "  [$(timestamp)]  $1"
    echo "════════════════════════════════════════════════════════════"
}

die() { echo "ERROR: $*" >&2; exit 1; }

ncpus() { nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4; }

download_tarball() {
    local url="$1" dest="$2"
    if [[ -f "$dest" ]]; then
        echo "  Tarball already present: $dest — skipping download"
    else
        echo "  Downloading: $url"
        curl -fsSL --retry 5 --retry-delay 3 -o "$dest" "$url" \
            || die "Failed to download $url"
    fi
}

echo "Installation log : ${INSTALL_LOG}"
echo "Install prefix   : ${INSTALL_PREFIX}"
echo "htslib           : ${HTSLIB_VERSION}"
echo "samtools         : ${SAMTOOLS_VERSION}"
echo "bcftools         : ${BCFTOOLS_VERSION}"
echo "score (liftover) : ref=${SCORE_REF}"

# =============================================================================
# 1. htslib
# =============================================================================
log_step "Building htslib ${HTSLIB_VERSION}"

HTSLIB_TARBALL="${SRC_DIR}/htslib-${HTSLIB_VERSION}.tar.bz2"
HTSLIB_URL="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
HTSLIB_SRC="${SRC_DIR}/htslib-${HTSLIB_VERSION}"
HTSLIB_LOG="${LOG_DIR}/htslib-${HTSLIB_VERSION}"

download_tarball "${HTSLIB_URL}" "${HTSLIB_TARBALL}"
[[ ! -d "${HTSLIB_SRC}" ]] && tar -xjf "${HTSLIB_TARBALL}" -C "${SRC_DIR}"

(
    cd "${HTSLIB_SRC}"
    echo "  Configuring htslib..."
    ./configure \
        --prefix="${INSTALL_PREFIX}" \
        --enable-libcurl \
        --without-libdeflate \
        2>&1 | tee "${HTSLIB_LOG}.configure"

    echo "  Compiling htslib ($(ncpus) threads)..."
    make -j"$(ncpus)" 2>&1 | tee "${HTSLIB_LOG}.make"

    echo "  Installing htslib..."
    make install 2>&1 | tee "${HTSLIB_LOG}.install"
)
echo "  htslib ${HTSLIB_VERSION} installed."

# =============================================================================
# 2. samtools
# =============================================================================
log_step "Building samtools ${SAMTOOLS_VERSION}"

SAMTOOLS_TARBALL="${SRC_DIR}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
SAMTOOLS_SRC="${SRC_DIR}/samtools-${SAMTOOLS_VERSION}"
SAMTOOLS_LOG="${LOG_DIR}/samtools-${SAMTOOLS_VERSION}"

download_tarball "${SAMTOOLS_URL}" "${SAMTOOLS_TARBALL}"
[[ ! -d "${SAMTOOLS_SRC}" ]] && tar -xjf "${SAMTOOLS_TARBALL}" -C "${SRC_DIR}"

(
    cd "${SAMTOOLS_SRC}"
    echo "  Configuring samtools..."
    ./configure \
        --prefix="${INSTALL_PREFIX}" \
        --with-htslib="${INSTALL_PREFIX}" \
        2>&1 | tee "${SAMTOOLS_LOG}.configure"

    echo "  Compiling samtools..."
    make -j"$(ncpus)" 2>&1 | tee "${SAMTOOLS_LOG}.make"

    echo "  Installing samtools..."
    make install 2>&1 | tee "${SAMTOOLS_LOG}.install"
)
echo "  samtools ${SAMTOOLS_VERSION} installed."

# =============================================================================
# 3. bcftools + bundled plugins
# =============================================================================
log_step "Building bcftools ${BCFTOOLS_VERSION} (bundled plugins)"

BCFTOOLS_TARBALL="${SRC_DIR}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"
BCFTOOLS_URL="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"
BCFTOOLS_SRC="${SRC_DIR}/bcftools-${BCFTOOLS_VERSION}"
BCFTOOLS_LOG="${LOG_DIR}/bcftools-${BCFTOOLS_VERSION}"

download_tarball "${BCFTOOLS_URL}" "${BCFTOOLS_TARBALL}"
[[ ! -d "${BCFTOOLS_SRC}" ]] && tar -xjf "${BCFTOOLS_TARBALL}" -C "${SRC_DIR}"

(
    cd "${BCFTOOLS_SRC}"
    echo "  Configuring bcftools..."
    # Try with optional deps; fall back cleanly if absent
    ./configure \
        --prefix="${INSTALL_PREFIX}" \
        --with-htslib="${INSTALL_PREFIX}" \
        2>&1 | tee "${BCFTOOLS_LOG}.configure"

    echo "  Compiling bcftools..."
    make -j"$(ncpus)" 2>&1 | tee "${BCFTOOLS_LOG}.make"

    echo "  Building bundled plugins..."
    make -j"$(ncpus)" plugins 2>&1 | tee "${BCFTOOLS_LOG}.plugins" || true

    echo "  Installing bcftools..."
    make install 2>&1 | tee "${BCFTOOLS_LOG}.install"

    # Copy any .so plugins that make install missed
    if ls plugins/*.so 2>/dev/null | grep -q .; then
        echo "  Copying bundled plugin *.so → ${PLUGIN_DIR}/"
        cp plugins/*.so "${PLUGIN_DIR}/"
    fi
)
echo "  bcftools ${BCFTOOLS_VERSION} installed."

# =============================================================================
# 4. +liftover plugin  (from github.com/freeseek/score)
# =============================================================================
log_step "Building bcftools +liftover plugin (freeseek/score)"

SCORE_SRC="${SRC_DIR}/score"
LIFTOVER_LOG="${LOG_DIR}/liftover"

if [[ -d "${SCORE_SRC}/.git" ]]; then
    echo "  Updating existing score clone..."
    git -C "${SCORE_SRC}" fetch origin 2>&1 | tee "${LIFTOVER_LOG}.fetch"
    git -C "${SCORE_SRC}" checkout "${SCORE_REF}" 2>&1 | tee -a "${LIFTOVER_LOG}.fetch" || true
else
    echo "  Cloning freeseek/score (ref=${SCORE_REF})..."
    git clone --depth 1 --branch "${SCORE_REF}" \
        https://github.com/freeseek/score.git "${SCORE_SRC}" \
        2>&1 | tee "${LIFTOVER_LOG}.clone" || {
        # branch name might not exist; try cloning default and checking out SHA
        echo "  Branch '${SCORE_REF}' not found, cloning default branch..."
        git clone --depth 1 \
            https://github.com/freeseek/score.git "${SCORE_SRC}" \
            2>&1 | tee "${LIFTOVER_LOG}.clone"
    }
fi

# The liftover.c lives in the repo root
LIFTOVER_C="${SCORE_SRC}/liftover.c"
[[ -f "${LIFTOVER_C}" ]] || die "liftover.c not found in ${SCORE_SRC} — check repo structure"

echo "  Compiling liftover.so via bcftools Makefile..."
# The plugin links against internal bcftools symbols, so it must be built
# through the bcftools Makefile (not a freestanding gcc -shared call).
cp "${LIFTOVER_C}" "${BCFTOOLS_SRC}/plugins/liftover.c"
(
    cd "${BCFTOOLS_SRC}"
    make -j"$(ncpus)" plugins/liftover.so 2>&1 | tee "${LIFTOVER_LOG}.compile"
    cp plugins/liftover.so "${PLUGIN_DIR}/liftover.so"
)

echo "  liftover.so installed to ${PLUGIN_DIR}/"

# =============================================================================
# 5. PATH + BCFTOOLS_PLUGINS export
# =============================================================================
log_step "Configuring PATH and BCFTOOLS_PLUGINS"

if [[ "${SHELL:-}" == */zsh ]]; then
    RC_FILE="$HOME/.zshrc"
else
    RC_FILE="$HOME/.bashrc"
fi

MARKER="sv_imputation — bcftools/samtools"
EXPORT_BLOCK="
# >>> ${MARKER} >>>
export PATH=\"${INSTALL_PREFIX}/bin:\$PATH\"
export BCFTOOLS_PLUGINS=\"${PLUGIN_DIR}\"
# <<< ${MARKER} <<<"

if grep -qF "${MARKER}" "${RC_FILE}" 2>/dev/null; then
    echo "  PATH block already present in ${RC_FILE} — skipping."
else
    printf '%s\n' "${EXPORT_BLOCK}" >> "${RC_FILE}"
    echo "  Appended export block to ${RC_FILE}"
fi

# Apply to current session too
export PATH="${INSTALL_PREFIX}/bin:$PATH"
export BCFTOOLS_PLUGINS="${PLUGIN_DIR}"

# =============================================================================
# 6. Smoke tests
# =============================================================================
log_step "Smoke tests"

pass() { echo "  [PASS] $*"; }
fail() { echo "  [FAIL] $*"; }

# htslib shared library
if ls "${INSTALL_PREFIX}/lib/libhts"* 2>/dev/null | grep -q .; then
    pass "htslib  → $(ls ${INSTALL_PREFIX}/lib/libhts*.dylib 2>/dev/null || ls ${INSTALL_PREFIX}/lib/libhts*.so 2>/dev/null | head -1)"
else
    fail "htslib shared library not found"
fi

# samtools
if "${INSTALL_PREFIX}/bin/samtools" --version >/dev/null 2>&1; then
    pass "samtools → $("${INSTALL_PREFIX}/bin/samtools" --version | head -1)"
else
    fail "samtools binary not found"
fi

# bcftools
if "${INSTALL_PREFIX}/bin/bcftools" --version >/dev/null 2>&1; then
    pass "bcftools → $("${INSTALL_PREFIX}/bin/bcftools" --version | head -1)"
else
    fail "bcftools binary not found"
fi

# liftover plugin
echo ""
echo "  Checking +liftover plugin..."
if BCFTOOLS_PLUGINS="${PLUGIN_DIR}" \
   "${INSTALL_PREFIX}/bin/bcftools" plugin liftover --help >/dev/null 2>&1; then
    pass "+liftover → available"
else
    fail "+liftover → check ${PLUGIN_DIR}/liftover.so"
    echo "  Plugin dir contents:"
    ls "${PLUGIN_DIR}/"
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  Installation complete!"
echo "  Master log : ${INSTALL_LOG}"
echo "  Per-tool   : ${LOG_DIR}/"
echo ""
echo "  Activate in the current shell:"
echo "    source ${RC_FILE}"
echo "════════════════════════════════════════════════════════════"
