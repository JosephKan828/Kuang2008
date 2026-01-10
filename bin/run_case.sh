#!/bin/sh

set -euo pipefail

# ==============================================
# Purpose:
# 	- Run Julia experiment with a selected CASE
# 	- Ensure `output/<CASE>/` exists
# 	- Pass `CASE` in a consistent way
# ==============================================

CASE_RAW="${1:-no_rad}"
RAD_SCALE="${2:-0.001}"

# folder-safe normalization
case "${CASE_RAW}" in
  "qt+cld_rad") CASE="qt_cld_rad" ;;
  *)           CASE="${CASE_RAW}" ;;
esac

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${ROOT}/output/${CASE}"

mkdir -p "${OUTDIR}"

JULIA_SCRIPT="${ROOT}/experiments/run_case.jl"
if [[ ! -f "${JULIA_SCRIPT}" ]]; then
  echo "[ERROR] Missing ${JULIA_SCRIPT}"
  echo "Create a single driver: experiments/run_case.jl"
  exit 1
fi

julia --project=. "${JULIA_SCRIPT}" "${CASE}" "${RAD_SCALE}" "${OUTDIR}" 
