#!/bin/sh

set -euo pipefail

# ============================================
# Purpose of this file
# 	- determine repo root
# 	- choose case
# 	- active python environment
# 	- call make_all.py with standardized paths
# 	- create output directories of missing
# =============================================

CASE_RAW="${1:-no_rad}"
RAD_SCALE="${2:-0.001}"

case "${CASE_RAW}" in
  "qt+cld_rad") CASE="qt_cld_rad" ;;
  *)           CASE="${CASE_RAW}" ;;
esac

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="/work/b11209013/Kuang2008"

if [[ "${CASE}" == "no_rad" ]] ; then
  IN_DIR="${WORK}/output/${CASE}"
  FIG_DIR="${ROOT}/figures/${CASE}"
  POST_DIR="${WORK}/output/${CASE}/post"

else
  IN_DIR="${WORK}/output/${CASE}/rad_scaling=${RAD_SCALE}"
  FIG_DIR="${ROOT}/figures/${CASE}/rad_scaling=${RAD_SCALE}"
  POST_DIR="${WORK}/output/${CASE}/post/rad_scaling=${RAD_SCALE}"
fi

mkdir -p "${FIG_DIR}" "${POST_DIR}"

VENV_ACT="${ROOT}/post/.venv/bin/activate"
if [[ ! -f "${VENV_ACT}" ]]; then
  echo "[ERROR] Python venv not found: ${VENV_ACT}"
  echo "Create it from repo root:"
  echo "  python -m venv post/.venv"
  echo "  source post/.venv/bin/activate"
  echo "  pip install -U pip"
  echo "  pip install -e post"
  exit 1
fi
# shellcheck disable=SC1090
source "${VENV_ACT}"

python -m kuang_post.make_all \
  --case "${CASE}" \
  --in-dir "${IN_DIR}" \
  --fig-dir "${FIG_DIR}" \
  --post-dir "${POST_DIR}"\
