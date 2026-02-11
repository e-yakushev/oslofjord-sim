#!/usr/bin/env bash
set -euo pipefail

############################################################
# Job ID + local output dirs
############################################################

# Batch sets BATCH_JOB_ID; locally we default to "local"
JOB_ID="${BATCH_JOB_ID:-local}"

LOCAL_INPUT_DIR="/app/data/input"
LOCAL_OUTPUT_ROOT="/app/data/output"
LOCAL_OUTPUT_DIR="${LOCAL_OUTPUT_ROOT}/${JOB_ID}"

mkdir -p "${LOCAL_INPUT_DIR}" "${LOCAL_OUTPUT_DIR}"

############################################################
# Logging: tee everything to a log file + console
############################################################

LOG_FILE="${LOCAL_OUTPUT_DIR}/simulation.log"
echo "📄 Logging to ${LOG_FILE}"
exec > >(tee -a "${LOG_FILE}") 2>&1

############################################################
# Basic configuration
############################################################

GCS_INPUT_PREFIX="${GCS_INPUT_PREFIX:-input}"
GCS_OUTPUT_PREFIX="${GCS_OUTPUT_PREFIX:-output}"

# Which Julia script to run (can be overridden via env)
SCRIPT_PATH="${SIMULATION_LAUNCHER:-/app/main.jl}"

echo "=== FjordSim / Oceananigans entrypoint ==="
echo "  JOB_ID          = ${JOB_ID}"
echo "  SCRIPT_PATH     = ${SCRIPT_PATH}"
echo "  LOCAL_INPUT_DIR = ${LOCAL_INPUT_DIR}"
echo "  LOCAL_OUTPUT_DIR= ${LOCAL_OUTPUT_DIR}"

############################################################
# Mode selection: local vs cloud
############################################################

if [ -z "${GCS_BUCKET:-}" ]; then
  echo "🌍 Local development mode (no GCS_BUCKET set)"
  echo "  Input  → ${LOCAL_INPUT_DIR}"
  echo "  Output → ${LOCAL_OUTPUT_DIR}"
else
  echo "☁️  Cloud mode (GCS_BUCKET=${GCS_BUCKET})"
  echo "  Input prefix:  gs://${GCS_BUCKET}/${GCS_INPUT_PREFIX}/"
  echo "  Output prefix: gs://${GCS_BUCKET}/${GCS_OUTPUT_PREFIX}/${JOB_ID}/"

  echo "=== Downloading input from GCS → ${LOCAL_INPUT_DIR}"
  gsutil -m cp -r "gs://${GCS_BUCKET}/${GCS_INPUT_PREFIX}/*" "${LOCAL_INPUT_DIR}" || \
    echo "⚠️ No input files found at gs://${GCS_BUCKET}/${GCS_INPUT_PREFIX}/ (continuing)."
fi

############################################################
# Verify input contents (if any)
############################################################

echo "=== Verifying input files in ${LOCAL_INPUT_DIR}"
if [ -d "${LOCAL_INPUT_DIR}" ] && [ "$(ls -A "${LOCAL_INPUT_DIR}" 2>/dev/null)" ]; then
  echo "📂 Contents of ${LOCAL_INPUT_DIR}:"
  find "${LOCAL_INPUT_DIR}" -maxdepth 2 -type f -exec ls -lh {} \;
  echo "✅ Input directory is ready."
else
  echo "⚠️ ${LOCAL_INPUT_DIR} is empty or missing. Check your local input mount or GCS paths."
fi


############################################################
# Run Julia simulation
############################################################

echo "=== Running Julia simulation (job: ${JOB_ID}) ==="
echo "Script: ${SCRIPT_PATH}"

julia --color=no --project=/app "${SCRIPT_PATH}" "$@"

status=$?

if [ "${status}" -ne 0 ]; then
  echo "❌ Simulation exited with status ${status}"
else
  echo "✅ Simulation completed successfully."
fi

############################################################
# Upload output in cloud mode
############################################################

if [ -n "${GCS_BUCKET:-}" ]; then
  echo "=== Uploading output from ${LOCAL_OUTPUT_DIR} → gs://${GCS_BUCKET}/${GCS_OUTPUT_PREFIX}/${JOB_ID}/"
  # The wildcard must be outside quotes to expand, but we still quote the prefix.
  gsutil -m cp -r "${LOCAL_OUTPUT_DIR}/"* "gs://${GCS_BUCKET}/${GCS_OUTPUT_PREFIX}/${JOB_ID}/" || \
    echo "⚠️ No output files to upload."
else
  echo "✅ Local run complete. Output saved to: ${LOCAL_OUTPUT_DIR}"
fi

echo "=== Done (job: ${JOB_ID}) ==="
exit "${status}"
