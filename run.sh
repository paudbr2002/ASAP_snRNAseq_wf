#!/bin/bash
# ASAP snRNA-seq Pipeline — COSMOS-SENS
#
# Usage:
#   bash run.sh                           # SoupX + Harmony, all samples
#   bash run.sh -resume                   # resume after failure
#   bash run.sh --skip-ambient            # no ambient correction
#   bash run.sh --scvi                    # scVI integration (GPU)
#   bash run.sh --males                   # integrate males only
#   bash run.sh --females                 # integrate females only
#   bash run.sh --scvi --males            # scVI + males only
#
# Extra metadata filter (edit EXTRA_FILTER below):
#   EXTRA_FILTER="region:SN,putamen"
#   EXTRA_FILTER="condition:PD;region:SN"
#   EXTRA_FILTER=""   # no filter (default)

set -e

PIPELINE_DIR=/home/paudbr/projects/mLOY_PD/scripts/nf_scpipeline
WORKDIR=/home/paudbr/projects/mLOY_PD
cd "$WORKDIR"

# ── Defaults ──────────────────────────────────────────────────────────────────
AMBIENT_TOOL="soupx"
INTEGRATION_TOOL="harmony"
SEX_FILTER="all"
EXTRA_FILTER=""
RESUME=""

# ── Parse flags ───────────────────────────────────────────────────────────────
for arg in "$@"; do
    case "$arg" in
        --skip-ambient)  AMBIENT_TOOL="skip"         ;;
        --scvi)          INTEGRATION_TOOL="scvi"     ;;
        --males)         SEX_FILTER="M"              ;;
        --females)       SEX_FILTER="F"              ;;
        -resume)         RESUME="-resume"            ;;
    esac
done

echo "ambient_tool     : ${AMBIENT_TOOL}"
echo "integration_tool : ${INTEGRATION_TOOL}"
echo "sex_filter       : ${SEX_FILTER}"
echo "extra_filter     : ${EXTRA_FILTER:-'(none)'}"
echo ""

module load Nextflow/24.10.0

nextflow run "${PIPELINE_DIR}/main.nf" \
    --samplesheet      "${WORKDIR}/data/ASAP_samplesheet_with_paths.xlsx" \
    --outdir           "${WORKDIR}/results_ASAP_jakobsson_cellbender_scivi" \
    --ambient_tool     "${AMBIENT_TOOL}" \
    --integration_tool "${INTEGRATION_TOOL}" \
    --sex_filter       "${SEX_FILTER}" \
    --extra_filter     "${EXTRA_FILTER}" \
    -c /scale/gr01/shared/nf-core/cosmos_sens_grp1_NFv24.config \
    -c "${PIPELINE_DIR}/nextflow.config" \
    -w "${WORKDIR}/work" \
    ${RESUME}
