#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ══════════════════════════════════════════════════════════════════════════════
//  ASAP Parkinson's snRNA-seq Pipeline
//
//  Step 1 — PREPARE_SAMPLES    : filter samplesheet, derive paths
//  Step 2 — ambient correction : SOUPX | skip
//  Step 3 — QC_PER_SAMPLE      : scrublet, MT%, adaptive gene filter
//  Step 4 — VALIDATE_SEX       : XIST / chrY heatmap + exclusion
//  Step 5 — integration        : INTEGRATE_HARMONY | INTEGRATE_SCVI
//
//  Key parameters
//    --ambient_tool    soupx (default) | skip
//    --integration_tool harmony (default) | scvi
//    --sex_filter      all (default) | M | F
//    --extra_filter    "col1:val1,val2;col2:val3"  (empty = none)
// ══════════════════════════════════════════════════════════════════════════════

log.info """
╔══════════════════════════════════════════════════════════════════╗
║   ASAP snRNA-seq Pipeline  —  mLOY / Parkinson's               ║
╠══════════════════════════════════════════════════════════════════╣
║  samplesheet      : ${params.samplesheet}
║  outdir           : ${params.outdir}
║  ambient_tool     : ${params.ambient_tool}
║  integration_tool : ${params.integration_tool}
║  sex_filter       : ${params.sex_filter}
║  extra_filter     : ${params.extra_filter ?: '(none)'}
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 1 — Samplesheet validation + path derivation
//  Applies hypoxia filter + any extra_filter.
//  Derives raw_matrix_path (for CellBender) from matrix_path.
// ══════════════════════════════════════════════════════════════════════════════
process PREPARE_SAMPLES {
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path samplesheet

    output:
    path "samples_valid.csv",   emit: samples_csv
    path "samples_missing.txt", emit: missing

    script:
    def ef = params.extra_filter ? "--extra_filter '${params.extra_filter}'" : ""
    """
    prepare_samples.py \\
        --samplesheet ${samplesheet} \\
        --out_csv     samples_valid.csv \\
        --out_missing samples_missing.txt \\
        ${ef}
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 2a — SoupX (R, CPU)
//  Reads outs/ via load10X(), estimates contamination with autoEstCont().
//  Output: per-sample MTX directory.
// ══════════════════════════════════════════════════════════════════════════════
process SOUPX {
    tag "${sample_id}"
    maxForks params.max_parallel
    publishDir "${params.outdir}/ambient_reports", mode: 'copy', pattern: "*_report.csv"

    input:
    tuple val(sample_id), val(meta)

    output:
    tuple val(sample_id), val(meta), path("${sample_id}_soupx"), emit: corrected
    path "${sample_id}_soupx_report.csv",                        emit: report

    script:
    """
    soupx.R \\
        --sample_id   ${sample_id} \\
        --matrix_path ${meta.matrix_path} \\
        --outdir      ${sample_id}_soupx \\
        --report      ${sample_id}_soupx_report.csv
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 3 — Per-sample QC (scanpy + scrublet)
//  --corrected accepts: MTX directory (SoupX), .h5 file (CellBender), or
//  the filtered_feature_bc_matrix dir itself (skip mode).
//  Detection is automatic inside qc_per_sample.py.
// ══════════════════════════════════════════════════════════════════════════════
process QC_PER_SAMPLE {
    tag "${sample_id}"
    maxForks params.max_parallel
    publishDir "${params.outdir}/per_sample_qc", mode: 'copy', pattern: "*.h5ad"
    publishDir "${params.outdir}/per_sample_qc", mode: 'copy', pattern: "*.EXCLUDED"

    input:
    tuple val(sample_id), val(meta), path(corrected_output)

    output:
    tuple val(sample_id), path("${sample_id}.h5ad"), optional: true, emit: h5ad
    path "*.EXCLUDED",                               optional: true, emit: excluded

    script:
    def meta_json = groovy.json.JsonOutput.toJson(meta)
    """
    cat > meta.json << 'METAEOF'
${meta_json}
METAEOF

    qc_per_sample.py \\
        --sample_id  ${sample_id} \\
        --corrected  ${corrected_output} \\
        --filtered   ${meta.matrix_path} \\
        --meta_json  meta.json
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 4 — Sex annotation validation
// ══════════════════════════════════════════════════════════════════════════════
process VALIDATE_SEX {
    publishDir "${params.outdir}/sex_validation", mode: 'copy'

    input:
    path h5ad_files
    path samplesheet

    output:
    path "validated_samples.csv",              emit: valid_csv
    path "excluded_samples.csv",               emit: excluded_csv
    path "sex_corrections.csv",                emit: corrections
    path "ASAP_samplesheet_sexcorrected.xlsx", emit: corrected_xlsx, optional: true
    path "heatmap_chrY_XIST_per_patient.png",  emit: heatmap

    script:
    """
    validate_sex.py \\
        --h5ad_dir    . \\
        --samplesheet ${samplesheet} \\
        --outdir      .
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 5a — Harmony integration (CPU)
// ══════════════════════════════════════════════════════════════════════════════
process INTEGRATE_HARMONY {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path h5ad_files
    val  sex_filter

    output:
    path "ASAP_integrated.h5ad", emit: integrated
    path "plots/*",              emit: plots
    path "tables/*",             emit: tables

    script:
    """
    integrate.py \\
        --h5ad_dir   . \\
        --outdir     . \\
        --sex_filter ${sex_filter}
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  PROCESS 6b — scVI-tools integration (GPU)
// ══════════════════════════════════════════════════════════════════════════════
process INTEGRATE_SCVI {
    publishDir "${params.outdir}/integration", mode: 'copy'

    input:
    path h5ad_files
    val  sex_filter

    output:
    path "ASAP_integrated.h5ad", emit: integrated
    path "scvi_model/",          emit: model
    path "plots/*",              emit: plots
    path "tables/*",             emit: tables

    script:
    """
    integrate_scvi.py \\
        --h5ad_dir   . \\
        --outdir     . \\
        --sex_filter ${sex_filter}
    """
}


// ══════════════════════════════════════════════════════════════════════════════
//  WORKFLOW
// ══════════════════════════════════════════════════════════════════════════════
workflow {

    // ── 1. Parse samplesheet ──────────────────────────────────────────────────
    PREPARE_SAMPLES(file(params.samplesheet))

    samples_ch = PREPARE_SAMPLES.out.samples_csv
        .splitCsv(header: true)
        .map { row ->
            def meta = row.collectEntries { k, v -> [k, v ?: ""] }
            tuple(row.sample_id, meta)
        }

    // ── 2. Ambient RNA correction — conditional on ambient_tool ───────────────
    def ambient_ch

    if (params.ambient_tool == "soupx") {
        SOUPX(samples_ch)
        ambient_ch = SOUPX.out.corrected
    } else {
        // "skip" — pass the filtered matrix directly; QC reads it as-is
        ambient_ch = samples_ch.map { sid, meta ->
            tuple(sid, meta, file(meta.matrix_path))
        }
    }

    // ── 3. QC per sample ──────────────────────────────────────────────────────
    QC_PER_SAMPLE(ambient_ch)

    // ── 4. Sex validation ─────────────────────────────────────────────────────
    h5ads_all = QC_PER_SAMPLE.out.h5ad
        .map { sid, h5ad -> h5ad }
        .collect()

    VALIDATE_SEX(h5ads_all, file(params.samplesheet))

    // ── 5. Filter h5ads to sex-validated samples ─────────────────────────────
    verified_h5ads = QC_PER_SAMPLE.out.h5ad
        .combine(VALIDATE_SEX.out.valid_csv)
        .filter { sid, h5ad, csv ->
            def ids = csv.text.split('\n').drop(1).collect { it.split(',')[0].trim() }
            sid in ids
        }
        .map { sid, h5ad, _csv -> h5ad }
        .collect()

    // ── 6. Integration ────────────────────────────────────────────────────────
    if (params.integration_tool == "scvi") {
        INTEGRATE_SCVI(verified_h5ads, params.sex_filter)
    } else {
        INTEGRATE_HARMONY(verified_h5ads, params.sex_filter)
    }
}


workflow.onComplete {
    log.info """
    Pipeline ${workflow.success ? 'COMPLETED ✓' : 'FAILED ✗'}
      ambient    : ${params.ambient_tool}
      integrate  : ${params.integration_tool}
      sex        : ${params.sex_filter}
      duration   : ${workflow.duration}
      output     : ${params.outdir}
    """.stripIndent()
}
