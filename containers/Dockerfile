FROM satijalab/seurat:latest

# ── System libraries ──────────────────────────────────────────────────────────
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libuv1-dev \
        libmagick++-dev \
        libglpk-dev \
        libgmp-dev \
        libsodium-dev \
        libgit2-dev \
        cmake \
        python3 python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ── R options ─────────────────────────────────────────────────────────────────
RUN echo 'options(repos=c(CRAN="https://cloud.r-project.org/"), timeout=600, Ncpus=4)' \
    >> /usr/local/lib/R/etc/Rprofile.site

# ── Bootstrap ─────────────────────────────────────────────────────────────────
RUN R -e "install.packages(c('BiocManager','remotes'), quiet=TRUE)"

# ── readxl ────────────────────────────────────────────────────────────────────
RUN R -e "install.packages('readxl', quiet=TRUE)"

# ── SoupX (Seurat already in base image → installs cleanly) ──────────────────
RUN R -e "install.packages('SoupX', quiet=TRUE)"
RUN R -e "if (!requireNamespace('SoupX', quietly=TRUE)) stop('SoupX FAILED')"

# ── Harmony ───────────────────────────────────────────────────────────────────
RUN R -e "install.packages('harmony', quiet=TRUE)"

# ── anndata (read .h5ad from scanpy) ─────────────────────────────────────────
RUN R -e "install.packages('anndata', quiet=TRUE)"

# ── Extra plotting + utils ────────────────────────────────────────────────────
RUN R -e "install.packages(c( \
        'ggrepel','ggridges','cowplot','gridExtra', \
        'RColorBrewer','viridis','pheatmap', \
        'openxlsx','writexl','clustree', \
        'reshape2','data.table','future','future.apply' \
    ), quiet=TRUE)"

# ── ggtangle (GitHub-only dep of clusterProfiler) ────────────────────────────
RUN R -e "remotes::install_github('YuLab-SMU/ggtangle', quiet=TRUE)"

# ── Bioconductor ──────────────────────────────────────────────────────────────
RUN R -e "BiocManager::install(c( \
        'SingleCellExperiment', \
        'scran','scater','scuttle', \
        'DESeq2','edgeR','limma', \
        'enrichplot', \
        'org.Hs.eg.db','org.Mm.eg.db', \
        'GSVA','GSEABase', \
        'BiocParallel','glmGamPoi' \
    ), ask=FALSE, update=FALSE, quiet=TRUE)"

# ── clusterProfiler after ggtangle ───────────────────────────────────────────
RUN R -e "BiocManager::install('clusterProfiler', ask=FALSE, update=FALSE)"
RUN R -e "if (!requireNamespace('clusterProfiler',quietly=TRUE)) stop('clusterProfiler FAILED')"

# ── DoubletFinder (non-critical) ──────────────────────────────────────────────
RUN R -e "tryCatch( \
        remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', \
            dependencies=FALSE, quiet=TRUE), \
        error=function(e) message('WARN DoubletFinder: ', e\$message) \
    )"

# ── SeuratDisk (non-critical) ─────────────────────────────────────────────────
RUN R -e "tryCatch({ \
        install.packages('hdf5r', quiet=TRUE); \
        remotes::install_github('mojaveazure/seurat-disk', \
            dependencies=FALSE, quiet=TRUE) \
    }, error=function(e) message('WARN SeuratDisk: ', e\$message))"

# ── Final verification ────────────────────────────────────────────────────────
RUN R -e " \
    pkgs <- c('readxl','SoupX','Seurat','harmony', \
              'DESeq2','clusterProfiler','enrichplot'); \
    ok   <- sapply(pkgs, requireNamespace, quietly=TRUE); \
    cat('\n=== Package verification ===\n'); \
    for (i in seq_along(ok)) \
        cat(sprintf('  %-25s %s\n', pkgs[i], if(ok[i]) 'OK' else 'MISSING')); \
    if (!all(ok)) stop('Critical packages missing!'); \
    cat('=== All OK ===\n\n') \
"

