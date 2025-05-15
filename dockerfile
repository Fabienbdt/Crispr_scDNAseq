#########################################################
# Dockerfile – conteneur pour run_infercnv.R
#########################################################

FROM bioconductor/bioconductor_docker:RELEASE_3_17

# -------- Dépendances système --------
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libcurl4-openssl-dev \   # <-- on garde seulement l’une des deux variantes
        libssl-dev \
        libxml2-dev \
        build-essential \
        gfortran \
        libgfortran5 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# -------- Packages CRAN --------
RUN install2.r --error \
        optparse dplyr ggplot2 Matrix data.table Seurat

# -------- inferCNV + rhdf5 --------
RUN R -e "BiocManager::install(c('rhdf5','infercnv'), ask = FALSE)"
