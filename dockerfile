#########################################################
# Dockerfile – conteneur pour run_infercnv.R
#########################################################

FROM bioconductor/bioconductor_docker:RELEASE_3_17

# -------- Dépendances système --------
RUN apt-get update && apt-get install -y \
    libhdf5-dev libcurl4-openssl-dev libssl-dev libcurl4-gnutls-dev \
    libxml2-dev build-essential gfortran libgfortran5 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# -------- Packages CRAN utiles --------
RUN install2.r --error \
    optparse dplyr ggplot2 Matrix data.table

# -------- Seurat 4 (compatibilité infercnv) --------
RUN R -e "install.packages('Seurat', version = '4.3.0', repos = 'https://cloud.r-project.org')"

# -------- Bioconductor : infercnv + rhdf5 --------
RUN R -e "BiocManager::install(c('rhdf5','infercnv'), ask = FALSE)"
