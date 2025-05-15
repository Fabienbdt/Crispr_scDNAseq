FROM bioconductor/bioconductor_docker:RELEASE_3_17

# -------- Dépendances système --------
RUN apt-get update && apt-get install -y \
    libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    build-essential gfortran libgfortran5 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# -------- Packages CRAN utiles --------
RUN install2.r --error \
    optparse dplyr ggplot2 Matrix data.table

# -------- Installer Seurat depuis CRAN --------
RUN R -e "install.packages('Seurat', repos = 'https://cloud.r-project.org')"

# -------- infercnv + rhdf5 + dépendances Bioconductor --------
RUN R -e "BiocManager::install(c('rhdf5', 'infercnv'), ask = FALSE)"
