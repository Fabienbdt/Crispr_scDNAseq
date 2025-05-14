
# → image officielle Bioconductor, basée sur Ubuntu 22.04
FROM bioconductor/bioconductor_docker:RELEASE_3_17

# Paquets système manquants pour rhdf5 / fastcluster / ape…
RUN apt-get update && apt-get install -y \
        libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
        build-essential gfortran libgfortran5 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Quelques utilitaires R CRAN classiques
RUN install2.r --error optparse dplyr ggplot2 Matrix data.table

# ---- installer infercnv + dépendances Bioconductor ----
RUN R -e "BiocManager::install('infercnv', ask = FALSE)"
