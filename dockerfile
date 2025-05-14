FROM rocker/r-ver:4.3.2

# Dépendances système (pour rhdf5, infercnv, etc.)
RUN apt-get update && apt-get install -y \
    libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev libxml2-utils

# Installer les packages R nécessaires
RUN install2.r --error \
    optparse dplyr ggplot2 Matrix data.table BiocManager

# Installer infercnv + bioconductor
RUN R -e "BiocManager::install(version = '3.17', ask = FALSE); \
          BiocManager::install(c('rhdf5','infercnv','IRanges','edgeR','locfit'))"
