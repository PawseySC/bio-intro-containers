FROM rocker/tidyverse:3.6.1

# Optional: document with labels
LABEL maintainer="Pawsey"

# Install required Bioconductor package
RUN R -e 'BiocManager::install("ggtree")'
