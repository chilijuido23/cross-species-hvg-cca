# =============================================================================
# Cross-Species HVG-CCA: Docker Container
# =============================================================================
# Build:  docker build -t cross-species-hvg-cca .
# Run:    docker run -v $(pwd)/data:/data -v $(pwd)/results:/results \
#             cross-species-hvg-cca snakemake --cores 8
# =============================================================================

FROM rocker/r-ver:4.3.0

LABEL maintainer="chilijuido23"
LABEL description="Cross-species scRNA-seq integration with adaptive ortholog alignment"

# System dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgsl-dev \
    cmake \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN install2.r --error --skipinstalled \
    Seurat \
    ggplot2 \
    patchwork \
    optparse \
    yaml \
    Matrix \
    rmarkdown \
    knitr \
    gridExtra \
    RColorBrewer \
    scales

# Install Snakemake
RUN pip3 install snakemake>=7.30.0

# Copy project
WORKDIR /app
COPY . .

# Default command
CMD ["snakemake", "--cores", "8"]
