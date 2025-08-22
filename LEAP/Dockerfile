ARG TARGETPLATFORM=linux/amd64
ARG BUILDPLATFORM=linux/amd64

FROM --platform=$TARGETPLATFORM debian:bullseye-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

RUN apt-get update && apt-get install -y \
    build-essential gcc g++ gfortran make \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
    libjpeg-dev libpng-dev libfreetype6-dev libtiff-dev \
    libx11-dev libxt-dev \
    libglu1-mesa-dev \
    libharfbuzz-dev libfribidi-dev \
    libglpk-dev \
    curl wget unzip git \
    python3 python3-pip python3-venv \
    software-properties-common \
    gnupg ca-certificates \
    zlib1g-dev \
    libv8-dev libnode-dev \
    libudunits2-dev libgdal-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY dependency/ /workspace/dependency/

RUN PLATFORM=${TARGETPLATFORM:-linux/amd64} && \
    case $PLATFORM in \
        "linux/amd64")  MICROMAMBA_ARCH=linux-64  ;; \
        "linux/arm64")  MICROMAMBA_ARCH=linux-aarch64  ;; \
        "linux/arm/v7") MICROMAMBA_ARCH=linux-armv7l  ;; \
        *) echo "Unsupported platform: $PLATFORM, defaulting to linux-64" && MICROMAMBA_ARCH=linux-64 ;; \
    esac && \
    mkdir -p /usr/local/bin && \
    curl -Ls "https://micro.mamba.pm/api/micromamba/${MICROMAMBA_ARCH}/latest" | tar -xvj -C /tmp && \
    cp /tmp/bin/micromamba /usr/local/bin/micromamba && \
    chmod +x /usr/local/bin/micromamba && \
    rm -rf /tmp/bin /tmp/info && \
    /usr/local/bin/micromamba --version

ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV MAMBA_EXE=/usr/local/bin/micromamba
ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$MAMBA_ROOT_PREFIX/bin:$PATH
ENV R_HOME=/opt/micromamba/envs/plugin_env/lib/R

RUN mkdir -p $MAMBA_ROOT_PREFIX && \
    mkdir -p $MAMBA_ROOT_PREFIX/envs && \
    mkdir -p $MAMBA_ROOT_PREFIX/pkgs && \
    mkdir -p $MAMBA_ROOT_PREFIX/etc/profile.d

RUN /usr/local/bin/micromamba create -y -n plugin_env python=3.10 -c conda-forge --root-prefix $MAMBA_ROOT_PREFIX

RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    pip install --no-cache-dir \
    'snakemake==7.14.0' \
    'pulp==2.7.0' \
    'tabulate==0.8.10'

RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    pip install --no-cache-dir -r /workspace/dependency/requirements.txt || true

RUN /usr/local/bin/micromamba install -y -n plugin_env \
    r-base=4.4.2 \
    r-essentials \
    r-renv \
    r-jsonlite \
    -c conda-forge \
    -r $MAMBA_ROOT_PREFIX && \
    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/R /usr/local/bin/R && \
    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/Rscript /usr/local/bin/Rscript

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"), download.file.method = "libcurl")' > /root/.Rprofile && \
    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "install.packages(c('renv', 'BiocManager'), dependencies = TRUE)"

RUN if [ -f "/workspace/dependency/renv.lock" ]; then \
    cd /workspace/dependency && \
    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "\
        library(renv); \
        options(renv.config.cache.enabled = FALSE); \
        renv::init(bare = TRUE, restart = FALSE); \
        renv::restore(lockfile = '/workspace/dependency/renv.lock', confirm = FALSE);" && \
    cp -r renv/library/*/R-*/* /opt/micromamba/envs/plugin_env/lib/R/library/ 2>/dev/null || true && \
    echo 'R package installation completed'; \
    fi

RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "if (!requireNamespace('optparse', quietly = TRUE)) install.packages('optparse', dependencies = TRUE)"

ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

RUN mkdir -p /workspace/logs && \
    chmod 777 /workspace

COPY Snakefile /workspace/Snakefile
COPY visualization_Snakefile /workspace/visualization_Snakefile

COPY scripts/ /scripts/

WORKDIR /workspace

RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'export MAMBA_ROOT_PREFIX=/opt/micromamba' >> /entrypoint.sh && \
    echo 'export MAMBA_EXE=/usr/local/bin/micromamba' >> /entrypoint.sh && \
    echo 'export PATH=$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$PATH' >> /entrypoint.sh && \
    echo '# Set R environment (Micromamba 기반)' >> /entrypoint.sh && \
    echo 'export R_HOME=/opt/micromamba/envs/plugin_env/lib/R' >> /entrypoint.sh && \
    echo 'export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /entrypoint.sh && \
    echo '# Activate micromamba environment' >> /entrypoint.sh && \
    echo 'eval "$($MAMBA_EXE shell activate -s bash -p $MAMBA_ROOT_PREFIX plugin_env)" 2>/dev/null || true' >> /entrypoint.sh && \
    echo 'cd /workspace' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD test -f /opt/micromamba/envs/plugin_env/bin/python || exit 1

CMD ["/bin/bash"]