ARG TARGETPLATFORM=linux/amd64
ARG BUILDPLATFORM=linux/amd64

FROM --platform=$TARGETPLATFORM debian:bullseye-slim

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

RUN apt-get update && apt-get install -y \
    build-essential gcc g++ gfortran make \
    libssl-dev libcurl4-openssl-dev libxml2-dev libxslt1-dev \
    libjpeg-dev libpng-dev libfreetype6-dev libtiff5-dev \
    libx11-dev libxt-dev xorg-dev \
    libglu1-mesa-dev \
    libharfbuzz-dev libfribidi-dev \
    libglpk-dev \
    curl wget unzip git \
    python3 python3-pip python3-venv \
    software-properties-common \
    gnupg ca-certificates \
    zlib1g-dev liblzma-dev \
    rsync \
    libv8-dev libnode-dev \
    libudunits2-dev libgdal-dev \
    libcairo2-dev \
    libfontconfig1-dev \
    libgeos-dev \
    libproj-dev \
    libmagick++-dev \
    libpoppler-cpp-dev \
    librsvg2-dev \
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
    zlib \
    xz \
    libpng \
    libjpeg-turbo \
    libcurl \
    libxml2 \
    libxslt \
    cairo \
    pango \
    harfbuzz \
    fribidi \
    freetype \
    fontconfig \
    libtiff \
    pkg-config \
    -c conda-forge \
    -r $MAMBA_ROOT_PREFIX && \
    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/R /usr/local/bin/R && \
    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/Rscript /usr/local/bin/Rscript

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"), download.file.method = "libcurl")' > /root/.Rprofile && \
    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "install.packages(c('renv', 'BiocManager'), dependencies = TRUE)"

RUN if [ -f "/workspace/dependency/renv.lock" ]; then \
    echo 'Installing R packages from renv.lock...' && \
    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "\
        # Set micromamba environment as primary library path \
        target_lib <- '/opt/micromamba/envs/plugin_env/lib/R/library'; \
        .libPaths(c(target_lib, .libPaths())); \
        cat('Library paths configured:', paste(.libPaths(), collapse=', '), '\\n'); \
        \
        # Change to dependency directory for renv operations \
        setwd('/workspace/dependency'); \
        \
        # Load renv and configure \
        library(renv); \
        options(renv.config.cache.enabled = FALSE); \
        options(renv.config.auto.snapshot = FALSE); \
        \
        # Initialize renv project \
        cat('Initializing renv project...\\n'); \
        renv::init(bare = TRUE, restart = FALSE); \
        \
        # Restore packages to renv library \
        cat('Restoring packages from renv.lock...\\n'); \
        tryCatch({ \
            renv::restore(lockfile = '/workspace/dependency/renv.lock', confirm = FALSE); \
            cat('✅ R packages successfully installed to renv library\\n'); \
        }, error = function(e) { \
            cat('⚠️ Some packages may have failed, but continuing...\\n'); \
            cat('Error details:', conditionMessage(e), '\\n'); \
        }); \
    " && \
    echo 'renv installation phase completed' ; \
fi

RUN if [ -d "/workspace/dependency/renv/library" ]; then \
    echo '=== Stage 1: Starting complete rsync copy ===' && \
    RENV_LIB_PATH="/workspace/dependency/renv/library/linux-debian-bullseye/R-4.4/x86_64-conda-linux-gnu" && \
    TARGET_LIB_PATH="/opt/micromamba/envs/plugin_env/lib/R/library" && \
    if [ -d "$RENV_LIB_PATH" ]; then \
        echo "Copy source: $RENV_LIB_PATH" && \
        echo "Copy target: $TARGET_LIB_PATH" && \
        echo "Packages to copy:" && \
        ls -1 "$RENV_LIB_PATH" | wc -l && echo " packages" && \
        echo "Starting rsync copy..." && \
        rsync -av --progress "$RENV_LIB_PATH/" "$TARGET_LIB_PATH/" && \
        echo "✅ rsync copy completed" && \
        echo "Post-copy package verification:" && \
        ls -1 "$TARGET_LIB_PATH" | wc -l && echo " packages" ; \
    else \
        echo "⚠️ Cannot find renv library path: $RENV_LIB_PATH" ; \
    fi ; \
else \
    echo "⚠️ renv library directory does not exist" ;fi

RUN echo '=== Stage 2: Cleanup temporary files for space optimization ===' && \
    echo '=== Stage 2: Starting temporary file cleanup ===' && \
    # Remove renv temporary library directory \
    if [ -d "/workspace/dependency/renv/library" ]; then \
        echo 'Removing renv temporary library...' && \
        RENV_LIB_SIZE=$(du -sh /workspace/dependency/renv/library 2>/dev/null | cut -f1 || echo "unknown") && \
        echo "renv library size to remove: $RENV_LIB_SIZE" && \
        rm -rf /workspace/dependency/renv/library && \
        echo '✅ renv temporary library removal completed' ; \
    fi && \
    # Remove renv cache \
    if [ -d "/root/.cache/R/renv" ]; then \
        echo 'Removing renv cache...' && \
        CACHE_SIZE=$(du -sh /root/.cache/R/renv 2>/dev/null | cut -f1 || echo "unknown") && \
        echo "Cache size to remove: $CACHE_SIZE" && \
        rm -rf /root/.cache/R/renv && \
        echo '✅ renv cache removal completed' ; \
    fi && \
    # Remove rsync package \
    echo 'Removing rsync package...' && \
    apt-get remove -y rsync && \
    apt-get autoremove -y && \
    echo '✅ rsync package removal completed' && \
    echo '=== Stage 2: Temporary file cleanup completed ==='

RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "if (!requireNamespace('optparse', quietly = TRUE)) install.packages('optparse', dependencies = TRUE)"

ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

RUN mkdir -p /workspace/logs && \
    chmod 777 /workspace

COPY Snakefile /workspace/Snakefile

COPY scripts/ /scripts/

WORKDIR /workspace

RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'export MAMBA_ROOT_PREFIX=/opt/micromamba' >> /entrypoint.sh && \
    echo 'export MAMBA_EXE=/usr/local/bin/micromamba' >> /entrypoint.sh && \
    echo 'export PATH=$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$PATH' >> /entrypoint.sh && \
    echo 'export R_HOME=/opt/micromamba/envs/plugin_env/lib/R' >> /entrypoint.sh && \
    echo 'export R_LIBS_USER=/opt/micromamba/envs/plugin_env/lib/R/library' >> /entrypoint.sh && \
    echo 'export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /entrypoint.sh && \
    echo 'eval "$($MAMBA_EXE shell activate -s bash -p $MAMBA_ROOT_PREFIX plugin_env)" 2>/dev/null || true' >> /entrypoint.sh && \
    echo 'cd /workspace' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD test -f /opt/micromamba/envs/plugin_env/bin/python || exit 1

CMD ["/bin/bash"]