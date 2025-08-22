# 멀티 플랫폼 지원
# BuildKit 자동 감지 ARG (fallback 기본값)
ARG TARGETPLATFORM=linux/amd64
ARG BUILDPLATFORM=linux/amd64

FROM --platform=$TARGETPLATFORM debian:bullseye-slim

# 비대화형 설치 설정
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

# 기본 패키지 설치 (V8, jsonvalidate 지원 패키지 포함)
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

# dependency 폴더 복사
COPY dependency/ /workspace/dependency/

# Micromamba 설치 (멀티 플랫폼, 기본값: linux/amd64)
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

# 환경 변수 설정 (Micromamba 기반으로 통일)
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV MAMBA_EXE=/usr/local/bin/micromamba
ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$MAMBA_ROOT_PREFIX/bin:$PATH
ENV R_HOME=/opt/micromamba/envs/plugin_env/lib/R

# Micromamba 환경 디렉토리 생성
RUN mkdir -p $MAMBA_ROOT_PREFIX && \
    mkdir -p $MAMBA_ROOT_PREFIX/envs && \
    mkdir -p $MAMBA_ROOT_PREFIX/pkgs && \
    mkdir -p $MAMBA_ROOT_PREFIX/etc/profile.d

# Python 환경 생성
RUN /usr/local/bin/micromamba create -y -n plugin_env python=3.10 -c conda-forge --root-prefix $MAMBA_ROOT_PREFIX

# Snakemake 및 필수 패키지 설치
RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    pip install --no-cache-dir \
    'snakemake==7.14.0' \
    'pulp==2.7.0' \
    'tabulate==0.8.10'

# Python 패키지 설치
RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    pip install --no-cache-dir -r /workspace/dependency/requirements.txt || true

# R 4.4.2 및 필수 개발 라이브러리 설치 (Micromamba 사용)
RUN /usr/local/bin/micromamba install -y -n plugin_env \
    r-base=4.4.2 \
    r-essentials \
    r-renv \
    r-jsonlite \
    -c conda-forge \
    -r $MAMBA_ROOT_PREFIX && \
    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/R /usr/local/bin/R && \
    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/Rscript /usr/local/bin/Rscript

# R 설정 및 초기 패키지 설치
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"), download.file.method = "libcurl")' > /root/.Rprofile && \
    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "install.packages(c('renv', 'BiocManager'), dependencies = TRUE)"

# renv를 사용한 패키지 설치
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

# 추가 R 패키지 설치
RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \
    Rscript -e "if (!requireNamespace('optparse', quietly = TRUE)) install.packages('optparse', dependencies = TRUE)"

ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

# 작업 디렉토리 생성 및 설정
RUN mkdir -p /workspace/logs && \
    chmod 777 /workspace

# Snakefile 복사
COPY Snakefile /workspace/Snakefile

# scripts 폴더 복사
COPY scripts/ /scripts/

WORKDIR /workspace

# Entrypoint 스크립트 생성 (올바른 R 경로로 수정)
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

# 헬스체크 설정
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD test -f /opt/micromamba/envs/plugin_env/bin/python || exit 1

# 기본 명령어 설정
CMD ["/bin/bash"]