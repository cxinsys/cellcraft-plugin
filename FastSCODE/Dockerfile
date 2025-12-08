FROM nvidia/cuda:12.1.0-cudnn8-devel-ubuntu20.04

# 비대화형 설치 설정
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Seoul

# 기본 패키지 설치
RUN apt-get update && apt-get install -y \
    build-essential gcc g++ gfortran make \
    libssl-dev libcurl4-openssl-dev libxml2-dev \
    libjpeg-dev libpng-dev libfreetype6-dev libtiff-dev \
    libx11-dev xorg-dev libxt-dev libglu1-mesa-dev \
    libharfbuzz-dev libfribidi-dev \
    libglpk-dev \
    curl wget unzip git \
    python3.10 python3-pip python3-venv \
    software-properties-common \
    gnupg ca-certificates \
    default-jdk \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Java 환경변수 설정
ENV JAVA_HOME=/usr/lib/jvm/default-java
ENV PATH=$JAVA_HOME/bin:$PATH

# Java 버전 확인
RUN java -version && javac -version

# R 설치 (특정 버전)
RUN apt-get update && \
    apt-get install -y dirmngr gpg-agent && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
    tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get update

# R 최신 버전 설치
RUN apt-get install -y r-base r-base-dev r-recommended && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# R 버전 확인
RUN R --version | head -1

# dependency 폴더 복사
COPY dependency/ /workspace/dependency/

# 디버깅: dependency 폴더 내용 확인
RUN echo '[DEBUG] Contents of /workspace/dependency folder:' && \
    ls -la /workspace/dependency/ && \
    echo '[DEBUG] Detailed file information:' && \
    find /workspace/dependency -type f -exec ls -lh {} \; && \
    echo '[DEBUG] File contents preview:' && \
    for file in /workspace/dependency/*.txt /workspace/dependency/*.yml /workspace/dependency/*.yaml /workspace/dependency/*.lock; do \
        if [ -f "$file" ]; then \
            echo "=== Contents of $(basename $file) ==="; \
            head -20 "$file" || echo "Could not read $file"; \
            echo ""; \
        fi; \
    done

# Micromamba 설치
RUN mkdir -p /usr/local/bin && \
    cd /tmp && \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj && \
    cp bin/micromamba /usr/local/bin/micromamba && \
    chmod +x /usr/local/bin/micromamba && \
    rm -rf /tmp/bin /tmp/info && \
    # 설치 확인
    /usr/local/bin/micromamba --version

# 환경 변수 설정
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV MAMBA_EXE=/usr/local/bin/micromamba
ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$MAMBA_ROOT_PREFIX/bin:$PATH

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

# R 패키지 설치 - 시스템 라이브러리에 직접 설치
RUN Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org'), download.file.method = 'libcurl')"

# renv 설치
RUN Rscript -e "install.packages('renv')"

# 추가 필수 패키지 확인 및 설치
RUN Rscript -e " \
    required_pkgs <- c('optparse', 'jsonlite', 'readr', 'dplyr', 'ggplot2', 'pheatmap', 'plotly'); \
    missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[,'Package']]; \
    if (length(missing_pkgs) > 0) { \
        cat('Installing missing packages:', paste(missing_pkgs, collapse=', '), '\\n'); \
        # BiocManager가 필요한 경우를 위해 먼저 설치 \
        if (!requireNamespace('BiocManager', quietly = TRUE)) { \
            install.packages('BiocManager') \
        }; \
        # 일반 CRAN 패키지 설치 시도 \
        cran_pkgs <- missing_pkgs[missing_pkgs %in% c('optparse', 'jsonlite', 'readr', 'dplyr', 'ggplot2', 'plotly')]; \
        if (length(cran_pkgs) > 0) { \
            install.packages(cran_pkgs, dependencies = TRUE) \
        }; \
        # pheatmap은 의존성이 많으므로 별도 처리 \
        if ('pheatmap' %in% missing_pkgs) { \
            install.packages('pheatmap', dependencies = TRUE) \
        } \
    }; \
    # 설치 확인 \
    cat('\\nChecking installed packages:\\n'); \
    for (pkg in required_pkgs) { \
        if (pkg %in% installed.packages()[,'Package']) { \
            cat(paste(pkg, 'is installed\\n')) \
        } else { \
            cat(paste('WARNING:', pkg, 'is NOT installed\\n')) \
        } \
    } \
    "

# renv 자동 활성화 비활성화를 위한 환경 변수
ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

# R 환경 변수 설정
ENV R_HOME=/usr/lib/R

# 작업 디렉토리 생성 및 설정
RUN mkdir -p /workspace/logs && \
    chmod 777 /workspace

# Snakefile 복사
COPY Snakefile /workspace/Snakefile
COPY visualization_Snakefile /workspace/visualization_Snakefile

# scripts 폴더 복사
COPY scripts/ /scripts/

WORKDIR /workspace

# Entrypoint 스크립트 생성
RUN echo '#!/bin/bash' > /entrypoint.sh
RUN echo 'export MAMBA_ROOT_PREFIX=/opt/micromamba' >> /entrypoint.sh
RUN echo 'export MAMBA_EXE=/usr/local/bin/micromamba' >> /entrypoint.sh
RUN echo 'export PATH=$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$PATH' >> /entrypoint.sh
RUN echo '# Set R environment' >> /entrypoint.sh
RUN echo 'export R_HOME=/usr/lib/R' >> /entrypoint.sh
RUN echo 'export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /entrypoint.sh
RUN echo '# Activate micromamba environment' >> /entrypoint.sh
RUN echo 'eval "$($MAMBA_EXE shell activate -s bash -p $MAMBA_ROOT_PREFIX plugin_env)" 2>/dev/null || true' >> /entrypoint.sh
RUN echo 'cd /workspace' >> /entrypoint.sh
RUN echo 'exec "$@"' >> /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

# 헬스체크 설정
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD test -f /opt/micromamba/envs/plugin_env/bin/python || exit 1

# 기본 명령어 설정
CMD ["/bin/bash"]