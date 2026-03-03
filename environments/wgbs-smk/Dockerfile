FROM ubuntu:22.04

# Use bash
SHELL ["/bin/bash", "-c"]

# Environment variables
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive
ENV BCFTOOLS_PLUGINS=/usr/local/bcftools/plugins
ENV PATH=/root/.local/bin:${PATH}
ENV PATH=/root/.cargo/bin:${PATH}
RUN echo -e $PATH | tr ":" "\n"

# Update and install packages
RUN apt-get -y update && \
    apt-get -y install \
        build-essential \
        python3-dev \
        python3-full \
        uuid-runtime \
        pkg-config \
        parallel \
        autoconf \
        gcc \
	      git \
        make \
        curl \
        wget \
        lz4 \
        jq \
        rsync \
        cmake \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        lbzip2 \
        libexpat1-dev \
        libncurses5-dev \
        libfreetype6-dev \
        libfontconfig1-dev

# Clean up
RUN apt-get -y clean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt-get/lists/*

# Install python deps with uv
ENV VIRTUAL_ENV=/python-venv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
RUN uv venv ${VIRTUAL_ENV}
ENV PATH=${VIRTUAL_ENV}/bin:${PATH}
RUN uv pip install --no-cache-dir \
    beautifulsoup4 \
    ipykernel \
    jinja2 \
    lxml \
    multiprocess \
    pandas \
    polars \
    pysam \
    qc-utils \
    snakemake

# Install rust tool chain
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

# Install gemBS-rs (WGBS-seq processing pipeline)
RUN git clone --recursive https://github.com/heathsc/gemBS-rs.git && \
    cd /gemBS-rs && \
    make gemBS_config.mk && \
    make install
RUN cp /gemBS-rs/rust/bin/{read_filter,bs_call,mextr,dbsnp_index,gemBS} /usr/local/bin
RUN rm -rf /gemBS-rs

# Install rust tools
RUN cargo install -j $(nproc) bigtools
# # Not needed anymore and takes a long time so commenting out for now
# RUN cargo install -j 32 qsv --version 4.0.0 --locked --bin qsvlite -F lite

# Environment variables (compiler)
ENV CFLAGS="-O3 -march=x86-64-v4"
ENV CXXFLAGS="-O3 -march=x86-64-v4"

# Install SuiteSparse v7.6.0 which contains CHOLMOD5
RUN apt-get install -y libopenblas-dev libgmp-dev libmpfr-dev && \
    git clone --branch v7.6.0 https://github.com/DrTimothyAldenDavis/SuiteSparse.git /tmp/SuiteSparse && \
    cd /tmp/SuiteSparse && \
    cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local && \
    cmake --build build -j$(nproc) && \
    cmake --install build && \
    ldconfig && \
    rm -rf /tmp/SuiteSparse

# Build variables
ARG GEM_MAPPER_VERSION=master

# Install methylation aware mapper, GEM3
RUN git clone --recursive https://github.com/smarco/gem3-mapper.git -b ${GEM_MAPPER_VERSION} && \
    cd gem3-mapper && \
    chmod +x configure && \
    ./configure --prefix=/usr/local --enable-cuda=no && \
    make FLAGS_OPT="-O3 -march=x86-64-v4"
RUN cp /gem3-mapper/bin/{gem-indexer,gem-retriever,gem-mapper} /usr/local/bin
RUN rm -rf /gem3-mapper

# Install pigz (faster gzip impl)
RUN curl -sSLf https://zlib.net/pigz/pigz.tar.gz | tar xzf - && \
    cd /pigz && \
    make && \
    cp pigz unpigz /usr/local/bin
RUN rm -rf /pigz

# Install seqkit
RUN curl -fsSL https://github.com/shenwei356/seqkit/releases/download/v2.11.0/seqkit_linux_amd64.tar.gz | tar xzf - && \
    cp seqkit /usr/local/bin && \
    chmod a+x /usr/local/bin/seqkit

# Install bedtools
RUN curl -s https://api.github.com/repos/arq5x/bedtools2/releases/latest | \
    jq -r '.assets[] | select(.name | test("bedtools-.*\\.tar\\.gz")) | .browser_download_url' | \
    xargs curl -fsSL | tar xzf - && \
    cd bedtools2 && \
    make static && \
    mv bin/bedtools.static /usr/local/bin/bedtools
RUN rm -rf bedtools2

# Install samtools stuff
RUN curl -s https://api.github.com/repos/samtools/htslib/releases/latest | \
    jq -r '.assets[] | select(.name | test("htslib-.*\\.tar\\.bz2")) | .browser_download_url' | \
    xargs curl -fsSL | tar xjf - && \
    cd /htslib-* && \
    ./configure --prefix=/usr && \
    make -j$(nproc) && \
    make install
RUN curl -s https://api.github.com/repos/samtools/samtools/releases/latest | \
    jq -r '.assets[] | select(.name | test("samtools-.*\\.tar\\.bz2")) | .browser_download_url' | \
    xargs curl -fsSL | tar xjf - && \
    cd /samtools-* && \
    ./configure && \
    make -j$(nproc) && \
    make install
RUN ln -sf /usr/include/suitesparse/cholmod.h /usr/include/cholmod.h && \
    ln -sf /usr/include/suitesparse/SuiteSparse_config.h /usr/include/SuiteSparse_config.h
RUN curl -s https://api.github.com/repos/samtools/bcftools/releases/latest | \
    jq -r '.assets[] | select(.name | test("bcftools-.*\\.tar\\.bz2")) | .browser_download_url' | \
    xargs curl -fsSL | tar xjf - && \
    cd /bcftools-* && \
    ./configure && \
    /bin/rm -f plugins/{score.{c,h},{munge,liftover,metal,blup}.c,pgs.{c,mk}} && \
    wget -P plugins http://raw.githubusercontent.com/freeseek/score/master/{score.{c,h},{munge,liftover,metal,blup}.c,pgs.{c,mk}} && \
    make CPPFLAGS="-I/usr/local/include/suitesparse" -j$(nproc) && \
    mkdir -p $BCFTOOLS_PLUGINS && \
    /bin/cp plugins/{munge,liftover,score,metal,pgs,blup}.so $BCFTOOLS_PLUGINS && \
    wget -P $BCFTOOLS_PLUGINS http://raw.githubusercontent.com/freeseek/score/master/assoc_plot.R && \
    chmod a+x $BCFTOOLS_PLUGINS/assoc_plot.R && \
    make install
RUN rm -rf /samtools-* /bcftools-* /htslib-*

# Install ucsc bedGraphToBigWig
RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedGraphToBigWig /usr/local/bin

# Test everything
RUN which gem-mapper && \
    which pigz && \
    which unpigz && \
    which bigtools && \
    which samtools && \
    which bcftools && \
    which gemBS && \
    which read_filter && \
    which bs_call && \
    which mextr

CMD ["/bin/bash"]
