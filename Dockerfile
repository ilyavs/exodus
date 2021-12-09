FROM ubuntu:16.04

# Set the locale. This works despite the perl warning.
ENV LANG=en_US.UTF-8
RUN apt-get clean && apt-get update && apt-get install -y locales && \
    sed -i -e "s/# $LANG.*/$LANG UTF-8/" /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=$LANG

# Install git
RUN apt-get update && \
    apt-get install -y git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log /tmp/* /var/tmp/*

# Install conda
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3.7 \
    && conda update  conda \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log /tmp/* /var/tmp/* \
    && conda clean --all --yes \
    && conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge

# Install mamba
RUN conda install -y mamba=0.19.1 && mamba clean --all --yes

# Install dependencies
RUN mamba install -y snakemake=5.20.1 fastqc=0.11.9 \
    bowtie2=2.4.2 samtools=1.11 \
    bwa=0.7.17 pysam=0.16.0.1 \
    pandas=1.1.3 numpy=1.19.2 \
    lxml=4.6.1 \
    && mamba clean --all --yes

# Clone the public repo
RUN git clone https://github.com/ilyavs/exodus.git