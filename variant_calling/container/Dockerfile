FROM mambaorg/micromamba:1.5.10

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    procps \
    findutils \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_SUBDIR=linux-64

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    truvari \
    bedtools \
    bcftools \
    htslib \
    tabix \
    pandas \
    numpy \
    && micromamba clean --all --yes

WORKDIR /opt
COPY needLR_v3.5_local /opt/needLR_v3.5_local

WORKDIR /opt/needLR_v3.5_local
RUN chmod +x needLR_v3.5_*.sh

ENV PATH=/opt/conda/bin:$PATH

CMD ["/bin/bash"]
