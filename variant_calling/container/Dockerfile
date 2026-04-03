FROM mambaorg/micromamba:1.5.8

# Set environment name
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Copy environment file
COPY environment.yml /tmp/environment.yml

# Install dependencies
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Set working directory
WORKDIR /opt

# Copy needLR code
COPY needLR_v3.5_local /opt/needLR_v3.5_local

# Set working directory inside needLR
WORKDIR /opt/needLR_v3.5_local

# Make scripts executable
RUN chmod +x needLR_v3.5_*.sh

# Default shell
SHELL ["/bin/bash", "-c"]
