# Navier-Stokes Solvers Docker Image
# Provides complete build environment with all dependencies

FROM ubuntu:22.04

LABEL maintainer="Diogo Ribeiro <dfr@esmad.ipp.pt>"
LABEL description="High-performance Navier-Stokes solvers with spectral and finite difference methods"
LABEL version="1.0.0"

# Prevent interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    cmake \
    git \
    libfftw3-dev \
    libfftw3-3 \
    libopenblas-dev \
    liblapack-dev \
    python3 \
    python3-pip \
    python3-numpy \
    python3-matplotlib \
    python3-scipy \
    valgrind \
    gdb \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Copy source code
COPY . /workspace/

# Build the solvers
RUN make clean && make all

# Set up Python environment for visualization
RUN pip3 install --no-cache-dir \
    numpy>=1.24.0 \
    scipy>=1.10.0 \
    matplotlib>=3.7.0 \
    pandas>=2.0.0 \
    seaborn>=0.12.0

# Create output directory
RUN mkdir -p /workspace/results

# Default command shows help
CMD ["make", "help"]

# Usage examples:
# Build: docker build -t navier-stokes-solvers .
# Run FD solver: docker run -v $(pwd)/results:/workspace/results navier-stokes-solvers ./bin/ns_fd_solver
# Run spectral: docker run -v $(pwd)/results:/workspace/results navier-stokes-solvers ./bin/ns_spectral_solver
# Interactive: docker run -it navier-stokes-solvers /bin/bash
