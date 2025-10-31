# Base Image
FROM biocontainers/biocontainers:latest

# Metadata
LABEL base_image="biocontainers:latest" \
        version="1" \
        software="ATaRVa" \
        software.version="0.4.0" \
        about.summary="ATaRVa is a tandem repeat genotyper, specially designed for long read data." \
        about.home="https://github.com/SowpatiLab/ATaRVa.git" \
        about.documentation="https://github.com/SowpatiLab/ATaRVa.git"
        about.license="MIT" \
        about.tag="Genomics, Next-Generation Sequencing, Bioinformatics, Tandem repeats, STR, VNTR, repeats, ONT, PacBio, microsatellites, long reads" \
        maintainers="Akshay Kumar Avvaru <avvaruakshay@gmail.com>, Abishek Kumar <abishekks@csirccmb.org>"

# Getting python from Docker Hub
FROM python:3.9.5

# Setting working directory
WORKDIR /app

# Install the dependencies
RUN apt-get update && \
        apt-get install -y --no-install-recommends \
        gcc \
        build-essential \
        libssl-dev \
        git

# Install ATaRVa
RUN git clone https://github.com/SowpatiLab/ATaRVa.git /app

RUN pip install --no-cache-dir . && \
        rm -rf /var/lib/apt/lists/*

# Define the default command to run your application
ENTRYPOINT ["python", "-m", "ATARVA.core"]
