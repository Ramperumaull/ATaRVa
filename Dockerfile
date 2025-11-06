#################################################################
# Dockerfile
#
# Software:         ATaRVa
# Software Version: 0.4.0
# Description:      ATaRVa image
# Summary:          ATaRVa is a tandem repeat genotyper, specially designed for long read data.
# Website:          https://github.com/SowpatiLab/ATaRVa
# License:          MIT
# Base Image:       python:3.13
# Tags:             Genomics, Next-Generation Sequencing, Bioinformatics, Tandem repeats, STR, VNTR, repeats, ONT, PacBio, microsatellites, long reads
# Maintainers:      Akshay Kumar Avvaru <avvaruakshay@gmail.com>, Abishek Kumar <abishekks@csirccmb.org>
# Build Cmd:        docker build -f Dockerfile -t atarva .
# Run Cmd:          docker run -i -t --rm atarva
#################################################################

# Getting python from Docker Hub
FROM python:3.13

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
