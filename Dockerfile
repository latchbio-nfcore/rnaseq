# DO NOT CHANGE
FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base-nextflow:v2.1.2

WORKDIR /tmp/docker-build/work/

SHELL [ \
    "/usr/bin/env", "bash", \
    "-o", "errexit", \
    "-o", "pipefail", \
    "-o", "nounset", \
    "-o", "verbose", \
    "-o", "errtrace", \
    "-O", "inherit_errexit", \
    "-O", "shift_verbose", \
    "-c" \
]
ENV TZ='Etc/UTC'
ENV LANG='en_US.UTF-8'

ARG DEBIAN_FRONTEND=noninteractive

# Install rig the R installation manager
RUN \
    curl \
        --location \
        --fail \
        --remote-name \
        https://github.com/r-lib/rig/releases/download/latest/rig-linux-latest.tar.gz && \
    tar \
        --extract \
        --gunzip \
        --file rig-linux-latest.tar.gz \
        --directory /usr/local/ && \
    rm rig-linux-latest.tar.gz

# Install R
RUN rig add release

COPY latch_environment.R /tmp/docker-build/work/latch_environment.R
RUN Rscript /tmp/docker-build/work/latch_environment.R

# Install needed packages
RUN pip install pandas && \
    apt-get update && apt-get install -y pigz

# Latch SDK
# DO NOT REMOVE
RUN pip install latch==2.53.0
RUN mkdir /opt/latch

# Copy workflow data (use .dockerignore to skip files)
COPY . /root/

# Latch workflow registration metadata
# DO NOT CHANGE
ARG tag
# DO NOT CHANGE
ENV FLYTE_INTERNAL_IMAGE $tag

WORKDIR /root
