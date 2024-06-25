# DO NOT CHANGE
from 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

workdir /tmp/docker-build/work/

shell [ \
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
env TZ='Etc/UTC'
env LANG='en_US.UTF-8'

arg DEBIAN_FRONTEND=noninteractive

# Get secure apt key for CRAN
run apt-get update --yes && \
    apt-get install --yes dirmngr apt-transport-https gnupg2 && \
    gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
    gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' | tee /etc/apt/trusted.gpg.d/cran_debian_key.asc

# install R requirements
run echo "deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/" | tee /etc/apt/sources.list.d/cran.list && \
    apt-get update --yes && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes r-base r-base-dev libxml2-dev libcurl4-openssl-dev libssl-dev wget

COPY latch_environment.R /tmp/docker-build/work/latch_environment.R
RUN Rscript /tmp/docker-build/work/latch_environment.R

# COPY clusterProfiler_4.12.0.tar.gz /tmp/docker-build/work/clusterProfiler_4.12.0.tar.gz
# RUN R -e "install.packages('/tmp/docker-build/work/clusterProfiler_4.12.0.tar.gz', repos = NULL, type='source')"

# COPY KEGG.db_2.8.0.tar.gz /tmp/docker-build/work/KEGG.db_2.8.0.tar.gz
# RUN R -e "install.packages('/tmp/docker-build/work/KEGG.db_2.8.0.tar.gz', repos = NULL, type='source')"

RUN pip install pandas

# Latch SDK
# DO NOT REMOVE
run pip install latch==2.47.3
run mkdir /opt/latch
run apt-get update && apt-get install -y default-jre-headless


# Copy workflow data (use .dockerignore to skip files)
copy . /root/

# Latch nextflow workflow entrypoint
# DO NOT CHANGE
run ln -s /root/.latch/bin/nextflow /root/nextflow
run ln -s /root/.latch/.nextflow /root/.nextflow


# Latch workflow registration metadata
# DO NOT CHANGE
arg tag
# DO NOT CHANGE
env FLYTE_INTERNAL_IMAGE $tag

workdir /root
