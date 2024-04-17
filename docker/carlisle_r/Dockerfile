FROM nciccbr/ccbr_ubuntu_base_20.04:v6

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install conda packages
COPY packages.txt /data2/
RUN mamba install \
    --no-channel-priority \
    -c bioconda -c conda-forge -c r \
    --file /data2/packages.txt
ENV PATH="/opt2/conda/bin:$PATH"
ENV R_LIBS_USER=/opt2/conda/lib/R/library/
# install ELBOW manually, fails with mamba
RUN wget https://bioconductor.riken.jp/packages/3.4/bioc/src/contrib/ELBOW_1.10.0.tar.gz && \
    R -e 'install.packages("ELBOW_1.10.0.tar.gz", repos = NULL, type="source", INSTALL_opts = "--no-lock")'

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
