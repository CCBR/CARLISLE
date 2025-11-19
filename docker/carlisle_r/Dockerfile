FROM nciccbr/ccbr_ubuntu_base_20.04:v7

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install conda packages
COPY environment.yml /data2/
ENV CONDA_ENV=carlisle
RUN mamba env create -n ${CONDA_ENV} -f /data2/environment.yml && \
    echo "conda activate ${CONDA_ENV}" > ~/.bashrc
ENV PATH="/opt2/conda/envs/${CONDA_ENV}/bin:$PATH"
ENV R_LIBS_USER=/opt2/conda/lib/R/library/

# install ELBOW manually, fails with mamba
RUN wget --no-check-certificate https://bioconductor.org/packages/3.4/bioc/src/contrib/ELBOW_1.10.0.tar.gz && \
    R -e 'install.packages("ELBOW_1.10.0.tar.gz", repos = NULL, type="source", INSTALL_opts = "--no-lock")'

# write R package versions
RUN R -s -e "readr::write_tsv(tibble::as_tibble(installed.packages()), '/data2/r-packages.tsv')"

# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}

# cleanup
WORKDIR /data2
RUN apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
