#!/usr/bin/env bash
# usage:
#   ./install.sh new/path/to/install
# examples
#   ./install.sh .v2.4.0
#   /data/CCBR_Pipeliner/Pipelines/CARLISLE/dev/install.sh /data/CCBR_Pipeliner/Pipelines/CARLISLE/.v2.5.0-a
set -euo pipefail

VERSION=$1
mkdir -p ${VERSION}/bin
INSTALL_PATH=$(readlink -f ${VERSION}/bin)
DIRNAME=$(readlink -f $(dirname $0))

if [ -n "$(ls -A $INSTALL_PATH 2>/dev/null)" ]
then
    echo "ERROR: directory not empty: ${INSTALL_PATH}"
    echo $(ls $INSTALL_PATH)
    exit 1
fi

# copy entire repo, including dotfiles
cp -r ${DIRNAME}/. ${INSTALL_PATH}

# export path
if [[ ":$PATH:" != *":${INSTALL_PATH}:"* ]];then
    export PATH="${PATH}:${INSTALL_PATH}"
fi

echo "${INSTALL_PATH}"
