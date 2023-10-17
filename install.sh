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

# copy essential files & directories

## carlisle CLI script
cp $DIRNAME/carlisle $INSTALL_PATH/

## all config & workflow scripts;
for subdir in config workflow/scripts
do  
    mkdir -p ${INSTALL_PATH}/$subdir
    cp ${DIRNAME}/$subdir/* ${INSTALL_PATH}/$subdir
done
cp ${DIRNAME}/workflow/Snakefile ${INSTALL_PATH}/workflow

## selected resources
mkdir -p ${INSTALL_PATH}/resources/
cp ${DIRNAME}/resources/*.yaml ${INSTALL_PATH}/resources/

# export path
if [[ ":$PATH:" != *":${INSTALL_PATH}:"* ]];then
    export PATH="${PATH}:${INSTALL_PATH}"
fi

echo "${INSTALL_PATH}"
