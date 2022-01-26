#!/usr/bin/env bash
# Author: Vishal Koparde, Ph.D.
# CCBR, NCI
# (c) 2021
#
# wrapper script to run the snakemake pipeline
# a) on an interactive node (runlocal) OR
# b) submit to the slurm load scheduler (run)
#
# DISCLAIMER: This wrapper only works on BIOWULF

PYTHON_VERSION="python/3.7"
SNAKEMAKE_VERSION="snakemake/5.24.1"

set -eo pipefail
module purge

SCRIPTNAME="$0"
SCRIPTBASENAME=$(readlink -f $(basename $0))

# set extra singularity bindings
EXTRA_SINGULARITY_BINDS="-B /data/CCBR_Pipeliner/:/data/CCBR_Pipeliner/"

function get_git_commitid_tag() {
# This function gets the latest git commit id and tag
# Input is PIPELINE_HOME folder which is a git folder
  cd $1
  gid=$(git rev-parse HEAD)
  tag=$(git describe --tags $gid 2>/dev/null)
  echo -ne "$gid\t$tag"
}

# ## setting PIPELINE_HOME
PIPELINE_HOME=$(readlink -f $(dirname "$0"))
echo "Pipeline Dir: $PIPELINE_HOME"
# set snakefile
SNAKEFILE="${PIPELINE_HOME}/workflow/Snakefile"

# get github commit tag
GIT_COMMIT_TAG=$(get_git_commitid_tag $PIPELINE_HOME)
echo "Git Commit/Tag: $GIT_COMMIT_TAG"

function usage() { 

# This function prints generic usage of the wrapper script.

    cat << EOF
${SCRIPTBASENAME}
--> run <Your Pipeline Short Name>

USAGE:
  bash ${SCRIPTNAME} -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
Required Arguments:
1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) run : run with slurm
    *) reset : DELETE workdir dir and re-init it
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
    *) runlocal : run without submitting to sbatch
2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
EOF
}

function err() { 

# This function print a error message (argument 1), the usage and then exits with non-zero status

cat <<< "
#
#
#
  $@
#
#
#
" && usage && exit 1 1>&2; 
}


function init() {

# This function initializes the workdir by:
# 1. creating the working dir
# 2. copying essential files like config.yaml and samples.tsv into the workdir
# 3. setting up logs and stats folders

# create output folder
if [ -d $WORKDIR ];then err "Folder $WORKDIR already exists!"; fi
mkdir -p $WORKDIR

# copy config and samples files
sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $WORKDIR/config.yaml
cp ${PIPELINE_HOME}/config/samples.tsv $WORKDIR/

#create log and stats folders
if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

echo "Done Initializing $WORKDIR. You can now edit $WORKDIR/config.yaml and $WORKDIR/samples.tsv"

}

function check_essential_files() {

# Checks if files essential to start running the pipeline exist in the workdir

  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!"; fi
  for f in config.yaml samples.tsv; do
    if [ ! -f $WORKDIR/$f ]; then err "Error: '${f}' file not found in workdir ... initialize first!";fi
  done

}

function reconfig(){
# Rebuild config file and replace the config.yaml in the WORKDIR
# this is only for dev purposes when new key-value pairs are being added to the config file

  check_essential_files
  sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/config/config.yaml > $WORKDIR/config.yaml
  echo "$WORKDIR/config.yaml has been updated!"

}

function runcheck(){
# Check "job-essential" files and load required modules

  check_essential_files
  module load $PYTHON_VERSION
  module load $SNAKEMAKE_VERSION
  SINGULARITY_BINDS="$EXTRA_SINGULARITY_BINDS -B ${PIPELINE_HOME}:${PIPELINE_HOME} -B ${WORKDIR}:${WORKDIR}"

}

function dryrun() {
# Dry-run

  runcheck
  run "--dry-run"
}

function unlock() {
# Unlock the workdir if previous snakemake run ended abruptly

  runcheck
  run "--unlock"  
}

function runlocal() {
# If the pipeline is fired up on an interactive node (with sinteractive), this function runs the pipeline

  runcheck
  if [ "$SLURM_JOB_ID" == "" ];then err "runlocal can only be done on an interactive node"; fi
  module load singularity
  run "local"
}

function runslurm() {
# Submit the execution of the pipeline to the biowulf job scheduler (slurm)

  runcheck
  run "slurm"
}

function preruncleanup() {
# Cleanup function to rename/move files related to older runs to prevent overwriting them.

  echo "Running..."

  # check initialization
  check_essential_files 

  cd $WORKDIR
  ## Archive previous run files
  if [ -f ${WORKDIR}/snakemake.log ];then 
    modtime=$(stat ${WORKDIR}/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/snakemake.log ${WORKDIR}/stats/snakemake.${modtime}.log
    if [ -f ${WORKDIR}/snakemake.log.HPC_summary.txt ];then 
      mv ${WORKDIR}/snakemake.log.HPC_summary.txt ${WORKDIR}/stats/snakemake.${modtime}.log.HPC_summary.txt
    fi
    if [ -f ${WORKDIR}/snakemake.stats ];then 
      mv ${WORKDIR}/snakemake.stats ${WORKDIR}/stats/snakemake.${modtime}.stats
    fi
  fi
  nslurmouts=$(find ${WORKDIR} -maxdepth 1 -name "slurm-*.out" |wc -l)
  if [ "$nslurmouts" != "0" ];then
    for f in $(ls ${WORKDIR}/slurm-*.out);do mv ${f} ${WORKDIR}/logs/;done
  fi

}


function run() {
# RUN function
# argument1 can be:
# 1. local or
# 2. dryrun or
# 3. unlock or
# 4. slurm

  if [ "$1" == "local" ];then

  preruncleanup

  snakemake -s $SNAKEFILE \
  --directory $WORKDIR \
  --printshellcmds \
  --use-singularity \
  --singularity-args "$SINGULARITY_BINDS" \
  --use-envmodules \
  --latency-wait 120 \
  --configfile ${WORKDIR}/config.yaml \
  --cores all \
  --stats ${WORKDIR}/snakemake.stats \
  2>&1|tee ${WORKDIR}/snakemake.log

  if [ "$?" -eq "0" ];then
    snakemake -s $SNAKEFILE \
    --report ${WORKDIR}/runlocal_snakemake_report.html \
    --directory $WORKDIR \
    --configfile ${WORKDIR}/config.yaml 
  fi

  elif [ "$1" == "slurm" ];then
  
  preruncleanup
# if QOS is other than "global" and is supplied in the cluster.json file then add " --qos={cluster.qos}" to the 
# snakemake command below
  cat > ${WORKDIR}/submit_script.sbatch << EOF
#!/bin/bash
#SBATCH --job-name="CCBRPipeline"
#SBATCH --mem=10g
#SBATCH --partition="norm"
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2

module load $PYTHON_VERSION
module load $SNAKEMAKE_VERSION
module load singularity

cd \$SLURM_SUBMIT_DIR

snakemake -s $SNAKEFILE \
--directory $WORKDIR \
--use-singularity \
--singularity-args "$SINGULARITY_BINDS" \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile ${WORKDIR}/config.yaml \
--cluster-config ${PIPELINE_HOME}/resources/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
-j 500 \
--rerun-incomplete \
--keep-going \
--stats ${WORKDIR}/snakemake.stats \
2>&1|tee ${WORKDIR}/snakemake.log

if [ "\$?" -eq "0" ];then
  snakemake -s $SNAKEFILE \
  --directory $WORKDIR \
  --report ${WORKDIR}/runslurm_snakemake_report.html \
  --configfile ${WORKDIR}/config.yaml 
fi

bash <(curl https://raw.githubusercontent.com/CCBR/Tools/master/Biowulf/gather_cluster_stats.sh 2>/dev/null) ${WORKDIR}/snakemake.log > ${WORKDIR}/snakemake.log.HPC_summary.txt

EOF

  sbatch ${WORKDIR}/submit_script.sbatch

  else # for unlock and dryrun 

snakemake $1 -s $SNAKEFILE \
--directory $WORKDIR \
--use-envmodules \
--printshellcmds \
--latency-wait 120 \
--configfile ${WORKDIR}/config.yaml \
--cluster-config ${PIPELINE_HOME}/config/cluster.json \
--cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
-j 500 \
--rerun-incomplete \
--keep-going \
--stats ${WORKDIR}/snakemake.stats

  fi

}

function reset() {
# Delete the workdir and re-initialize it

  echo "Working Dir: $WORKDIR"
  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!";fi
  echo "Deleting $WORKDIR"
  rm -rf $WORKDIR
  echo "Re-Initializing $WORKDIR"
  init
}


function main(){
# Main function which parses all arguments

  if [ $# -eq 0 ]; then usage; exit 1; fi

  for i in "$@"
  do
  case $i in
      -m=*|--runmode=*)
        RUNMODE="${i#*=}"
      ;;
      -w=*|--workdir=*)
        WORKDIR="${i#*=}"
      ;;
      *)
        err "Unknown argument!"    # unknown option
      ;;
  esac
  done
  WORKDIR=$(readlink -f "$WORKDIR")
  echo "Working Dir: $WORKDIR"

  case $RUNMODE in
    init) init && exit 0;;
    dryrun) dryrun && exit 0;;
    unlock) unlock && exit 0;;
    run) runslurm && exit 0;;
    runlocal) runlocal && exit 0;;
    reset) reset && exit 0;;
    dry) dryrun && exit 0;;                      # hidden option
    local) runlocal && exit 0;;                  # hidden option
    reconfig) reconfig && exit 0;;               # hidden option for debugging
    *) err "Unknown RUNMODE \"$RUNMODE\"";;
  esac
}

# call the main function

main "$@"
