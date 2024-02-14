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

PYTHON_VERSION="python/3.9"
SNAKEMAKE_VERSION="snakemake/7.19.1"
SINGULARITY_VERSION="singularity/3.10.5"

set -eo pipefail

SCRIPTNAME="$0"
SCRIPTBASENAME=$(readlink -f $(basename $0))

#define cluster, partitions dependent on host
hostID=`echo $HOSTNAME`
PARTITIONS="norm,ccr"
cluster_specific_yaml="cluster_biowulf.yaml"
tools_specific_yaml="tools_biowulf.yaml"

# essential files
# these are relative to the workflows' base folder
# these are copied into the WORKDIR
ESSENTIAL_FILES="config/config.yaml config/samples.tsv config/contrasts.tsv config/fqscreen_config.conf config/multiqc_config.yaml config/rpackages.csv"
ESSENTIAL_FOLDERS="workflow/scripts annotation"
# set extra singularity bindings
EXTRA_SINGULARITY_BINDS="-B /data/CCBR_Pipeliner/,/lscratch"

# ## setting PIPELINE_HOME
PIPELINE_HOME=$(readlink -f $(dirname "$0"))
echo "Pipeline Dir: $PIPELINE_HOME"
# set snakefile
SNAKEFILE="${PIPELINE_HOME}/workflow/Snakefile"

function get_version_from_path() {
  PIPELINE_HOME=$1
  VERSION="$(cat $PIPELINE_HOME/VERSION)"
  echo -ne $VERSION
}

function get_version (){
echo "Version: $(get_version_from_path $PIPELINE_HOME)"
}

function usage() { 
# This function prints generic usage of the wrapper script.
  echo "${SCRIPTBASENAME}
  --> run CARLISLE
  Cut And Run anaLysIS pipeLinE

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
      *) runtest: run on cluster with included hg38 test dataset
  2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
  
  Optional Arguments: 
      *) -f / --force: Force flag will re-initialize a previously initialized workdir
  "
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
  if [[ -d $WORKDIR ]] && [[ -z $FORCEFLAG ]]; then 
    err "Folder $WORKDIR already exists! If you'd like to re-initialize use the -f/--force flag $FORCEFLAG|"
  fi
  mkdir -p $WORKDIR

  # copy essential files
  if [[ ! -d $WORKDIR/config ]]; then mkdir $WORKDIR/config; fi
  for f in $ESSENTIAL_FILES; do
    echo "Copying essential file: $f"
    fbn=$(basename $f)
    sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" -e "s/WORKDIR/${WORKDIR//\//\\/}/g" ${PIPELINE_HOME}/$f > $WORKDIR/config/$fbn
  done

  # rename config dependent on partition used
  cp ${PIPELINE_HOME}/resources/$cluster_specific_yaml $WORKDIR/config/cluster.yaml
  cp ${PIPELINE_HOME}/resources/$tools_specific_yaml $WORKDIR/config/tools.yaml

  # copy essential folders
  for f in $ESSENTIAL_FOLDERS;do
    rsync -avz --no-perms --no-owner --no-group --progress $PIPELINE_HOME/$f $WORKDIR/
  done

  #create log and stats folders
  if [ ! -d $WORKDIR/logs ]; then mkdir -p $WORKDIR/logs;echo "Logs Dir: $WORKDIR/logs";fi
  if [ ! -d $WORKDIR/stats ];then mkdir -p $WORKDIR/stats;echo "Stats Dir: $WORKDIR/stats";fi

  #create links for rose file
  if [[ ! -d $WORKDIR/annotation ]]; then 
    mkdir -p $WORKDIR/annotation
    cp ${PIPELINE_HOME}/annotation/* $WORKDIR/annotation/
  fi

  echo "Done Initializing $WORKDIR. You can now edit $WORKDIR/config/config.yaml and $WORKDIR/config/samples.tsv"

}

function check_essential_files() {
# Checks if files essential to start running the pipeline exist in the workdir

  if [ ! -d $WORKDIR ];then err "Folder $WORKDIR does not exist!"; fi
  for f in $ESSENTIAL_FILES; do
    fbn=$(basename $f)
    if [ ! -f $WORKDIR/config/$fbn ]; then err "Error: '${fbn}' file not found in $WORKDIR ... initialize first!";fi
  done
  for f in $ESSENTIAL_FOLDERS;do
    fbn=$(basename $f)
    if [ ! -d $WORKDIR/$fbn ]; then err "Error: '${fbn}' folder not found in $WORKDIR ... initialize first!";fi
  done
}

function set_singularity_binds(){
# this functions tries find what folders to bind
# biowulf specific
  echo "$PIPELINE_HOME" > ${WORKDIR}/tmp1
  echo "$WORKDIR" >> ${WORKDIR}/tmp1
  grep -o '\/.*' <(cat ${WORKDIR}/config/config.yaml ${WORKDIR}/config/samples.tsv)|tr '\t' '\n'|grep -v ' \|\/\/'|sort|uniq >> ${WORKDIR}/tmp1
  grep gpfs ${WORKDIR}/tmp1|awk -F'/' -v OFS='/' '{print $1,$2,$3,$4,$5}' |sort|uniq > ${WORKDIR}/tmp2
  grep -v gpfs ${WORKDIR}/tmp1|awk -F'/' -v OFS='/' '{print $1,$2,$3}'|sort|uniq > ${WORKDIR}/tmp3
  while read a;do readlink -f $a;done < ${WORKDIR}/tmp3 > ${WORKDIR}/tmp4
  binds=$(cat ${WORKDIR}/tmp2 ${WORKDIR}/tmp3 ${WORKDIR}/tmp4|sort|uniq |tr '\n' ',')
  rm -f ${WORKDIR}/tmp?
  binds=$(echo $binds|awk '{print substr($1,1,length($1)-1)}')
  SINGULARITY_BINDS="-B $EXTRA_SINGULARITY_BINDS,$binds"
}

function rescript(){
# recopy the scripts folder to WORKDIR
  check_essential_files
  rsync -avz --no-perms --no-owner --no-group --progress ${PIPELINE_HOME}/workflow/scripts/ $WORKDIR/scripts/
  echo "$WORKDIR/scripts folder has been updated!"
}

function runcheck(){
# Check "job-essential" files and load required modules
  check_essential_files
  module load $PYTHON_VERSION
  module load $SNAKEMAKE_VERSION
  # SINGULARITY_BINDS="$EXTRA_SINGULARITY_BINDS -B ${PIPELINE_HOME}:${PIPELINE_HOME} -B ${WORKDIR}:${WORKDIR}"
}

function controlcheck(){
# check controls are not listed in comparisons
  control_list=`awk '$3 ~ /Y/' ${WORKDIR}/config/samples.tsv | awk '{print $1}'`
  check1=`awk '{print $1}' ${WORKDIR}/config/contrasts.tsv`
  check2=`awk '{print $2}' ${WORKDIR}/config/contrasts.tsv`
  
  for sample_id in ${control_list[@]}; do
    if [[ $check1 == $sample_id || $check2 == $sample_id ]]; then 
      echo "Controls ($sample_id) cannot be listed in contrast.csv - update and re-run"
      echo "$check1 okkk $check2"
      exit 0
    fi
  done
}

function dryrun() {
# Dry-run
  runcheck
  controlcheck
  
  if [ ! -d ${WORKDIR}/logs/dryrun/ ]; then mkdir ${WORKDIR}/logs/dryrun/; fi

  if [ -f ${WORKDIR}/dryrun.log ]; then
    modtime=$(stat ${WORKDIR}/dryrun.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/dryrun.log ${WORKDIR}/logs/dryrun/dryrun.${modtime}.log
  fi
  run "--dry-run -r" | tee ${WORKDIR}/dryrun.log
}

function unlock() {
# Unlock the workdir if previous snakemake run ended abruptly
  runcheck
  run "--unlock"  
}

function dag() {
  runcheck
  module load graphviz
  snakemake -s $SNAKEFILE --configfile ${WORKDIR}/config/config.yaml --forceall --dag |dot -Teps > ${WORKDIR}/dag.eps
}

function runlocal() {
# If the pipeline is fired up on an interactive node (with sinteractive), this function runs the pipeline
  runcheck
  set_singularity_binds
  if [ "$SLURM_JOB_ID" == "" ];then err "runlocal can only be done on an interactive node"; fi
  module load $SINGULARITY_VERSION
  run "--dry-run" && echo "Dry-run was successful .... starting local execution" && run "local"
}

function runtest() {
  module load $PYTHON_VERSION
  module load $SNAKEMAKE_VERSION
  module load $SINGULARITY_VERSION
  check_essential_files
  sed -e "s/PIPELINE_HOME/${PIPELINE_HOME//\//\\/}/g" ${PIPELINE_HOME}/.test/samples.test.tsv > $WORKDIR/config/samples.tsv
  cp ${PIPELINE_HOME}/.test/contrasts.test.tsv $WORKDIR/config/contrasts.tsv
  check_essential_files
  run "--dry-run" && echo "Dry-run was successful .... starting slurm execution"
  runslurm
}

function runslurm() {
# Submit the execution of the pipeline to the biowulf job scheduler (slurm)
  runcheck
  set_singularity_binds
  run "--dry-run" && echo "Dry-run was successful .... submitting jobs to job-scheduler" && run "slurm"
}

function preruncleanup() {
# Cleanup function to rename/move files related to older runs to prevent overwriting them.
  echo "Running..."

  # check initialization
  check_essential_files 

  cd $WORKDIR
  ## Archive previous run files
  if [ -f ${WORKDIR}/logs/snakemake.log ];then 
    modtime=$(stat ${WORKDIR}/logs/snakemake.log |grep Modify|awk '{print $2,$3}'|awk -F"." '{print $1}'|sed "s/ //g"|sed "s/-//g"|sed "s/://g")
    mv ${WORKDIR}/logs/snakemake.log ${WORKDIR}/stats/snakemake.${modtime}.log

    if [ -f ${WORKDIR}/logs/snakemake.log.HPC_summary.txt ];then 
      mv ${WORKDIR}/logs/snakemake.log.HPC_summary.txt ${WORKDIR}/stats/snakemake.${modtime}.log.HPC_summary.txt
    fi
    if [ -f ${WORKDIR}/logs/snakemake.stats ];then 
      mv ${WORKDIR}/logs/snakemake.stats ${WORKDIR}/stats/snakemake.${modtime}.stats
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

  if [ "$1" == "local" ]; then
    preruncleanup

    snakemake -s $SNAKEFILE \
    --directory $WORKDIR \
    --printshellcmds \
    --use-singularity \
    --singularity-args "$SINGULARITY_BINDS" \
    --use-envmodules \
    --latency-wait 120 \
    --configfile ${WORKDIR}/config/config.yaml \
    --cores all \
    --stats ${WORKDIR}/logs/snakemake.stats \
    2>&1|tee ${WORKDIR}/logs/snakemake.log

    if [ "$?" -eq "0" ];then
      snakemake -s $SNAKEFILE \
      --report ${WORKDIR}/logs/runlocal_snakemake_report.html \
      --directory $WORKDIR \
      --configfile ${WORKDIR}/config/config.yaml 
    fi

  elif [ "$1" == "slurm" ];then
  
    preruncleanup

    cat > ${WORKDIR}/submit_script.sbatch << EOF
#!/bin/bash
#SBATCH --job-name="CARLISLE"
#SBATCH --mem=40g
#SBATCH --partition="$PARTITIONS"
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=2

    module load $PYTHON_VERSION
    module load $SNAKEMAKE_VERSION
    module load $SINGULARITY_VERSION

    cd \$SLURM_SUBMIT_DIR

    snakemake -s $SNAKEFILE \
    --directory $WORKDIR \
    --use-singularity \
    --singularity-args "$SINGULARITY_BINDS" \
    --use-envmodules \
    --printshellcmds \
    --latency-wait 120 \
    --configfile ${WORKDIR}/config/config.yaml \
    --cluster-config ${WORKDIR}/config/cluster.yaml \
    --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
    -j 500 \
    --rerun-incomplete \
    --keep-going \
    --restart-times 1 \
    --stats ${WORKDIR}/logs/snakemake.stats \
    2>&1|tee ${WORKDIR}/logs/snakemake.log

    if [ "\$?" -eq "0" ];then
      snakemake -s $SNAKEFILE \
      --directory $WORKDIR \
      --report ${WORKDIR}/logs/runslurm_snakemake_report.html \
      --configfile ${WORKDIR}/config/config.yaml 
    fi

    bash <(curl https://raw.githubusercontent.com/CCBR/Tools/master/Biowulf/gather_cluster_stats.sh 2>/dev/null) ${WORKDIR}/logs/snakemake.log > ${WORKDIR}/logs/snakemake.log.HPC_summary.txt
EOF

    sbatch ${WORKDIR}/submit_script.sbatch

  else # for unlock and dryrun 
    snakemake $1 -s $SNAKEFILE \
    --directory $WORKDIR \
    --use-envmodules \
    --printshellcmds \
    --latency-wait 120 \
    --configfile ${WORKDIR}/config/config.yaml \
    --cluster-config ${WORKDIR}/config/cluster.yaml \
    --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}" \
    -j 500 \
    --rerun-incomplete \
    --keep-going \
    --touch \
    --stats ${WORKDIR}/logs/snakemake.stats
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

function printbinds(){
  set_singularity_binds
  echo $SINGULARITY_BINDS
}

function main(){
# Main function which parses all arguments

  if [ $# -eq 0 ]; then usage; exit 1; fi

  for i in "$@"; do
    case $i in
        -m=*|--runmode=*)
          RUNMODE="${i#*=}"
        ;;
        -w=*|--workdir=*)
          WORKDIR="${i#*=}"
        ;;
        -f|--force)
          FORCEFLAG="ON"
        ;;
        -v|--version)
          get_version && exit 0
        ;;
        -h|--help)
          usage && exit 0;;
        *)
          err "Unknown argument!"    # unknown option
        ;;
    esac
  done

  WORKDIR=$(readlink -m "$WORKDIR")
  echo "Working Dir: $WORKDIR"

  case $RUNMODE in
    init) init && exit 0;;
    dryrun) dryrun && exit 0;;
    unlock) unlock && exit 0;;
    run) runslurm && exit 0;;
    runlocal) runlocal && exit 0;;
    reset) reset && exit 0;;
    runtest) runtest && exit 0;;
    rescript) rescript && exit 0;;               # hidden option for debugging
    printbinds) printbinds && exit 0;;           # hidden option
    *) err "Unknown RUNMODE \"$RUNMODE\"";;
  esac
}

# call the main function
main "$@"
