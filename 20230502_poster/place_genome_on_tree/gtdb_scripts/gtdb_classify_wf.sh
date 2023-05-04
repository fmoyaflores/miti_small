#!/bin/bash -x
# An 'r5.12xlarge' instance was used to run this script.
# In order to get the genomes palces on the full tree, the --full_tree flag was used.
# This requires, per the documentation, at least 320GB of memory.

set -euo pipefail
# update this to reflect the path to the GTDB database in your system
GTDBTK_DATA_PATH='/mnt/efs/databases/GTDB/release207_v2' 
# full path to the directory containing the assemblies
assembly_dir=$1 
# project name
project="TransposonLibrary_20210331"
# number of threads to use for pplacer
pplacer_threads=30
# number of threads to use for GTDB-Tk
cpus=45
# file extension of the assemblies
ext="fna"

# Docker specific variables
docker_image=quay.io/biocontainers/gtdbtk
docker_image_version=2.1.1--pyhdfd78af_1

## Only change the variables above this line ##

workdir=$(pwd)
username=$(whoami)

docker container run --rm \
    -v $workdir:$workdir \
    -w $workdir \
    -u $(id -u $username):$(id -g $username) \
    -v $GTDBTK_DATA_PATH:$GTDBTK_DATA_PATH \
    -e GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH \
    $docker_image:$docker_image_version \
        gtdbtk classify_wf \
        --genome_dir ${assembly_dir} \
        --extension ${ext} \
        --prefix gtdb.${project} \
        --out_dir gtdbtk-results-classify \
        --cpus $cpus \
        --pplacer_cpus ${pplacer_threads} \
        --full_tree \
        --keep_intermediates
