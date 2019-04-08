#!/bin/bash
module load compilers/icc/intel-cc
module load compilers/gcc/7.4.1-low
module load librarys/glog
PARAMETERS="-xont -g144m -t16"
INPUT_FILE="/home/gpu_ubuntu/public/data/SRR6702603.fastq.gz"
OUTPUT_FILE=test
./build/wtdbg $PARAMETERS -i $INPUT_FILE -fo $OUTPUT_FILE |& tee ./log/run.log
