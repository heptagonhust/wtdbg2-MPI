#!/bin/bash
source /usr/share/modules/init/bash

module load compilers/icc/intel-cc
module load compilers/gcc/7.4.1
module load librarys/glog
PARAMETERS="-xont -g144m -t16"
INPUT_FILE="/home/gpu_ubuntu/public/data/SRR6702603.fastq.gz"
OUTPUT_FILE=test
/home/gpu_ubuntu/public/software/intel/vtune_amplifier_2019.1.0.579888/bin64/amplxe-cl -collect hotspots \
    -target-duration-type long -knob enable-stack-collection=true \
    -result-dir /home/gpu_ubuntu/zhanglichen/profiling/wtdbg/$OUTPUT_FILE \
    -- ./build/wtdbg $PARAMETERS -i $INPUT_FILE -fo $OUTPUT_FILE |& tee ./log/profiling.log
