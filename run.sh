#!/bin/bash
source /usr/share/modules/init/bash

module load compilers/icc/intel-cc
module load compilers/gcc/7.4.1
module load librarys/glog
module load mpi/intel-mpi

PARAMETERS="-xont -g144m -t16"
INPUT_FILE="/home/gpu_ubuntu/public/data/SRR6702603.fastq"
OUTPUT_FILE=test
export FI_SOCKETS_IFACE=eth0
export FI_PROVIDER=tcp
mpirun -gdb -machinefile ./hostfile -genv I_MPI_DEBUG=2 ./build/wtdbg $PARAMETERS -i $INPUT_FILE -fo $OUTPUT_FILE |& tee ./log/run.log
