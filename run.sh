#!/bin/bash
source /usr/share/modules/init/bash

module load compilers/icc/intel-cc
module load compilers/gcc/7.4.1
module load librarys/glog
module load mpi/intel-mpi

PARAMETERS="-xont -g144m -t16"
PARAMETERS="-x sq -g125m -t32"
INPUT_FILE="/home/gpu_ubuntu/public/data/SRR6702603.fastq"
INPUT_FILE="/home/gpu_ubuntu/public/data/Arabidopsis_assembly.fasta"
OUTPUT_FILE=test
export FI_SOCKETS_IFACE=eth0
export FI_PROVIDER=tcp
mpirun -gdb -machinefile ./hostfile -env I_MPI_DEBUG=5 ./build/wtdbg $PARAMETERS -i $INPUT_FILE -fo $OUTPUT_FILE |& tee ./log/run.log
