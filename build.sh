#!/bin/bash
source /usr/share/modules/init/bash
module load compilers/gcc/7.4.1
module load compilers/icc/intel-cc
module load librarys/glog

rm -rf build log
mkdir build
mkdir log
cd build

CC=icc CXX=icpc cmake ../src
make |& tee ../log/build.log
