export MODULEPATH=/home/gpu_ubuntu/public/modulefiles
module load compilers/icc/intel-cc
module load compilers/gcc/7.4.1

rm CMakeCache.txt
rm -r ./CMakeFiles/

CC=icc CXX=icpc cmake .

make > ./icc.log 2>&1
less ./icc.log
