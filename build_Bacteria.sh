#!/bin/bash

set -ie

#if [[ -z $1 ]]; then echo -e "USAGE\n\t$0 BUILD_NAME" && exit 256; fi
#BUILD_NAME=$1


#module load cmake
#module load cuda/9.0
#module load gcc/6.3.0
#module load extra
#module load openmpi/2.0.1-slurm-16.05.4
#rm -rf $BUILD_NAME
#mkdir -p $BUILD_NAME && cd $BUILD_NAME

cmake \
-DCMAKE_C_COMPILER=$(which gcc)  \
-DCMAKE_CXX_COMPILER=$(which g++)  \
-DCUDA_INCLUDE_DIRS=${CUDA_HOME}/include \
-DBUILD_TESTING=OFF \
ENABLE_TESTS_COMPILATION=False \
./

make -j 20
