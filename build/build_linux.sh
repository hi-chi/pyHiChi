#!/bin/sh

C_COMPILER=gcc
CXX_COMPILER=g++
LINKER=ld

# Default options are
# -DUSE_OPENMP=OFF -DUSE_GTEST=ON

CPU_OPTIONS=
MIC_OPTIONS=

CPU_VENDOR=$(grep -m 1 'vendor_id' /proc/cpuinfo | sed 's/^.*: //')
NUM_CORES=$(grep processor /proc/cpuinfo | wc -l)
CPU_FLAGS=$(grep flags /proc/cpuinfo | uniq)

command_exists ()
{
    type $1 > /dev/null 2>&1;
}

if command_exists icc && command_exists icpc ; then
    C_COMPILER=icc
    CXX_COMPILER=icpc
    LINKER=icpc
fi

if [ "$CPU_VENDOR" = "GenuineIntel" ] && [ "$CPU_FLAGS" = *"avx"* ] && [ "$CPU_OPTIONS" != *"-DUSE_AVX="* ]; then
    CPU_OPTIONS="-DUSE_AVX=ON $CPU_OPTIONS"
fi

clean=false
python_path="python"
USE_FFTW="OFF"
USE_MKL="OFF"
USE_OMP="OFF"

script=$0

while [ $# -gt 0 ]
do
key="$1"

case $key in
    -openmp)
    USE_OMP="ON"
    shift # past argument
    ;;
    -python)
    python_path="$2"
    shift # past argument
    shift # past value
    ;;
    -fftw)
    USE_FFTW="ON"
    shift # past argument
    ;;
    -mkl_fft)
    USE_MKL="ON"
    shift # past argument
    ;;
    *)    # unknown option
    shift # past argument
    ;;
esac
done

if [ $USE_OMP = "ON" ]; then
    CPU_OPTIONS="$CPU_OPTIONS -DUSE_OMP=ON"
fi
if [ $USE_FFTW = "ON" ]; then
    CPU_OPTIONS="$CPU_OPTIONS -DUSE_FFTW=ON"
fi
if [ $USE_MKL = "ON" ]; then
    CPU_OPTIONS="$CPU_OPTIONS -DUSE_MKL=ON"
fi
CPU_OPTIONS="$CPU_OPTIONS -DPYTHON_EXECUTABLE:FILEPATH=$python_path"


BUILD_DIR="unix_makefiles"
if [ ! -d $BUILD_DIR ]; then
    mkdir -p $BUILD_DIR
fi
cd $BUILD_DIR

CXX=$CXX_COMPILER CC=$C_COMPILER LD=$LINKER cmake -G "Unix Makefiles" $CPU_OPTIONS ../..
make -j $NUM_CORES -k 2> /dev/null
if [ $? -ne 0 ]; then
    make
fi
mkdir -p ../../bin
cp src/pyHiChi/pyHiChi.so ../../bin/pyHiChi.so
cd ..

if $clean ; then
    BUILD_DIR="unix_makefiles"
    if [ -d $BUILD_DIR ]; then
        rm -rf $BUILD_DIR
    fi
fi
