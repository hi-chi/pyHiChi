#!/bin/bash

GENERATOR="Unix Makefiles"
USE_MPI=OFF
USE_OPENMP=OFF
USE_TESTS=OFF
USE_PTESTS=OFF
USE_FFTW=OFF
USE_MKL=OFF
FFTW_DIR=
PYTHON=
COMMAND="build"
CONFIG="Release"
COMPILER_OPTIONS=""

C_COMPILER="gcc"
CXX_COMPILER="g++"
LINKER="ld"

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

detect_compiler() {
    if command_exists icx && command_exists icpx; then
        C_COMPILER="icx"
        CXX_COMPILER="icpx"
        LINKER="icpx"
        echo "Using IntelLLVM Compiler"
    elif command_exists icc && command_exists icpc; then
        C_COMPILER="icc"
        CXX_COMPILER="icpc"
        LINKER="icpc"
        echo "Using Intel Compiler Classic"
    else
        echo "Using GCC Compiler"
    fi
    echo "Compiler version:"
    $CXX_COMPILER --version
}

build_fftw_project() {
    BUILD_DIR=$1
    if [[ ! -d $BUILD_DIR ]]; then
        mkdir -p $BUILD_DIR
    fi
    cd $BUILD_DIR
    
    FFTW_INSTALL_DIR=$(pwd)/../../3rdparty/fftw
    
    if [[ ! -d 3rdparty ]]; then
        mkdir -p 3rdparty
    fi
    cd ./3rdparty
    if [[ ! -d fftw ]]; then
        mkdir -p fftw
    fi
    cd ./fftw
    
    CMAKE_OPTIONS=(
        -G "$GENERATOR"
        -DUSE_OMP=$USE_OPENMP
        -DCMAKE_BUILD_TYPE=$CONFIG
        -DCMAKE_C_FLAGS="$COMPILER_OPTIONS"
        -DCMAKE_CXX_FLAGS="$COMPILER_OPTIONS"
        -DCMAKE_EXE_LINKER_FLAGS="$COMPILER_OPTIONS"
    )
    
    if ! CXX="$CXX_COMPILER" CC="$C_COMPILER" LD="$LINKER" \
        cmake "${CMAKE_OPTIONS[@]}" $FFTW_INSTALL_DIR; then
        echo "ERROR: CMake configuration of FFTW failed"
        exit 1
    fi
    
    if ! cmake --build ./ --config=$CONFIG -j; then
        echo "ERROR: FFTW build is failed"
        exit 1
    fi
    
    FFTW_DIR=$FFTW_INSTALL_DIR/fftw-install
    
    cd ../../..
}

build_project() {
    BUILD_DIR=$1
    if [[ ! -d $BUILD_DIR ]]; then
        mkdir -p $BUILD_DIR
    fi
    cd $BUILD_DIR
    
    CMAKE_OPTIONS=(
        -G "$GENERATOR"
        -DUSE_MPI=$USE_MPI
        -DUSE_OMP=$USE_OPENMP
        -DUSE_TESTS=$USE_TESTS
        -DUSE_PTESTS=$USE_PTESTS
        -DPython_EXECUTABLE:FILEPATH=$PYTHON
        -DCMAKE_BUILD_TYPE=$CONFIG
        -DCMAKE_C_FLAGS="$COMPILER_OPTIONS"
        -DCMAKE_CXX_FLAGS="$COMPILER_OPTIONS"
        -DCMAKE_EXE_LINKER_FLAGS="$COMPILER_OPTIONS"
    )
    
    if [[ $USE_FFTW == "ON" ]]; then
        CMAKE_OPTIONS+=(-DUSE_FFTW=$USE_FFTW)
        CMAKE_OPTIONS+=(-DFFTW_DIR=$FFTW_DIR)
    fi
    if [[ $USE_MKL == "ON" ]]; then
        CMAKE_OPTIONS+=(-DUSE_MKL=$USE_MKL)
    fi
    
    if ! CXX="$CXX_COMPILER" CC="$C_COMPILER" LD="$LINKER" \
        cmake "${CMAKE_OPTIONS[@]}" ../..; then
        echo "ERROR: CMake configuration failed"
        exit 1
    fi
    
    if ! cmake --build ./ --config=$CONFIG -j; then
        echo "ERROR: build is failed"
        exit 1
    fi
    
    mkdir -p ../../bin
    cp src/pyHiChi/pyHiChi.*.so ../../bin/pyHiChi.so
    
    cd ..
}

clean_project() {
    if [[ -d "unix_makefiles" ]]; then
        echo "Removing unix_makefiles"
        rm -rf unix_makefiles
    fi
}

show_help() {
    cat << EOF
========== Help ==========

Usage: $0 [OPTIONS] [COMMAND]

Commands:
    build       Build the project (default)
    build_fftw  Build the FFTW project only
    clean       Remove build directory

Options:
    -g "value", --generator="value", -G "value"            Set generator (default: Unix Makefiles) 
    -c "value", --config="value"                           Set build configuration (default: Release) 
    -f "value", --flags="value"                            Set additional compiler flags
    -m, --mpi, --with-mpi, -DUSE_MPI=ON                    Enable MPI support (default: OFF)
    --without-mpi, -DUSE_MPI=OFF                           Disable MPI support
    -o, --openmp, --with-openmp, -DUSE_OPENMP=ON           Enable OpenMP support (default: OFF)
    --without-openmp, -DUSE_OPENMP=OFF                     Disable OpenMP support
    --tests, --with-tests, -DUSE_TESTS=ON                  Enable Unit Tests (default: OFF)
    --without-tests, -DUSE_TESTS=OFF                       Disable Unit Tests
    --ptests, --with-ptests, -DUSE_PTESTS=ON               Enable Performance Tests (default: OFF)
    --without-ptests, -DUSE_PTESTS=OFF                     Disable Performance Tests
    --fftw, --with-fftw, -DUSE_FFTW=ON                     Enable FFTW (default: OFF)
    --without-fftw, -DUSE_FFTW=OFF                         Disable FFTW
    --mkl-fft, --with-mkl, -DUSE_MKL=ON                    Enable Intel MKL FFT (default: OFF)
    --without-mkl, -DUSE_MKL=OFF                           Disable Intel MKL FFT
    --fftw-dir, --fftw-dir                                 Set FFTW pre-installed path
    --python                                               Set python path
    -h, --help                                             Show this help message

Examples:
    $0                                                     # Build with default options
    $0 build                                               # Build with default options
    $0 build_fftw                                          # Build FFTW only
    $0 clean                                               # Remove build directory
    $0 -m -o                                               # Build with MPI and OpenMP (1)
    $0 -mpi -openmp                                        # Build with MPI and OpenMP (2)
    $0 --mpi --openmp                                      # Build with MPI and OpenMP (3)
    $0 -DUSE_MPI=ON -DUSE_OPENMP=ON                        # Build with MPI and OpenMP (4), CMake style

Full examples, bash style:
    $0 -g "Unix Makefiles" -m -o --tests --ptests --fftw --config=Release
    $0 --generator="Unix Makefiles" --with-mpi --with-openmp --with-tests --mkl-fft --with-ptests
    $0 --generator="Unix Makefiles" --without-mpi --without-openmp --without-tests --without-ptests

Full examples, CMake style:    
    $0 -G "Unix Makefiles" -DUSE_MPI=ON -DUSE_OPENMP=ON -DUSE_TESTS=ON -DUSE_PTESTS=ON -DUSE_FFTW=ON
    $0 -G "Unix Makefiles" -DUSE_MPI=OFF -DUSE_OPENMP=OFF -DUSE_TESTS=OFF -DUSE_PTESTS=OFF -DUSE_FFTW=OFF
EOF
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            # Commands
            build)
                COMMAND="build"
                shift
                ;;
            build_fftw)
                COMMAND="build_fftw"
                shift
                ;;
            clean)
                COMMAND="clean"
                shift
                ;;
            
            # Generator
            -g|-G|-generator|--generator)
                GENERATOR="$2"
                shift 2
                ;;
            --generator=*)
                GENERATOR="${1#*=}"
                shift
                ;;
            -generator=*)
                GENERATOR="${1#*=}"
                shift
                ;;
                
            # Configuration
            -c|-config|--config)
                CONFIG="$2"
                shift 2
                ;;
            --config=*)
                CONFIG="${1#*=}"
                shift
                ;;
            -config=*)
                CONFIG="${1#*=}"
                shift
                ;;
                
            # Compiler flags
            -f|-flags|--flags)
                COMPILER_OPTIONS="$2"
                shift 2
                ;;
            --flags=*)
                COMPILER_OPTIONS="${1#*=}"
                shift
                ;;
            -flags=*)
                COMPILER_OPTIONS="${1#*=}"
                shift
                ;;
            
            # MPI
            -m|-mpi|--mpi|--with-mpi)
                USE_MPI="ON"
                shift
                ;;
            -DUSE_MPI=ON)
                USE_MPI="ON"
                shift
                ;;
            --without-mpi)
                USE_MPI="OFF"
                shift
                ;;
            -DUSE_MPI=OFF)
                USE_MPI="OFF"
                shift
                ;;
            
            # OpenMP
            -o|-openmp|--openmp|--with-openmp)
                USE_OPENMP="ON"
                shift
                ;;
            -DUSE_OPENMP=ON)
                USE_OPENMP="ON"
                shift
                ;;
            --without-openmp)
                USE_OPENMP="OFF"
                shift
                ;;
            -DUSE_OPENMP=OFF)
                USE_OPENMP="OFF"
                shift
                ;;
            
            # Tests
            -tests|--tests|--with-tests)
                USE_TESTS="ON"
                shift
                ;;
            -DUSE_TESTS=ON)
                USE_TESTS="ON"
                shift
                ;;
            --without-tests)
                USE_TESTS="OFF"
                shift
                ;;
            -DUSE_TESTS=OFF)
                USE_TESTS="OFF"
                shift
                ;;

            # Performance Tests
            -ptests|--ptests|--with-ptests)
                USE_PTESTS="ON"
                shift
                ;;
            -DUSE_PTESTS=ON)
                USE_PTESTS="ON"
                shift
                ;;
            --without-ptests)
                USE_PTESTS="OFF"
                shift
                ;;
            -DUSE_PTESTS=OFF)
                USE_PTESTS="OFF"
                shift
                ;;
                
            # FFTW
            -fftw|--fftw|--with-fftw)
                USE_FFTW="ON"
                shift
                ;;
            -DUSE_FFTW=ON)
                USE_FFTW="ON"
                shift
                ;;
            --without-fftw)
                USE_FFTW="OFF"
                shift
                ;;
            -DUSE_FFTW=OFF)
                USE_FFTW="OFF"
                shift
                ;;

            # Intel MKL FFT
            -mkl_fft|--mkl_fft|-mkl-fft|--mkl-fft|-mkl|--mkl|--with-mkl)
                USE_MKL="ON"
                shift
                ;;
            -DUSE_MKL=ON)
                USE_MKL="ON"
                shift
                ;;
            --without-mkl)
                USE_MKL="OFF"
                shift
                ;;
            -DUSE_MKL=OFF)
                USE_MKL="OFF"
                shift
                ;;

            # FFTW Path
            -fftw_dir|--fftw_dir|-fftw-dir|--fftw-dir)
                FFTW_DIR="$2"
                shift 2
                ;;
            --fftw_dir=*)
                FFTW_DIR="${1#*=}"
                shift
                ;;
            --fftw-dir=*)
                FFTW_DIR="${1#*=}"
                shift
                ;;
            -fftw_dir=*)
                FFTW_DIR="${1#*=}"
                shift
                ;;
            -fftw-dir=*)
                FFTW_DIR="${1#*=}"
                shift
                ;;

            # Python Path
            -python|--python)
                PYTHON="$2"
                shift 2
                ;;
            --python=*)
                PYTHON="${1#*=}"
                shift
                ;;
            -python=*)
                PYTHON="${1#*=}"
                shift
                ;;
            
            # Help
            -h|-help|--help)
                show_help
                exit 0
                ;;
            
            *)
                echo "ERROR: Unexpected argument: $1" >&2
                show_help
                exit 1
                ;;
        esac
    done
}

check_params() {
    if [[ $USE_MKL == "ON" ]]; then
        if [[ $USE_FFTW == "ON" ]]; then
            echo "WARNING: options USE_FFTW and USE_MKL are both enabled, Intel MKL will be used"
        fi
    fi
    if [[ $USE_FFTW == "ON" ]]; then
        echo "NOTE: FFTW will be downloaded. You can set your own FFTW path with the -fftw_dir option"
    fi
    # TODO: MPI mode with FFT
    if [[ $USE_FFTW == "ON" ]]; then
        if [[ $USE_MPI == "ON" ]]; then
            echo "ERROR: FFT is not supported in MPI mode"
            exit 1
        fi
    fi
    if [[ $USE_MKL == "ON" ]]; then
        if [[ $USE_MPI == "ON" ]]; then
            echo "ERROR: FFT is not supported in MPI mode"
            exit 1
        fi
    fi
}

show_status() {
    cat << EOF
============ Configuration ============
Generator        : $GENERATOR
Configuration    : $CONFIG
USE_MPI          : $USE_MPI
USE_OPENMP       : $USE_OPENMP
USE_TESTS        : $USE_TESTS
USE_PTESTS       : $USE_PTESTS
USE_FFTW         : $USE_FFTW
USE_MKL          : $USE_MKL
FFTW_DIR         : $FFTW_DIR
PYTHON           : $PYTHON
Additional flags : $COMPILER_OPTIONS
echo =======================================
EOF
}

main() {
    detect_compiler
    parse_args "$@"
    check_params
    show_status
    
    BUILD_DIR="unix_makefiles"
    
    case $COMMAND in
        build_fftw)
            build_fftw_project $BUILD_DIR
            ;;
        build)
            if [[ $USE_FFTW == "ON" && -z "$FFTW_DIR" ]]; then
                build_fftw_project $BUILD_DIR
            fi
            build_project $BUILD_DIR
            ;;
        clean)
            clean_project
            ;;
        *)
            echo "ERROR: unknown command -- $COMMAND"
            ;;
    esac
}

main "$@"
