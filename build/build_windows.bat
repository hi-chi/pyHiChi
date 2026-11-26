@echo off
setlocal enabledelayedexpansion

if "%~1"=="/h" goto :ShowHelp
if "%~1"=="/help" goto :ShowHelp
if "%~1"=="/?" goto :ShowHelp

set "GENERATOR=Visual Studio 18 2026"
set "TOOLSET="
set "USE_MPI=OFF"
set "USE_OPENMP=OFF"
set "USE_TESTS=OFF"
set "USE_PTESTS=OFF"
set "USE_FFTW=OFF"
set "USE_MKL=OFF"
set "FFTW_DIR="
set "PYTHON="

:Options
if "%~1"=="" goto :OptionsEnd
if "%~1"=="/g" (
    set "GENERATOR=%~2"
    shift
    shift
    goto :Options
)
if "%~1"=="/t" (
    set "TOOLSET=%~2"
    shift
    shift
    goto :Options
)
if "%~1"=="/mpi" (
    set "USE_MPI=ON"
    shift
    goto :Options
)
if "%~1"=="/openmp" (
    set "USE_OPENMP=ON"
    shift
    goto :Options
)
if "%~1"=="/tests" (
    set "USE_TESTS=ON"
    shift
    goto :Options
)
if "%~1"=="/ptests" (
    set "USE_PTESTS=ON"
    shift
    goto :Options
)
if "%~1"=="/fftw" (
    set "USE_FFTW=ON"
    shift
    goto :Options
)
if "%~1"=="/mkl_fft" (
    set "USE_MKL=ON"
    shift
    goto :Options
)
if "%~1"=="/mkl-fft" (
    set "USE_MKL=ON"
    shift
    goto :Options
)
if "%~1"=="/fftw_dir" (
    set "FFTW_DIR=%~2"
    shift
    shift
    goto :Options
)
if "%~1"=="/fftw-dir" (
    set "FFTW_DIR=%~2"
    shift
    shift
    goto :Options
)
if "%~1"=="/python" (
    set PYTHON=%~2
    shift
    shift
    goto :Options
)
echo [WARNING] Unknown argument: %~1
shift
goto :Options
:OptionsEnd

REM ========== Check parameters ==========
if "%USE_MKL%"=="ON" (
    if "%USE_FFTW%"=="ON" (
        echo WARNING: options USE_FFTW and USE_MKL are both enabled, Intel MKL will be used
    )
)
if "%USE_FFTW%"=="ON" (
    echo NOTE: FFTW will be downloaded. You can set your own FFTW path with the /fftw_dir option
)
REM TODO: MPI mode with FFT
if "%USE_FFTW%"=="ON" (
    if "%USE_MPI%"=="ON" (
        echo ERROR: FFT is not supported in MPI mode
        goto :eof
    )
)
if "%USE_MKL%"=="ON" (
    if "%USE_MPI%"=="ON" (
        echo ERROR: FFT is not supported in MPI mode
        goto :eof
    )
)

REM ========== Сontrol output ==========
echo ============ Configuration ============
echo Generator       : %GENERATOR%
if defined TOOLSET (
    echo Toolset         : %TOOLSET%
)
echo USE_MPI         : %USE_MPI%
echo USE_OPENMP      : %USE_OPENMP%
echo USE_TESTS       : %USE_TESTS%
echo USE_PTESTS      : %USE_PTESTS%
echo USE_FFTW        : %USE_FFTW%
echo USE_MKL         : %USE_MKL%
echo FFTW_DIR        : %FFTW_DIR%
echo PYTHON          : %PYTHON%
echo =======================================
REM ==============================

set "GENERATOR_TOOLSET_CMAKE_LINE="
if defined GENERATOR (
    set GENERATOR_TOOLSET_CMAKE_LINE=-G "!GENERATOR!"
)
if defined TOOLSET (
    set GENERATOR_TOOLSET_CMAKE_LINE=!GENERATOR_TOOLSET_CMAKE_LINE! -T "!TOOLSET!"
)

md visual_studio
cd visual_studio

REM ========== Download FFTW ==========
if "%USE_FFTW%"=="ON" (
    if "%FFTW_DIR%"=="" (
        echo Installing FFTW...
        
        set FFTW_INSTALL_DIR=%cd%\..\..\3rdparty\fftw
        
        md 3rdparty\fftw
        cd 3rdparty\fftw
        
        cmake !GENERATOR_TOOLSET_CMAKE_LINE! -A x64 -DUSE_OMP=%USE_OPENMP% !FFTW_INSTALL_DIR!
        cmake --build . --config Release -j
        
        cd ..\..
        
        set "FFTW_DIR=!FFTW_INSTALL_DIR!\fftw-install"
    )
)

REM ========== Build Hi-Chi ==========
set "PYTHON_CMAKE_LINE="
if defined PYTHON (
    set "PYTHON_CMAKE_LINE=-DPython_EXECUTABLE:FILEPATH=!PYTHON!"
)

set "FFT_CMAKE_LINE="
if "%USE_MKL%"=="ON" (
    set "FFT_CMAKE_LINE=!FFT_CMAKE_LINE! -DUSE_MKL=%USE_MKL%"
)
if "%USE_FFTW%"=="ON" (
    set "FFT_CMAKE_LINE=!FFT_CMAKE_LINE! -DUSE_FFTW=%USE_FFTW% -DFFTW_DIR=!FFTW_DIR!"
)

set "TEST_CMAKE_LINE="
if "%USE_TESTS%"=="ON" (
    set "TEST_CMAKE_LINE=!TEST_CMAKE_LINE! -DUSE_TESTS=%USE_TESTS% -Dgtest_force_shared_crt=ON"
)
if "%USE_PTESTS%"=="ON" (
    set "TEST_CMAKE_LINE=!TEST_CMAKE_LINE! -DUSE_PTESTS=%USE_PTESTS% -DBENCHMARK_ENABLE_TESTING=OFF"
)

cmake !GENERATOR_TOOLSET_CMAKE_LINE! -A x64 -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=%USE_MPI% -DUSE_OMP=%USE_OPENMP% !TEST_CMAKE_LINE! !FFT_CMAKE_LINE! !PYTHON_CMAKE_LINE! ../..

cmake --build . --config Release -j

xcopy /y src\pyHiChi\Release\* ..\..\bin\* > nul

cd ..
goto :eof
REM ==============================

REM ========== Help output ==========
:ShowHelp
echo.
echo ========== Help ==========
echo.
echo Usage: %~nx0 [OPTIONS]
echo.
echo OPTIONS:
echo   /g "generator"         Set generator (default: Visual Studio 18 2026)
echo   /t "toolset"           Set toolset (default: empty)
echo   /mpi                   Turn on MPI
echo   /openmp                Turn on OpenMP
echo   /tests                 Turn on unit tests
echo   /ptests                Turn on performance tests
echo   /fftw                  Include the FFTW library to perform FFT
echo   /mkl_fft, /mkl-fft     Include the Intel MKL library to perform FFT
echo   /fftw_dir, /fftw-dir   Set the FFTW pre-installed directory
echo   /h, /help              Show help
echo.
echo Examples:
echo   %~nx0
echo   %~nx0 /g "Visual Studio 18 2026"
echo   %~nx0 /g "Visual Studio 18 2026" /t "Intel C++ Compiler 2026" /mpi /openmp
echo.
echo ========================================
goto :eof
REM ==============================
