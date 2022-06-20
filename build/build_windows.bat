@echo off

set argCount=0
for %%x in (%*) do (
   set /A argCount+=1
)

set USE_OPENMP="OFF"
set GENERATOR="Visual Studio 15 2017"
set TOOLSET=""
set USE_FFTW="OFF"
set FFTW_DIR=""
set USE_MKL="OFF"
set USE_TESTS="ON"
set USE_PTESTS="ON"
set PYTHON="python"

:Options
if "%1"=="/openmp" (
  set USE_OPENMP="ON"
  shift
  goto Options
)
if "%1"=="/fftw" (
  set USE_FFTW="ON"
  shift
  goto Options
)
if "%1"=="/fftw_dir" (
  set FFTW_DIR=%2
  shift
  shift
  goto Options
)
if "%1"=="/mkl_fft" (
  set USE_MKL="ON"
  shift
  goto Options
)
if "%1"=="/g" (
  set GENERATOR=%2
  shift
  shift
  goto Options
)
if "%1"=="/t" (
  set TOOLSET=%2
  shift
  shift
  goto Options
)
if "%1"=="/python" (
  set PYTHON=%2
  shift
  shift
  goto Options
)
if "%1" NEQ "" (
  echo Usage: %~n0%~x0 [/python python_path] [/g generator_cmake] [/openmp] [/fftw]
  echo Error: Unknown Option: %1
  goto :OptionsEnd
)
:OptionsEnd

set GENERATOR_TOOLSET_CMAKE_LINE=-G %GENERATOR%
if %TOOLSET% NEQ "" (
  set GENERATOR_TOOLSET_CMAKE_LINE=-G %GENERATOR% -T %TOOLSET%
)


md visual_studio
cd visual_studio


REM fftw
if %USE_FFTW% EQU "ON" (
  if %FFTW_DIR% EQU "" (
    echo Installing FFTW...
    
    set FFTW_INSTALL_DIR=%cd%\..\..\3rdparty\fftw
    
    md 3rdparty\fftw
    cd 3rdparty\fftw
    
    cmake %GENERATOR_TOOLSET_CMAKE_LINE% -A x64 -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=%USE_OPENMP% %FFTW_INSTALL_DIR%
    cmake --build . --config Release
    
    cd ..\..
    
    set FFTW_DIR=%FFTW_INSTALL_DIR%\fftw-install
  )
)


REM hi-chi
cmake %GENERATOR_TOOLSET_CMAKE_LINE% -A x64 -DCMAKE_BUILD_TYPE=Release -DUSE_TESTS=%USE_TESTS% -Dgtest_force_shared_crt=ON -DUSE_PTESTS=%USE_PTESTS% -DBENCHMARK_ENABLE_TESTING=OFF -DRUN_HAVE_STD_REGEX=0 -DUSE_OMP=%USE_OPENMP% -DUSE_FFTW=%USE_FFTW% -DFFTW_DIR=%FFTW_DIR% -DUSE_MKL=%USE_MKL% -DPYTHON_EXECUTABLE:FILEPATH=%PYTHON% ../..
cmake --build . --config Release

xcopy /y src\pyHiChi\Release\* ..\..\bin\* > nul

cd ..
