@echo off

set argCount=0
for %%x in (%*) do (
   set /A argCount+=1
)

set USE_OPENMP="ON"
set USE_MPI="ON"
set GENERATOR="Visual Studio 16 2019"
set USE_FFTW="OFF"
set USE_MKL="OFF"
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

md visual_studio
cd visual_studio
cmake -G %GENERATOR% -A x64 -DCMAKE_BUILD_TYPE=Release -DUSE_TESTS=ON -Dgtest_force_shared_crt=ON -DUSE_PTESTS=ON -DBENCHMARK_ENABLE_TESTING=OFF -DUSE_OMP=%USE_OPENMP% -DUSE_MPI=%USE_MPI% -DUSE_FFTW=%USE_FFTW% -DUSE_MKL=%USE_MKL% -DPYTHON_EXECUTABLE:FILEPATH=%PYTHON% ../..
cmake --build . --config Release

xcopy /y src\pyHiChi\Release\* ..\..\bin\* > nul

cd ..
