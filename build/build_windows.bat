@echo off

set argCount=0
for %%x in (%*) do (
   set /A argCount+=1
)

set USE_OPENMP=ON
set GENERATOR="Visual Studio 15 2017 Win64"
if not "%argCount%"=="0" (
	set GENERATOR=%1
)

set PYTHON=python

if "%2"=="-DPYTHON" (
		set PYTHON=%3
	if "%4"=="-DUSE_OPENMP" (
		if "%5"=="OFF" (
			set USE_OPENMP=OFF
		)
	)		
)
if "%2"=="-DUSE_OPENMP" (
	if "%3"=="OFF" (
		set USE_OPENMP=OFF
	)
	if "%4"=="-DPYTHON" (
		set PYTHON=%5
	)
)

md visual_studio
cd visual_studio
cmake -G %GENERATOR% -DUSE_GTEST=ON -DUSE_OPENMP=%USE_OPENMP% -DPYTHON_EXECUTABLE:FILEPATH=%PYTHON% ../..
cmake --build . --config Release

xcopy /y src\pyHiChi\Release\* ..\..\bin\* > nul

cd ..