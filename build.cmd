@echo off
setlocal
cd /d "%~dp0"
if defined VSCMD_VER goto compile
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not exist "%VSWHERE%" exit /b 1
for /f "usebackq delims=" %%I in (`"%VSWHERE%" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VCROOT=%%I"
if not defined VCROOT exit /b 1
call "%VCROOT%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
:compile
if not exist build mkdir build
pushd build
cl /nologo /O2 /std:c++20 /EHsc /openmp:llvm ..\src\implicert.cpp /Fe:implicert.exe /link Psapi.lib /STACK:268435456
if errorlevel 1 exit /b 1
cl /nologo /O2 /std:c++20 /EHsc /openmp:llvm ..\src\export_checkpoint.cpp /Fe:export_checkpoint.exe /link Psapi.lib /STACK:268435456
if errorlevel 1 exit /b 1
cl /nologo /O2 /std:c++20 /EHsc /openmp:llvm ..\src\trace.cpp /Fe:trace.exe /link Psapi.lib /STACK:268435456
if errorlevel 1 exit /b 1
cl /nologo /O2 /std:c++20 /fp:strict /EHsc /openmp:llvm ..\src\verify_k4.cpp /Fe:verify_k4.exe
if errorlevel 1 exit /b 1
popd
