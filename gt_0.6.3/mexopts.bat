@echo off
rem -------------------------------------------------------------------------------
rem  [system] : GpsTools
rem  [module] : mexopts.bat
rem  [func]   : compile and link options used for building mex-files
rem  [argin]  : none
rem  [argout] : none
rem  [note]   : based upon matlab default mexopts
rem             environment
rem               linker  : ms visual studio.net (vc7)
rem               C/C++   : intel C/C++ compiler 8.0
rem               Fortran : intel visual fortran 8.0
rem               library : intel MKL 7.0.1
rem  [version]: $Revision: 6 $ $Date: 04/11/19 8:03 $
rem  [history]: 04/03/29  0.1  new for Intel C/C++ Compiler, Intel MKL
rem             04/10/13  0.2  MKL6.1->MKL7.0.1
rem             05/03/05  0.3  MKL7.0.1->MKL7.2
rem -------------------------------------------------------------------------------
rem ********************************************************************
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%
set MSVCDir=C:\Program Files\Microsoft Visual Studio .NET\Vc7
set MSDevDir=%MSVCDir%\..\Common7
set IntelCDir=C:\Program Files\Intel\CPP\Compiler80\Ia32
set IntelFDir=C:\Program Files\Intel\Fortran\compiler80\Ia32
set IntelMklDir=C:\Program Files\Intel\MKL72
set PATH=%MSVCDir%\BIN;%MSDevDir%\bin;%MSDevDir%\IDE;%IntelCDir%\BIN;%IntelFDir%\BIN;%PATH%
set INCLUDE=%MSVCDir%\INCLUDE;%MSVCDir%\MFC\INCLUDE;%IntelMklDir%\INCLUDE;%INCLUDE%
set LIB=%MSVCDir%\LIB;%MSVCDir%\PlatformSDK\LIB;%IntelFDir%\LIB;%IntelCDir%\LIB;%IntelMklDir%\ia32\LIB;%LIB%

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
rem set COMPILER=cl
set COMPILER=icl
set COMPFLAGS=-c -Zp8 -W2 -DMATLAB_MEX_FILE -nologo
rem set OPTIMFLAGS=-O2 -Oy- -G5
set OPTIMFLAGS=-O3 -G7
set DEBUGFLAGS=-Zi
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Library creation command
rem ********************************************************************
set PRELINK_CMDS1=lib /def:"%MATLAB%\extern\include\matlab.def" /machine:ix86 /OUT:%LIB_NAME%1.lib 
set PRELINK_CMDS2=lib /def:"%MATLAB%\extern\include\libmatlbmx.def" /machine:ix86 /OUT:%LIB_NAME%2.lib 
set PRELINK_DLLS=lib /def:"%MATLAB%\extern\include\%DLL_NAME%.def" /machine:ix86 /OUT:%DLL_NAME%.lib	

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LINKER=link
set LINKFLAGS=/dll /export:mexFunction /MAP %LIB_NAME%1.lib %LIB_NAME%2.lib /implib:%LIB_NAME%.lib
set LINKDEBUGFLAGS=/debug
set LINK_FILE=mkl_c.lib libifcore.lib
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%.dll"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del "%OUTDIR%%MEX_NAME%.map"
