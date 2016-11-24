@echo off
rem -------------------------------------------------------------------------------
rem  [system] : GpsTools
rem  [module] : mexopts.bat
rem  [func]   : compile and link options used to build mex-files on windows 32bit
rem  [note]   : environment
rem             Matlab : R2006b (7.3)
rem             Linker : MS Visual Studio 2008 (9.0)
rem             C/C++  : MS Visual Studio 2008 (9.0)
rem             Fortran: Intel Visual Fortran 10.1
rem  [version]: $Revision: 6 $ $Date: 04/11/19 8:03 $
rem  [history]: 04/03/29  0.1  new for Intel C/C++ Compiler, Intel MKL
rem             04/10/13  0.2  MKL6.1->MKL7.0.1
rem             05/03/05  0.3  MKL7.0.1->MKL7.2
rem             08/11/20  0.4  MS VS 2008 and Intel Fortran 10.1 (gt_0.6.4)
rem -------------------------------------------------------------------------------
rem ********************************************************************
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%
set VSINSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio 9.0
set VCINSTALLDIR=%VSINSTALLDIR%\VC
set SDKDIR=C:\Program Files\Microsoft SDKs\Windows\v6.0A
set IFCDIR=C:\Program Files (x86)\Intel\Compiler\Fortran\10.1.029\IA32
set PATH=%VCINSTALLDIR%\BIN\;%VCINSTALLDIR%\PlatformSDK\bin;%VSINSTALLDIR%\Common7\IDE;%SDKDIR%\bin;%VSINSTALLDIR%\Common7\Tools;%VSINSTALLDIR%\Common7\Tools\bin;%VCINSTALLDIR%\VCPackages;%MATLAB_BIN%;%PATH%
set INCLUDE=%VCINSTALLDIR%\ATLMFC\INCLUDE;%VCINSTALLDIR%\INCLUDE;%SDKDIR%\include;%INCLUDE%
set LIB=%IFCDIR%\LIB;%VCINSTALLDIR%\ATLMFC\LIB;%VCINSTALLDIR%\LIB;%SDKDIR%\lib;%MATLAB%\extern\lib\win32;%LIB%
set MW_TARGET_ARCH=win32

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=cl
set COMPFLAGS=/c /Zp8 /GR /W3 /EHsc- /Zc:wchar_t- /DMATLAB_MEX_FILE /nologo /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0
rem set OPTIMFLAGS=/MD /O2 /Oy- /DNDEBUG
rem set DEBUGFLAGS=/MD /Zi /Fd"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set OPTIMFLAGS=/MT /O2 /Oy- /DNDEBUG
set DEBUGFLAGS=/MT /Zi /Fd"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win32\microsoft
set LINKER=link
set LINKFLAGS=/dll /export:%ENTRYPOINT% /MAP /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.x /MACHINE:X86 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
set LINKDEBUGFLAGS=/DEBUG /PDB:"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%%MEX_EXT%"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del "%OUTDIR%%MEX_NAME%.map"
set POSTLINK_CMDS1=del %LIB_NAME%.x
set POSTLINK_CMDS2=mt -outputresource:"%OUTDIR%%MEX_NAME%%MEX_EXT%";2 -manifest "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest"
set POSTLINK_CMDS3=del "%OUTDIR%%MEX_NAME%%MEX_EXT%.manifest" 
