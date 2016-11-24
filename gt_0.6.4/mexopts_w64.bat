@echo off
rem -------------------------------------------------------------------------------
rem  [system] : GpsTools
rem  [module] : mexopts_win64.bat
rem  [func]   : compile and link options used to build mex-files on windows 64bit
rem  [note]   : environment
rem             Matlab : R2007a (7.3)
rem             Linker : MS Visual Studio 2008 (9.0)
rem             C/C++  : MS Visual Studio 2008 (9.0)
rem             Fortran: Intel Visual Fortran 10.1
rem  [version]: $Revision: 6 $ $Date: 04/11/19 8:03 $
rem  [history]: 09/04/30  0.1  new
rem -------------------------------------------------------------------------------
rem ********************************************************************
rem General parameters
rem ********************************************************************
set MATLAB=%MATLAB%
set VSINSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio 9.0
set VCINSTALLDIR=%VSINSTALLDIR%\VC
set SDKDIR=C:\Program Files\Microsoft SDKs\Windows\v6.0A
set IFCDir=C:\Program Files (x86)\Intel\Compiler\Fortran\10.1.029\em64t
set PATH=%VCINSTALLDIR%\BIN\amd64;%VCINSTALLDIR%\PlatformSDK\bin\amd64;%VSINSTALLDIR%\Common7\IDE;%SDKDIR%\bin;%VSINSTALLDIR%\Common7\Tools;%VSINSTALLDIR%\Common7\Tools\bin;%VCINSTALLDIR%\VCPackages;%MATLAB_BIN%;%PATH%
set INCLUDE=%VCINSTALLDIR%\ATLMFC\INCLUDE;%VCINSTALLDIR%\INCLUDE;%SDKDIR%\include;%INCLUDE%
set LIB=%IFCDIR%\LIB;%VCINSTALLDIR%\ATLMFC\LIB\amd64;%VCINSTALLDIR%\LIB\amd64;%SDKDIR%\lib\x64;%MATLAB%\extern\lib\win64;%LIB%
set MW_TARGET_ARCH=win64

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=cl
set COMPFLAGS=-c -Zp8 -GR -W3 -EHs -D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE -D_SECURE_SCL=0 -DMATLAB_MEX_FILE -nologo
rem set OPTIMFLAGS=/MD -O2 -Oy- -DNDEBUG
rem set DEBUGFLAGS=/MD -Zi -Fd"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set OPTIMFLAGS=/MT -O2 -Oy- -DNDEBUG
set DEBUGFLAGS=/MT -Zi -Fd"%OUTDIR%%MEX_NAME%%MEX_EXT%.pdb"
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win64\microsoft
set LINKER=link
set LINKFLAGS=/dll /export:%ENTRYPOINT% /MAP /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.x /MACHINE:AMD64 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib
set LINKOPTIMFLAGS=
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
