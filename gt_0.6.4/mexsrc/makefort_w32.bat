@echo off
rem 
rem GpsTools: compile fortran sources for windows 32bit
rem
rem compiler: Intel Visual Fortran 10.1
rem
rem 2008/12/01  0.1  new
rem

set IFCDIR=C:\Program Files (x86)\Intel\Compiler\Fortran\10.1.029\IA32\BIN
set PATH=%IFCDIR%;%PATH%
set COMPILER=ifort
set OPTFLAGS=/c /MD /O3 /nologo

@echo on
%COMPILER% %OPTFLAGS% ceppred.f /object:ceppred_w32.obj
%COMPILER% %OPTFLAGS% gmf.f /object:gmf_w32.obj
%COMPILER% %OPTFLAGS% gpt.f /object:gpt_w32.obj
%COMPILER% %OPTFLAGS% vmf1_ht.f /object:vmf1_ht_w32.obj
