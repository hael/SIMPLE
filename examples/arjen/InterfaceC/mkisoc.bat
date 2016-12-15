@echo off
rem
rem Example belonging to "Modern Fortran in Practice" by Arjen Markus
rem
rem This work is licensed under the Creative Commons Attribution 3.0 Unported License.
rem To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
rem or send a letter to:
rem Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
rem
rem Simple batch file to build test_iso_c program
rem
rem Note: both gfortran and Intel Fortran
rem Use "intel" as the argument to compile via Intel Fortran
rem

if /%1 == /intel goto :intel

rem gfortran ...
echo Compiling with gfortran

gcc -c sqlite3_iso_c.c
gfortran -o test_iso_c test_iso_c.f90 sqlite3_iso_c.o

goto :end

rem Intel Fortran ...
:intel
echo Compiling with Intel Fortran

cl /c sqlite3_iso_c.c
ifort test_iso_c.f90 sqlite3_iso_c.obj

:end
