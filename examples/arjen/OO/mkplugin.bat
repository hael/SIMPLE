rem @echo off
rem mkplugin.bat
rem Batchfile to build the plugin example
rem (using gfortran or Intel Fortran)
rem
rem Example belonging to "Modern Fortran in Practice" by Arjen Markus
rem
rem This work is licensed under the Creative Commons Attribution 3.0 Unported License.
rem To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
rem or send a letter to:
rem Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
rem
rem

if /%1 == /intel goto :intel

rem gfortran
echo Compiling with gfortran

gfortran -g -c prng_class.f90
rem gfortran -g -o prng_uniform.dll     prng_uniform_plugin.f90     prng_class.o prng_uniform.def     -shared
rem gfortran -g -o prng_exponential.dll prng_exponential_plugin.f90 prng_class.o prng_exponential.def -shared
gfortran -g -o prng_uniform.dll     prng_uniform_plugin.f90     prng_class.o -shared
gfortran -g -o prng_exponential.dll prng_exponential_plugin.f90 prng_class.o -shared

gfortran -g -o prng_factory.dll     win_gnu_dynlib.f90 dynlib.f90 prng_factory_plugin.f90 prng_class.o -shared

gfortran -g -o prng_use_plugin.exe  prng_use_plugin.f90 prng_factory.dll

goto :end

:intel
rem
rem Note the various compile and link options. They are required!

echo Compiling with Intel Fortran

ifort /c prng_class.f90
ifort /exe:prng_uniform.dll     prng_uniform_plugin.f90     prng_class.obj /dll
ifort /exe:prng_exponential.dll prng_exponential_plugin.f90 prng_class.obj /dll

ifort /exe:prng_factory.dll     win_ifort_dynlib.f90 dynlib.f90 prng_factory_plugin.f90 prng_class.obj /dll

ifort /exe:prng_use_plugin.exe  prng_use_plugin.f90 prng_factory.lib /MD

:end
echo Done
