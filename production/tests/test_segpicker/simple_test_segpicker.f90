program simple_test_segpicker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters
use simple_pickseg
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: mic       = '/home/elmlundho/cache/NanoX/NP_Minyoung_dry/intgs/FoilHole_6388490_Data_6383078_31_20240207_130359_EER_intg.mrc'
type(parameters), target    :: params
type(pickseg) :: picker

params_glob => params
params_glob%pcontrast = 'black'
params_glob%lp        = 10.
params_glob%nsig      = 1.5
call picker%pick( mic )

end program simple_test_segpicker