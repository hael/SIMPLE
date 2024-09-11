program AFM_File_IO
include 'simple_lib.f08'
use simple_AFM_image

!AFM IO
type(AFM_image) ::  Cob16, Cob16_avg
type(parameters), target    :: params
real                        :: start, finish 

params_glob => params
params_glob%pcontrast = 'black'
params_glob%lp        = 10.
params_glob%nsig      = 1.5

call cpu_time(start)
call read_ibw1(Cob16, '/Users/atifao/Downloads/IBW/Cob_450016.ibw')
call zero_padding1(Cob16)
call align_avg1(Cob16, Cob16_avg)
call pick_valid1(Cob16_avg)
call cpu_time(finish)
print *, finish - start 
end program AFM_File_IO