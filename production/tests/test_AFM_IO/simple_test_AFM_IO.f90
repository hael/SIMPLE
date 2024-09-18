program AFM_File_IO
include 'simple_lib.f08'
use simple_AFM_image
use simple_syslib
use iso_fortran_env
use iso_c_binding

!AFM IO
type(AFM_image), allocatable :: stack_stack(:)
type(AFM_image), allocatable :: stack_avg(:)
type(parameters), target     :: params
real                         :: start, finish 
integer                      :: clip_len, file_num, file_iter
character(len=LONGSTRLEN), allocatable  :: file_list(:)
character(len = 255)    :: directory = '/Users/atifao/Downloads/IBW/'

params_glob => params
params_glob%pcontrast = 'black'
params_glob%lp        = 10.
params_glob%nsig      = 1.5

! call system('ls -1 -d /Users/atifao/Downloads/IBW/*.ibw > temp; wc -l temp > list.txt; cat temp >> list.txt; rm -f temp')
call cpu_time(start)
call simple_list_files(trim(directory) // '*.ibw', file_list)
allocate(stack_stack(size(file_list)))
allocate(stack_avg(size(file_list)))

do file_iter = 1, size(file_list)
    call read_ibw(stack_stack(file_iter), file_list(file_iter))
    call zero_padding(stack_stack(file_iter))
    call align_avg(stack_stack(file_iter), stack_avg(file_iter))
    if(file_list(file_iter) == '/Users/atifao/Downloads/IBW/Cob_450007.ibw') then 
        do clip_len = 1, size(stack_avg(file_iter)%img_array)
            call stack_avg(file_iter)%img_array(clip_len)%clip_inplace([1024, 900, 1])
        end do
    end if
    call pick_valid(stack_avg(file_iter))
end do 

call cpu_time(finish)
print *, finish - start 
end program AFM_File_IO