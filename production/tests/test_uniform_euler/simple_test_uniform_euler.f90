program simple_test_uniform_euler
include 'simple_lib.f08'
use simple_sym
use simple_ori
implicit none
integer, parameter  :: N_SAMPLES = 2000
type(sym)   :: pgrpsyms
type(ori)   :: o
integer     :: i
call pgrpsyms%new('c1')
call o%new(.true.)    
call pgrpsyms%rnd_euler(o)
do i = 1, N_SAMPLES
    call pgrpsyms%rnd_euler(o)
    print "(f20.15, f20.15, f20.15)", o%get_normal()*100
enddo
call pgrpsyms%kill
call o%kill
end program simple_test_uniform_euler