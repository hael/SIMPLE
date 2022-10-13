program simple_test_euler
include 'simple_lib.f08'
use simple_oris, only: oris
use simple_sym,  only: sym
use simple_ori,  only: ori
implicit none
type(oris) :: oris_obj
type(sym)  :: pgrpsyms               !< symmetry elements object
type(ori)  :: o1, o2
real       :: eullims(3,2), diff_max, norm_diff, min_e(3)
integer    :: athres
call pgrpsyms%new('c1')
oris_obj = oris(1, is_ptcl=.false.)
call oris_obj%set_euler(1, [232.743408, 173.557648, 205.224075])
call oris_obj%get_ori(1, o1)
eullims = pgrpsyms%get_eullims()
do athres = 1, 15
    call oris_obj%get_ori(1, o2)
    print *, '-- athres = ', athres
    call pgrpsyms%rnd_euler(o1, real(athres), o2)
    min_e(1)  = minval( [abs(o1%e1get() - o2%e1get()), abs(o1%e1get() - (o2%e1get() - 360)), abs(o1%e1get() - (o2%e1get() + 360))] )
    min_e(2)  = minval( [abs(o1%e2get() - o2%e2get()), abs(o1%e2get() - (o2%e2get() - 180)), abs(o1%e2get() - (o2%e2get() + 180))] )
    min_e(3)  = minval( [abs(o1%e3get() - o2%e3get()), abs(o1%e3get() - (o2%e3get() - 360)), abs(o1%e3get() - (o2%e3get() + 360))] )
    diff_max  = maxval(min_e)
    norm_diff = sqrt(min_e(1)**2 + min_e(2)**2 +min_e(3)**2)
    print *, ' max  diff (e1,e2,e3) = ', diff_max
    print *, ' norm diff (o1 vs o2) = ', norm_diff
enddo
print *, eullims
call o1%print_ori()
call o2%print_ori()
end program simple_test_euler