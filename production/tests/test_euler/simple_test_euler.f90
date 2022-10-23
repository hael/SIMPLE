program simple_test_euler
include 'simple_lib.f08'
implicit none
type(oris) :: oris_obj
type(sym)  :: pgrpsyms               !< symmetry elements object
type(ori)  :: o1, o2
real       :: eullims(3,2), diff_max, norm_diff, min_e(3), euls(3)
integer    :: athres, isample
integer, parameter :: NSAMPLE = 100
call pgrpsyms%new('c4')
oris_obj = oris(1, is_ptcl=.false.)
call oris_obj%set_euler(1, [45., 90., 180.])
call oris_obj%get_ori(1, o1)
eullims = pgrpsyms%get_eullims()
print *, 'lim1: ', eullims(1,1), eullims(1,2)
print *, 'lim2: ', eullims(2,1), eullims(2,2)
print *, 'lim3: ', eullims(3,1), eullims(3,2)
print *, 'RND EULS'
call o2%new(.true.)
do isample = 1,NSAMPLE
    call pgrpsyms%rnd_euler(o2)
    euls = o2%get_euler()
    print *, euls(1), euls(2), euls(3)
end do




! do athres = 1, 15
!     call oris_obj%get_ori(1, o2)
!     print *, '-- athres = ', athres
!     call pgrpsyms%rnd_euler(o1, real(athres), o2)
!     min_e(1)  = minval( [abs(o1%e1get() - o2%e1get()), abs(o1%e1get() - (o2%e1get() - 360)), abs(o1%e1get() - (o2%e1get() + 360))] )
!     min_e(2)  = minval( [abs(o1%e2get() - o2%e2get()), abs(o1%e2get() - (o2%e2get() - 180)), abs(o1%e2get() - (o2%e2get() + 180))] )
!     min_e(3)  = minval( [abs(o1%e3get() - o2%e3get()), abs(o1%e3get() - (o2%e3get() - 360)), abs(o1%e3get() - (o2%e3get() + 360))] )
!     diff_max  = maxval(min_e)
!     norm_diff = sqrt(min_e(1)**2 + min_e(2)**2 +min_e(3)**2)
!     print *, ' max  diff (e1,e2,e3) = ', diff_max
!     print *, ' norm diff (o1 vs o2) = ', norm_diff
! enddo
! print *, eullims
! call o1%print_ori()
! call o2%print_ori()
end program simple_test_euler
