program simple_test_euler
include 'simple_lib.f08'
implicit none
type(oris) :: oris_obj
type(sym)  :: pgrpsyms               !< symmetry elements object
type(ori)  :: o1, o2
real       :: eullims(3,2), euls(3)
integer    :: isample
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
end program simple_test_euler
